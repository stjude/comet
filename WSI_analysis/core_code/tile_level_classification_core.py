### Core TILE-LEVEL classification code: a per-tile classifier (logistic regression or a
### 2-layer MLP, both torch, CPU/GPU) trained on individual tumor tiles — every kept tile
### inherits its slide's label — with slide-level prediction = mean tile-softmax (argmax).
### Mirrors NB_MYCN_tile/rmst_tile_models.py + rmst_tile_common.py (same production design
### as RMSsubtype_tile), but self-contained and generic (binary or multi-class).
###
### Takes a train CSV, a test CSV (columns: slide,label), and a path to a tile-features
### directory (one pickle per slide: {"tiles": (n_tiles,fdim) float16, "locs": (n_tiles,2),
### ...}; see github_repo/data/gen_rms_tiles.py / gen_nb_tiles.py). By default bags are
### restricted to tumor tiles only, using a precomputed indicator (see
### core_code/tile_tum_detector.py and github_repo/data/gen_rms_tumor_indicator.py /
### gen_nb_tumor_indicator.py). Hyperparameters are module-level defaults, overridable via
### an optional params dict.

import os
import pickle

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.metrics import roc_auc_score, accuracy_score, balanced_accuracy_score, f1_score, confusion_matrix

### cap torch intra-op threads: on shared many-core nodes, unbounded threading across
### concurrently running fits causes contention/thrashing instead of speedup
torch.set_num_threads(int(os.environ.get('OMP_NUM_THREADS','8')))

SCRIPT_DIR=os.path.dirname(os.path.abspath(__file__))
REPO_DIR=os.path.dirname(SCRIPT_DIR)

###############################################################################
# Default hyperparameters (tile-level classification)
###############################################################################
tile_classifiers=['logreg','mlp2']
weight_decays=[1e-5,1e-3]
mlp_hidden_dims=[128,256]
mlp_dropouts=[0.0,0.3]
standardize=True
n_epochs=20
batch_size=8192
learning_rate=1e-3
class_weight_mode='balanced'  # 'balanced' or None
random_state=0

### tumor-tile filtering: restrict training/eval to tumor tiles only (see
### tile_tum_detector.py and data/gen_rms_tumor_indicator.py / gen_nb_tumor_indicator.py),
### mirroring the production tile-level pipelines. slides with fewer than
### min_tumor_tiles kept tiles are dropped.
use_tumor_filter=True
tum_indicator_path=os.path.join(REPO_DIR,'data','RMS_tumor_indicator.p')
min_tumor_tiles=50


###############################################################################
# Models (per-tile; standardization folded into forward pass)
###############################################################################

class Standardizer(nn.Module):
    def __init__(self,mean,std,enabled=True):
        super().__init__()
        self.enabled=bool(enabled)
        self.register_buffer('mean',mean)
        self.register_buffer('std',std)

    def forward(self,x):
        if not self.enabled:
            return x
        return (x-self.mean)/self.std


class TileLogReg(nn.Module):
    def __init__(self,in_dim,n_classes,standardizer):
        super().__init__()
        self.standardizer=standardizer
        self.linear=nn.Linear(in_dim,n_classes)

    def forward(self,x):
        return self.linear(self.standardizer(x))


class TileMLP2(nn.Module):
    def __init__(self,in_dim,n_classes,hidden_dim,dropout,standardizer):
        super().__init__()
        self.standardizer=standardizer
        self.net=nn.Sequential(
            nn.Linear(in_dim,hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim,n_classes),
        )

    def forward(self,x):
        return self.net(self.standardizer(x))


def build_tile_model(classifier,in_dim,n_classes,standardizer,hidden_dim=128,dropout=0.0):
    if classifier=='logreg':
        return TileLogReg(in_dim,n_classes,standardizer)
    if classifier=='mlp2':
        return TileMLP2(in_dim,n_classes,hidden_dim,dropout,standardizer)
    raise ValueError('unknown tile classifier: '+classifier)


def compute_standardizer(X,enabled):
    d=X.shape[1]
    if not enabled:
        return Standardizer(torch.zeros(d),torch.ones(d),enabled=False)
    mean=X.mean(axis=0)
    std=X.std(axis=0)
    std=np.where(std<1e-6,1.0,std)
    return Standardizer(
        torch.from_numpy(mean.astype(np.float32)),
        torch.from_numpy(std.astype(np.float32)),
        enabled=True,
    )


def class_weights_from_labels(y,n_classes):
    counts=np.bincount(y,minlength=n_classes).astype(np.float64)
    counts=np.maximum(counts,1.0)
    n=float(y.shape[0])
    w=n/(float(n_classes)*counts)
    return torch.tensor(w,dtype=torch.float32)


###############################################################################
# Data loading
###############################################################################

###load a train/test split csv (columns: slide,label) — same schema as the other cores
def load_split_csv(csv_file):
    df=pd.read_csv(csv_file)
    slides=list(df['slide'])
    labs=list(df['label'])
    return slides,labs


###load a tumor-indicator pickle: {slide_noext: {(x,y): bool_is_tumor}}, or {} if missing
def load_tum_indicator(path):
    if not path or not os.path.isfile(path):
        return {}
    return pickle.load(open(path,'rb'))


###load one slide's tiles: {"tiles": (n_tiles,fdim) float16, "locs": (n_tiles,2), ...} ->
###(n_tiles,fdim) float32, or None. If tum_indicator has an entry for this slide, keep
###only tiles flagged as tumor (matched by loc).
def load_slide_tiles(tiles_dir,slide,tum_indicator=None):
    slide_noext=slide.replace('.svs','')
    path=os.path.join(tiles_dir,slide_noext+'.p')
    if not os.path.isfile(path):
        return None
    d=pickle.load(open(path,'rb'))
    feats=np.asarray(d['tiles'],dtype=np.float32)
    if feats.ndim!=2 or feats.shape[0]==0:
        return None

    if tum_indicator and slide_noext in tum_indicator and d.get('locs') is not None:
        slide_ind=tum_indicator[slide_noext]
        locs=np.asarray(d['locs'])
        keep=np.array([bool(slide_ind.get((int(x),int(y)),False)) for x,y in locs])
        feats=feats[keep]

    return feats


###build a flat (n_tiles_total,fdim) tile matrix + labels from a list of slides/labels,
###also returning each kept slide's row range (for slide-level test aggregation)
def build_tile_matrix(tiles_dir,slides,labs,tum_indicator,min_tumor_tiles):
    feats_list=[]
    labs_list=[]
    slide_rows={}
    slide_y={}
    pos=0
    for slide,lab in zip(slides,labs):
        if lab<0:
            continue
        feats=load_slide_tiles(tiles_dir,slide,tum_indicator=tum_indicator)
        if feats is None or feats.shape[0]<min_tumor_tiles:
            continue
        n=feats.shape[0]
        feats_list.append(feats)
        labs_list.append(np.full(n,int(lab),dtype=np.int64))
        slide_noext=slide.replace('.svs','')
        slide_rows[slide_noext]=(pos,pos+n)
        slide_y[slide_noext]=int(lab)
        pos=pos+n

    if not feats_list:
        return None,None,{},{}

    X=np.concatenate(feats_list,axis=0)
    y=np.concatenate(labs_list,axis=0)
    return X,y,slide_rows,slide_y


###############################################################################
# Train / eval
###############################################################################

def resolve_device(device):
    if device in (None,'auto'):
        return torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    return torch.device(device)


###metrics: AUC (binary: prob of positive class; multi-class: OVR macro), accuracy,
###balanced accuracy, macro F1, confusion matrix
def classify_metrics(y_true,y_pred,probs,n_classes):
    labels=list(range(n_classes))
    try:
        if n_classes==2:
            cauc=roc_auc_score(y_true,probs[:,1])
        else:
            cauc=roc_auc_score(y_true,probs,multi_class='ovr',average='macro',labels=labels)
    except (ValueError,TypeError):
        cauc=float('nan')
    acc=accuracy_score(y_true,y_pred)
    bal_acc=balanced_accuracy_score(y_true,y_pred)
    f1m=f1_score(y_true,y_pred,average='macro',zero_division=0)
    cm=confusion_matrix(y_true,y_pred,labels=labels)

    M={}
    M['ovr_auc']=float(cauc)
    M['accuracy']=float(acc)
    M['balanced_accuracy']=float(bal_acc)
    M['f1_macro']=float(f1m)
    M['confusion_matrix']=cm
    return M


###train one per-tile model (Adam + CrossEntropyLoss, optional class weights), CPU or GPU
def train_tile_model(model,Xtrain,ytrain,device,hp):
    torch.manual_seed(random_state)
    model=model.to(device)
    model.train()
    opt=torch.optim.Adam(model.parameters(),lr=hp['learning_rate'],weight_decay=hp['weight_decay'])

    class_weight=None
    if hp['class_weight_mode']=='balanced':
        n_classes=int(ytrain.max())+1
        class_weight=class_weights_from_labels(ytrain,n_classes).to(device)
    loss_fn=nn.CrossEntropyLoss(weight=class_weight)

    Xtr=torch.from_numpy(Xtrain)
    ytr=torch.from_numpy(ytrain)
    n=Xtr.shape[0]
    bs=min(hp['batch_size'],max(1,n))
    rng=np.random.RandomState(random_state)

    for _ep in range(hp['n_epochs']):
        perm=rng.permutation(n)
        for start in range(0,n,bs):
            sel=perm[start:start+bs]
            xb=Xtr[sel].to(device)
            yb=ytr[sel].to(device)
            opt.zero_grad()
            loss=loss_fn(model(xb),yb)
            loss.backward()
            opt.step()

    return model


###evaluate a trained model on test tiles, aggregating to slide level (mean softmax)
@torch.no_grad()
def eval_tile_model(model,Xtest,test_slide_rows,test_slide_y,device,n_classes,batch_size):
    model.eval()
    Xte=torch.from_numpy(Xtest)
    n=Xte.shape[0]
    probs_all=np.empty((n,n_classes),dtype=np.float32)
    for start in range(0,n,batch_size):
        xb=Xte[start:start+batch_size].to(device)
        p=torch.softmax(model(xb),dim=1).cpu().numpy()
        probs_all[start:start+xb.shape[0]]=p

    y_true=[]
    slide_probs=[]
    for slide_noext,(a,b) in test_slide_rows.items():
        slide_probs.append(probs_all[a:b].mean(axis=0))
        y_true.append(test_slide_y[slide_noext])

    y_true=np.asarray(y_true,dtype=np.int64)
    slide_probs=np.asarray(slide_probs,dtype=np.float64)
    y_pred=np.argmax(slide_probs,axis=1).astype(np.int64)
    return classify_metrics(y_true,y_pred,slide_probs,n_classes)


###train+eval one (classifier,hyperparam) combo; on failure return nan metrics + reason
###instead of raising, so one bad combo does not stop the whole run
def fit_and_score(classifier,hp,Xtrain,ytrain,Xtest,test_slide_rows,test_slide_y,n_classes,device):
    try:
        standardizer=compute_standardizer(Xtrain,hp['standardize'])
        model=build_tile_model(classifier,Xtrain.shape[1],n_classes,standardizer,
                                hidden_dim=hp.get('hidden_dim',128),dropout=hp.get('dropout',0.0))
        model=train_tile_model(model,Xtrain,ytrain,device,hp)
        M=eval_tile_model(model,Xtest,test_slide_rows,test_slide_y,device,n_classes,hp['batch_size'])
        M['ok']=True
        M['reason']=None
    except Exception as exc:
        M={}
        M['ovr_auc']=float('nan')
        M['accuracy']=float('nan')
        M['balanced_accuracy']=float('nan')
        M['f1_macro']=float('nan')
        M['confusion_matrix']=None
        M['ok']=False
        M['reason']=str(exc)

    return M


###merge user-supplied params dict with the module defaults (only overrides matching keys)
def resolve_params(params):
    P={}
    P['tile_classifiers']=tile_classifiers
    P['weight_decays']=weight_decays
    P['mlp_hidden_dims']=mlp_hidden_dims
    P['mlp_dropouts']=mlp_dropouts
    P['standardize']=standardize
    P['n_epochs']=n_epochs
    P['batch_size']=batch_size
    P['learning_rate']=learning_rate
    P['class_weight_mode']=class_weight_mode
    P['device']='auto'
    P['use_tumor_filter']=use_tumor_filter
    P['tum_indicator_path']=tum_indicator_path
    P['min_tumor_tiles']=min_tumor_tiles

    if params is not None:
        for k in params.keys():
            P[k]=params[k]

    return P


###main entry point: run all (classifier,hyperparam) combos given a train csv + a test csv
def classify_tile_level(train_csv,test_csv,tiles_dir,params=None):
    P=resolve_params(params)
    device=resolve_device(P['device'])

    tum_indicator=load_tum_indicator(P['tum_indicator_path']) if P['use_tumor_filter'] else {}

    train_slides,train_labs=load_split_csv(train_csv)
    test_slides,test_labs=load_split_csv(test_csv)

    Xtrain,ytrain,_train_rows,_train_y=build_tile_matrix(
        tiles_dir,train_slides,train_labs,tum_indicator,P['min_tumor_tiles'])
    Xtest,_ytest_tiles,test_slide_rows,test_slide_y=build_tile_matrix(
        tiles_dir,test_slides,test_labs,tum_indicator,P['min_tumor_tiles'])

    results=[]
    if Xtrain is None or Xtest is None or len(test_slide_rows)<1 or len(np.unique(ytrain))<2:
        return results

    n_classes=int(max(ytrain.max(),max(test_slide_y.values()))+1)

    for classifier in P['tile_classifiers']:
        if classifier=='logreg':
            for wd in P['weight_decays']:
                hp={
                    'weight_decay':wd,'standardize':P['standardize'],'n_epochs':P['n_epochs'],
                    'batch_size':P['batch_size'],'learning_rate':P['learning_rate'],
                    'class_weight_mode':P['class_weight_mode'],
                }
                M=fit_and_score('logreg',hp,Xtrain,ytrain,Xtest,test_slide_rows,test_slide_y,n_classes,device)
                R={'classifier':'logreg','weight_decay':wd}
                R.update(M)
                results.append(R)

        elif classifier=='mlp2':
            for hdim in P['mlp_hidden_dims']:
                for do in P['mlp_dropouts']:
                    for wd in P['weight_decays']:
                        hp={
                            'weight_decay':wd,'hidden_dim':hdim,'dropout':do,
                            'standardize':P['standardize'],'n_epochs':P['n_epochs'],
                            'batch_size':P['batch_size'],'learning_rate':P['learning_rate'],
                            'class_weight_mode':P['class_weight_mode'],
                        }
                        M=fit_and_score('mlp2',hp,Xtrain,ytrain,Xtest,test_slide_rows,test_slide_y,n_classes,device)
                        R={'classifier':'mlp2','hidden_dim':hdim,'dropout':do,'weight_decay':wd}
                        R.update(M)
                        results.append(R)

    return results
