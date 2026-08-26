### Core MIL (multiple-instance-learning) tile-level classification code: 3 MIL
### architectures (ABMIL, CLAM-SB, DSMIL) x 3 hyperparameter variants each = 9
### model/variant combinations, trained bag-by-bag (batch size 1 slide) on tile-level
### foundation-model features. Simple by design (no padded-batch DataLoader, no
### parallelization/sharding) — one slide at a time, so it runs unmodified on CPU or GPU.
###
### Takes a train CSV, a test CSV (columns: slide,label), and a path to a tile-features
### directory (one pickle per slide: {"tiles": (n_tiles,fdim) float16, "locs": (n_tiles,2),
### ...}; see github_repo/data/gen_rms_tiles.py). By default bags are restricted to tumor
### tiles only, using a precomputed indicator (see tile_tum_detector.py and
### github_repo/data/gen_rms_tumor_indicator.py -> data/RMS_tumor_indicator.p), mirroring
### RMSsubtype_tile's tumor-tile-only training. Model/variant hyperparameters are
### module-level defaults (mil_variant_specs), overridable via an optional params dict.

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
# Default hyperparameters (MIL tile-level classification)
###############################################################################
mil_models=['abmil','clam_sb','dsmil']
mil_variants=['v1','v2','v3']
max_tiles=1000
early_stop_patience=5
random_state=0

### tumor-tile filtering: restrict bags to tumor tiles only (see tile_tum_detector.py
### and data/gen_rms_tumor_indicator.py), mirroring RMSsubtype_tile's tumor-tile-only
### training. slides with fewer than min_tumor_tiles kept tiles are dropped.
use_tumor_filter=True
tum_indicator_path=os.path.join(REPO_DIR,'data','RMS_tumor_indicator.p')
min_tumor_tiles=50

### one entry per (model,variant) combo — 3 models x 3 variants = 9 combinations.
### edit any of these, or pass an override through params['mil_variant_specs'].
mil_variant_specs={
    ('abmil','v1'):{'hidden_dim':512,'attn_dim':256,'dropout':0.25,'lr':2e-4,'epochs':50,'weight_decay':1e-5},
    ('abmil','v2'):{'hidden_dim':768,'attn_dim':384,'dropout':0.15,'lr':1e-4,'epochs':50,'weight_decay':1e-5},
    ('abmil','v3'):{'hidden_dim':256,'attn_dim':128,'dropout':0.35,'lr':5e-4,'epochs':40,'weight_decay':1e-4},
    ('clam_sb','v1'):{'hidden_dim':512,'attn_dim':256,'dropout':0.25,'lr':1e-4,'epochs':20,'weight_decay':1e-5},
    ('clam_sb','v2'):{'hidden_dim':512,'attn_dim':384,'dropout':0.25,'lr':1e-4,'epochs':30,'weight_decay':1e-5},
    ('clam_sb','v3'):{'hidden_dim':384,'attn_dim':128,'dropout':0.30,'lr':2e-4,'epochs':25,'weight_decay':1e-5},
    ('dsmil','v1'):{'hidden_dim':512,'attn_dim':256,'dropout':0.25,'lr':2e-4,'epochs':50,'weight_decay':1e-5},
    ('dsmil','v2'):{'hidden_dim':768,'attn_dim':384,'dropout':0.20,'lr':1e-4,'epochs':50,'weight_decay':1e-5},
    ('dsmil','v3'):{'hidden_dim':384,'attn_dim':192,'dropout':0.30,'lr':3e-4,'epochs':35,'weight_decay':1e-4},
}


###############################################################################
# MIL model architectures (single unbatched bag: features is (n_tiles,in_dim))
###############################################################################

class ABMILClassifier(nn.Module):
    """Gated attention MIL (Ilse et al.) + linear head."""
    def __init__(self,in_dim,n_classes,hidden_dim,attn_dim,dropout):
        super().__init__()
        self.feature_proj=nn.Sequential(nn.Linear(in_dim,hidden_dim),nn.GELU(),nn.Dropout(dropout))
        self.attention_V=nn.Sequential(nn.Linear(hidden_dim,attn_dim),nn.Tanh())
        self.attention_U=nn.Sequential(nn.Linear(hidden_dim,attn_dim),nn.Sigmoid())
        self.attention_w=nn.Linear(attn_dim,1)
        self.head=nn.Linear(hidden_dim,n_classes)

    def forward(self,features):
        h=self.feature_proj(features)
        a_v=self.attention_V(h)
        a_u=self.attention_U(h)
        scores=self.attention_w(a_v*a_u).squeeze(-1)
        attn=torch.softmax(scores,dim=0)
        slide_rep=torch.sum(h*attn.unsqueeze(-1),dim=0)
        logits=self.head(slide_rep)
        return logits


class CLAMSB(nn.Module):
    """Single-branch CLAM (Lu et al.) — bag-level cross-entropy only."""
    def __init__(self,in_dim,n_classes,hidden_dim,attn_dim,dropout):
        super().__init__()
        self.fc=nn.Sequential(nn.Linear(in_dim,hidden_dim),nn.ReLU(),nn.Dropout(dropout))
        self.attention_a=nn.Sequential(nn.Linear(hidden_dim,attn_dim),nn.Tanh(),nn.Dropout(dropout))
        self.attention_b=nn.Sequential(nn.Linear(hidden_dim,attn_dim),nn.Sigmoid(),nn.Dropout(dropout))
        self.attention_c=nn.Linear(attn_dim,1)
        self.classifier=nn.Linear(hidden_dim,n_classes)

    def forward(self,features):
        h=self.fc(features)
        a=self.attention_a(h)
        b=self.attention_b(h)
        scores=self.attention_c(a*b).squeeze(-1)
        attn=torch.softmax(scores,dim=0)
        slide_rep=torch.sum(h*attn.unsqueeze(-1),dim=0)
        logits=self.classifier(slide_rep)
        return logits


class DSMIL(nn.Module):
    """Dual-stream MIL (Li et al. 2021) — critical-instance max + attention bag, averaged."""
    def __init__(self,in_dim,n_classes,hidden_dim,attn_dim,dropout):
        super().__init__()
        self.i_classifier=nn.Sequential(nn.Linear(in_dim,hidden_dim),nn.ReLU(),nn.Dropout(dropout),nn.Linear(hidden_dim,n_classes))
        self.feature_proj=nn.Sequential(nn.Linear(in_dim,hidden_dim),nn.ReLU(),nn.Dropout(dropout))
        self.attention_V=nn.Linear(hidden_dim,attn_dim)
        self.attention_U=nn.Linear(hidden_dim,attn_dim)
        self.attention_w=nn.Linear(attn_dim,1)
        self.b_classifier=nn.Linear(hidden_dim,n_classes)

    def forward(self,features):
        inst_logits=self.i_classifier(features)
        h=self.feature_proj(features)
        a_v=torch.tanh(self.attention_V(h))
        a_u=torch.sigmoid(self.attention_U(h))
        scores=self.attention_w(a_v*a_u).squeeze(-1)
        attn=torch.softmax(scores,dim=0)
        slide_rep=torch.sum(h*attn.unsqueeze(-1),dim=0)
        bag_logits=self.b_classifier(slide_rep)
        max_logits=torch.max(inst_logits,dim=0).values
        logits=0.5*(bag_logits+max_logits)
        return logits


def build_mil_model(mil_model,in_dim,n_classes,hp):
    if mil_model=='abmil':
        return ABMILClassifier(in_dim,n_classes,hp['hidden_dim'],hp['attn_dim'],hp['dropout'])
    if mil_model=='clam_sb':
        return CLAMSB(in_dim,n_classes,hp['hidden_dim'],hp['attn_dim'],hp['dropout'])
    if mil_model=='dsmil':
        return DSMIL(in_dim,n_classes,hp['hidden_dim'],hp['attn_dim'],hp['dropout'])
    raise ValueError('unknown MIL model: '+mil_model)


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


###load one slide's tile bag: {"tiles": (n_tiles,fdim) float16, "locs": (n_tiles,2), ...}
###-> (n_tiles,fdim) float32, or None. If tum_indicator has an entry for this slide, keep
###only tiles flagged as tumor (matched by loc).
def load_tile_bag(tiles_dir,slide,tum_indicator=None):
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
        if feats.shape[0]==0:
            return None

    return feats


###randomly subsample a bag down to max_tiles rows (fixed seed per slide for reproducibility)
def subsample_tiles(feats,max_tiles,seed):
    if feats.shape[0]<=max_tiles:
        return feats
    rng=np.random.RandomState(seed)
    idx=rng.choice(feats.shape[0],size=max_tiles,replace=False)
    idx.sort()
    return feats[idx]


###build (slide,feats,label) bags from a list of slides/labels, skipping missing tiles,
###label<0, or (after tumor-tile filtering) fewer than min_tumor_tiles kept tiles
def build_bags(tiles_dir,slides,labs,max_tiles,seed_base,tum_indicator=None,min_tumor_tiles=0):
    bags=[]
    for i,(slide,lab) in enumerate(zip(slides,labs)):
        if lab<0:
            continue
        feats=load_tile_bag(tiles_dir,slide,tum_indicator=tum_indicator)
        if feats is None or feats.shape[0]<min_tumor_tiles:
            continue
        feats=subsample_tiles(feats,max_tiles,seed_base+i)
        bags.append((slide,feats,int(lab)))
    return bags


###############################################################################
# Train / eval
###############################################################################

def resolve_device(device):
    if device in (None,'auto'):
        return torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    return torch.device(device)


###multi-class metrics: OVR macro AUC, accuracy, balanced accuracy, macro F1, confusion matrix
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


###train one MIL model bag-by-bag (batch size 1 slide/step), CPU or GPU
def train_mil(model,train_bags,hp,device):
    optimizer=torch.optim.Adam(model.parameters(),lr=hp['lr'],weight_decay=hp['weight_decay'])
    criterion=nn.CrossEntropyLoss()

    rng=np.random.RandomState(random_state)
    best_loss=float('inf')
    patience=0
    epochs_run=0

    for _ep in range(hp['epochs']):
        order=np.arange(len(train_bags))
        rng.shuffle(order)

        model.train()
        total_loss=0.0
        for idx in order:
            _slide,feats,lab=train_bags[idx]
            x=torch.from_numpy(feats).to(device=device,dtype=torch.float32)
            y=torch.tensor([lab],device=device,dtype=torch.long)

            optimizer.zero_grad()
            logits=model(x)
            loss=criterion(logits.unsqueeze(0),y)
            loss.backward()
            optimizer.step()
            total_loss+=float(loss.item())

        epochs_run=epochs_run+1
        avg_loss=total_loss/max(len(train_bags),1)
        if avg_loss<best_loss-1e-4:
            best_loss=avg_loss
            patience=0
        else:
            patience=patience+1
            if patience>=early_stop_patience:
                break

    return epochs_run


###evaluate a trained MIL model bag-by-bag, return classify_metrics dict
def eval_mil(model,test_bags,n_classes,device):
    model.eval()
    y_true=[]
    y_pred=[]
    probs=[]

    with torch.no_grad():
        for _slide,feats,lab in test_bags:
            x=torch.from_numpy(feats).to(device=device,dtype=torch.float32)
            logits=model(x)
            p=torch.softmax(logits,dim=0).cpu().numpy()
            probs.append(p)
            y_pred.append(int(np.argmax(p)))
            y_true.append(lab)

    probs=np.stack(probs,axis=0)
    return classify_metrics(np.array(y_true),np.array(y_pred),probs,n_classes)


###train+eval one (mil_model,mil_variant) combo; on failure return nan metrics + reason
###instead of raising, so one bad combo does not stop the whole run
def fit_and_score(mil_model,hp,train_bags,test_bags,in_dim,n_classes,device):
    try:
        torch.manual_seed(random_state)
        model=build_mil_model(mil_model,in_dim,n_classes,hp).to(device)
        epochs_run=train_mil(model,train_bags,hp,device)
        M=eval_mil(model,test_bags,n_classes,device)
        M['ok']=True
        M['reason']=None
        M['epochs_run']=epochs_run
    except Exception as exc:
        M={}
        M['ovr_auc']=float('nan')
        M['accuracy']=float('nan')
        M['balanced_accuracy']=float('nan')
        M['f1_macro']=float('nan')
        M['confusion_matrix']=None
        M['ok']=False
        M['reason']=str(exc)
        M['epochs_run']=0

    return M


###merge user-supplied params dict with the module defaults (only overrides matching keys)
def resolve_params(params):
    P={}
    P['mil_models']=mil_models
    P['mil_variants']=mil_variants
    P['mil_variant_specs']=mil_variant_specs
    P['max_tiles']=max_tiles
    P['device']='auto'
    P['use_tumor_filter']=use_tumor_filter
    P['tum_indicator_path']=tum_indicator_path
    P['min_tumor_tiles']=min_tumor_tiles

    if params is not None:
        for k in params.keys():
            P[k]=params[k]

    return P


###main entry point: run all 3x3 MIL model/variant combos given a train csv + a test csv
def classify_mil(train_csv,test_csv,tiles_dir,params=None):
    P=resolve_params(params)
    device=resolve_device(P['device'])

    tum_indicator=load_tum_indicator(P['tum_indicator_path']) if P['use_tumor_filter'] else {}

    train_slides,train_labs=load_split_csv(train_csv)
    test_slides,test_labs=load_split_csv(test_csv)

    train_bags=build_bags(tiles_dir,train_slides,train_labs,P['max_tiles'],random_state,
                           tum_indicator=tum_indicator,min_tumor_tiles=P['min_tumor_tiles'])
    test_bags=build_bags(tiles_dir,test_slides,test_labs,P['max_tiles'],random_state+10_000,
                          tum_indicator=tum_indicator,min_tumor_tiles=P['min_tumor_tiles'])

    results=[]

    train_labs_kept=[lab for _s,_f,lab in train_bags]
    if len(train_bags)<2 or len(test_bags)<1 or len(set(train_labs_kept))<2:
        return results

    in_dim=train_bags[0][1].shape[1]
    n_classes=int(max(max(train_labs_kept),max(lab for _s,_f,lab in test_bags))+1)

    for mil_model in P['mil_models']:
        for mil_variant in P['mil_variants']:
            hp=P['mil_variant_specs'][(mil_model,mil_variant)]
            M=fit_and_score(mil_model,hp,train_bags,test_bags,in_dim,n_classes,device)
            R={'mil_model':mil_model,'mil_variant':mil_variant}
            R.update(M)
            results.append(R)

    return results
