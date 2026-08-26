### Core methylation-based tumor-classification code: PCA + {LDA, Logistic Regression,
### 2-layer MLP} on top-N-variable-probe methylation representations. Same models and
### architecture as the real methylation (top-20K) classification pipeline: LDA, a
### single-layer torch LR (softmax), and a 2-layer torch MLP — no QDA/RF here, matching
### that reference (unlike the WSI core, which also has QDA and RF).
###
### Takes a train CSV, a test CSV (columns: slide,label), and a path to a methylation
### representations pickle (see github_repo/data/gen_meth_top5000_sample_reps.py):
###   {'probe_ids','sentrix_ids','beta':(N,fdim) float16,'metadata':{sentrix_id:{'comet_id','slide'}}}
### The slide -> Sentrix mapping needed to join the CSVs to the beta matrix is read from
### that pickle's own metadata, so the train/test CSVs can be the exact same ones used by
### the WSI core (same splits, same slide,label schema).
###
### Integrative (WSI+methylation) classification is separate and not implemented here.

import os
import pickle

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import roc_auc_score, accuracy_score, balanced_accuracy_score, f1_score, confusion_matrix

### cap torch intra-op threads: on shared many-core nodes, unbounded threading across
### concurrently running fits causes contention/thrashing instead of speedup
torch.set_num_threads(int(os.environ.get('OMP_NUM_THREADS','8')))

###############################################################################
# Default hyperparameters (methylation classification)
###############################################################################
PC_dims=[500,1000]
classifiers=['lda','lr','mlp']
random_state=0

### torch training (shared by 'lr' and 'mlp'; matches the WSI torch MLP hyperparameters)
torch_batch_size=64
torch_n_epochs=50
torch_reg_l1=1e-5
torch_reg_l2=1e-5
torch_learning_rate=0.05

### 'mlp' architecture only ('lr' is a single Linear(in,C) layer, no hidden layer)
mlp_hidden_dim=512
mlp_dropout=0.2


###load methylation reps pickle: {'probe_ids','sentrix_ids','beta','metadata'}
def load_methylation_reps(methylation_reps_file):
    methD=pickle.load(open(methylation_reps_file,'rb'))
    return methD


###load a train/test split csv (columns: slide,label) — same schema as the WSI core
def load_split_csv(csv_file):
    df=pd.read_csv(csv_file)
    slides=list(df['slide'])
    labs=list(df['label'])
    return slides,labs


###slide -> Sentrix id, read from the methylation pickle's own metadata
def build_slide_to_sentrix(metadata):
    slide_to_sentrix={}
    for sx,m in metadata.items():
        slide=m.get('slide')
        if slide is not None:
            slide_to_sentrix[slide]=sx
    return slide_to_sentrix


###build X,Y from the methylation dict + list of slides/labels (skip slides with no Sentrix match or label<0)
def build_XY(methD,slides,labs):
    beta=methD['beta']
    sentrix_ids=methD['sentrix_ids']
    sentrix_index={sx:i for i,sx in enumerate(sentrix_ids)}
    slide_to_sentrix=build_slide_to_sentrix(methD['metadata'])

    snum=len(slides)
    fdim=beta.shape[1]
    X=np.zeros((snum,fdim),dtype=np.float32)
    Y=np.zeros((snum,),dtype=np.int32)

    cnt=0
    for i in range(snum):
        cslide=slides[i]
        clab=labs[i]
        sx=slide_to_sentrix.get(cslide)
        if (sx is not None) and (clab>-1):
            X[cnt,:]=beta[sentrix_index[sx],:].astype(np.float32)
            Y[cnt]=clab
            cnt=cnt+1

    X=X[:cnt,:]
    Y=Y[:cnt]
    return X,Y


###train-fold-only PCA (fit on train, applied to train+test)
def pca_reduce(Xtrain,Xtest,n_components):
    ncomp=min(n_components,Xtrain.shape[0],Xtrain.shape[1])
    pca=PCA(n_components=ncomp,svd_solver='full',random_state=random_state).fit(Xtrain)
    pXtrain=pca.transform(Xtrain)
    pXtest=pca.transform(Xtest)
    return pXtrain,pXtest


###multi-class metrics: OVR macro AUC, accuracy, balanced accuracy, macro F1, confusion matrix
def classify_metrics(Ytest,probs):
    ypred=np.argmax(probs,axis=1)
    try:
        if probs.shape[1]==2:
            cauc=roc_auc_score(Ytest,probs[:,1])
        else:
            cauc=roc_auc_score(Ytest,probs,multi_class='ovr',average='macro')
    except (ValueError,TypeError):
        cauc=float('nan')
    acc=accuracy_score(Ytest,ypred)
    bal_acc=balanced_accuracy_score(Ytest,ypred)
    f1m=f1_score(Ytest,ypred,average='macro',zero_division=0)
    cm=confusion_matrix(Ytest,ypred)

    M={}
    M['ovr_auc']=float(cauc)
    M['accuracy']=float(acc)
    M['balanced_accuracy']=float(bal_acc)
    M['f1_macro']=float(f1m)
    M['confusion_matrix']=cm
    return M


###torch net: Linear(in,C) if hidden_dim is None (LR), else Linear(in,hidden)->ReLU->Dropout->Linear(hidden,C) (MLP)
def build_torch_net(n_in,n_classes,hidden_dim,dropout_p):
    if hidden_dim is None:
        net=nn.Linear(n_in,n_classes,bias=True)
    else:
        net=nn.Sequential(
            nn.Linear(n_in,hidden_dim,bias=True),
            nn.ReLU(),
            nn.Dropout(p=dropout_p),
            nn.Linear(hidden_dim,n_classes,bias=True),
        )
    return net.to(dtype=torch.float32)


###train a torch net on CPU (Adam + cross-entropy + L1/L2 on all parameters), return test-set probs
###shared by 'lr' (hidden_dim=None) and 'mlp' (hidden_dim=torch_params['hidden_dim'])
def torch_predict_proba(Xtrain,Ytrain,Xtest,Ytest,torch_params):
    device=torch.device('cpu')
    n_classes=int(max(Ytrain.max(),Ytest.max())+1)
    n_in=int(Xtrain.shape[1])

    torch.manual_seed(random_state)
    net=build_torch_net(n_in,n_classes,torch_params.get('hidden_dim'),torch_params.get('dropout'))

    Xtr=torch.from_numpy(np.ascontiguousarray(Xtrain)).to(device=device,dtype=torch.float32)
    Ytr=torch.from_numpy(np.ascontiguousarray(Ytrain.astype(np.int64))).to(device=device,dtype=torch.long)
    Xte=torch.from_numpy(np.ascontiguousarray(Xtest)).to(device=device,dtype=torch.float32)

    bs=min(torch_params['batch_size'],max(1,Xtr.shape[0]))
    loader=DataLoader(TensorDataset(Xtr,Ytr),batch_size=bs,shuffle=True,drop_last=False)
    opt=torch.optim.Adam(net.parameters(),lr=torch_params['learning_rate'])
    ce_loss=nn.CrossEntropyLoss()

    net.train()
    for _ep in range(torch_params['n_epochs']):
        for xb,yb in loader:
            opt.zero_grad()
            logits=net(xb)
            loss=ce_loss(logits,yb)
            l1=torch.zeros((),dtype=torch.float32)
            l2=torch.zeros((),dtype=torch.float32)
            for p in net.parameters():
                l1=l1+p.abs().sum()
                l2=l2+(p*p).sum()
            loss=loss+torch_params['reg_l1']*l1+torch_params['reg_l2']*l2
            loss.backward()
            opt.step()

    net.eval()
    with torch.no_grad():
        logits_te=net(Xte)
        probs=torch.softmax(logits_te,dim=1).cpu().numpy().astype(np.float64)

    return probs


###fit one classifier on (Xtrain,Ytrain), score on (Xtest,Ytest)
###on failure, return nan metrics + reason instead of raising, so one bad combo does not stop the whole run
def fit_and_score(kind,params,Xtrain,Ytrain,Xtest,Ytest):
    try:
        if kind=='lda':
            clf=LinearDiscriminantAnalysis(solver='svd')
            clf.fit(Xtrain,Ytrain)
            probs=clf.predict_proba(Xtest)
        elif kind=='lr':
            probs=torch_predict_proba(Xtrain,Ytrain,Xtest,Ytest,params)
        elif kind=='mlp':
            probs=torch_predict_proba(Xtrain,Ytrain,Xtest,Ytest,params)
        else:
            raise ValueError('unknown classifier kind: '+kind)

        M=classify_metrics(Ytest,probs)
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
    P['PC_dims']=PC_dims
    P['classifiers']=classifiers
    P['torch_batch_size']=torch_batch_size
    P['torch_n_epochs']=torch_n_epochs
    P['torch_reg_l1']=torch_reg_l1
    P['torch_reg_l2']=torch_reg_l2
    P['torch_learning_rate']=torch_learning_rate
    P['mlp_hidden_dim']=mlp_hidden_dim
    P['mlp_dropout']=mlp_dropout

    if params is not None:
        for k in params.keys():
            P[k]=params[k]

    return P


###main entry point: classify methylation reps given a train csv + a test csv
def classify_methylation(train_csv,test_csv,methylation_reps_file,params=None):
    P=resolve_params(params)

    methD=load_methylation_reps(methylation_reps_file)
    train_slides,train_labs=load_split_csv(train_csv)
    test_slides,test_labs=load_split_csv(test_csv)

    Xtrain,Ytrain=build_XY(methD,train_slides,train_labs)
    Xtest,Ytest=build_XY(methD,test_slides,test_labs)

    results=[]

    if Xtrain.shape[0]<2 or Xtest.shape[0]<1 or len(np.unique(Ytrain))<2:
        return results

    for cdim in P['PC_dims']:
        pXtrain,pXtest=pca_reduce(Xtrain,Xtest,cdim)

        if 'lda' in P['classifiers']:
            M=fit_and_score('lda',{},pXtrain,Ytrain,pXtest,Ytest)
            R={'classifier':'lda','PC_dim':cdim}
            R.update(M)
            results.append(R)

        if 'lr' in P['classifiers']:
            lr_params={
                'hidden_dim':None,
                'dropout':None,
                'batch_size':P['torch_batch_size'],
                'n_epochs':P['torch_n_epochs'],
                'reg_l1':P['torch_reg_l1'],
                'reg_l2':P['torch_reg_l2'],
                'learning_rate':P['torch_learning_rate'],
            }
            M=fit_and_score('lr',lr_params,pXtrain,Ytrain,pXtest,Ytest)
            R={'classifier':'lr','PC_dim':cdim}
            R.update(M)
            results.append(R)

        if 'mlp' in P['classifiers']:
            mlp_params={
                'hidden_dim':P['mlp_hidden_dim'],
                'dropout':P['mlp_dropout'],
                'batch_size':P['torch_batch_size'],
                'n_epochs':P['torch_n_epochs'],
                'reg_l1':P['torch_reg_l1'],
                'reg_l2':P['torch_reg_l2'],
                'learning_rate':P['torch_learning_rate'],
            }
            M=fit_and_score('mlp',mlp_params,pXtrain,Ytrain,pXtest,Ytest)
            R={'classifier':'mlp','PC_dim':cdim}
            R.update(M)
            results.append(R)

    return results
