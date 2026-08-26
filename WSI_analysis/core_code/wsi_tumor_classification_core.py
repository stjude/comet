### Core WSI tumor-classification code: PCA + {LDA, QDA, Random Forest, Logistic
### Regression, 2-layer MLP} on SAMPLER (percentile) representations. Takes a train
### CSV, a test CSV (columns: slide,label), and a path to a sampler-representations
### pickle (slide -> (10,fdim) array; see github_repo/data). WSI/sampler-only for now —
### methylation-based and integrative (WSI+methylation) classification are separate
### and not implemented here.

import os
import pickle

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, accuracy_score, balanced_accuracy_score, f1_score, confusion_matrix

### cap torch intra-op threads: on shared many-core nodes, unbounded threading across
### concurrently running fits causes contention/thrashing instead of speedup
torch.set_num_threads(int(os.environ.get('OMP_NUM_THREADS','8')))

###############################################################################
# Default hyperparameters (WSI sampler classification)
###############################################################################
PC_dims=[500,1000]
Cvec=[0.1,1,10,100]
num_trees=[50,100,200]
rf_max_depths=[2,3,5,10]
QDA_reg_vals=[0,0.25,0.5]
classifiers=['lda','qda','rf','lr','mlp']
random_state=0

### 2-layer MLP architecture + training (matches the model used for WSI tumor
### classification: Linear(in,hidden) -> ReLU -> Dropout -> Linear(hidden,C) logits,
### softmax at eval, CPU only). Single fixed setting (not swept like the other
### classifiers) — edit these to try a different architecture.
mlp_hidden_dim=512
mlp_dropout=0.2
mlp_batch_size=64
mlp_n_epochs=50
mlp_reg_l1=1e-5
mlp_reg_l2=1e-5
mlp_learning_rate=0.05


###load sampler reps pickle: {slide: (10,fdim) array}
def load_sampler_reps(sampler_reps_file):
    sampD=pickle.load(open(sampler_reps_file,'rb'))
    return sampD


###load a train/test split csv (columns: slide,label)
def load_split_csv(csv_file):
    df=pd.read_csv(csv_file)
    slides=list(df['slide'])
    labs=list(df['label'])
    return slides,labs


###build X,Y from a sampler dict + list of slides/labels (skip slides missing from sampD or label<0)
def build_XY(sampD,slides,labs):
    snum=len(slides)
    sampDslides=list(sampD.keys())
    fdim=len(np.ravel(sampD[sampDslides[0]]))

    X=np.zeros((snum,fdim),dtype=np.float32)
    Y=np.zeros((snum,),dtype=np.int32)

    cnt=0
    for i in range(snum):
        cslide=slides[i]
        clab=labs[i]
        if (cslide in sampDslides) and (clab>-1):
            X[cnt,:]=np.ravel(sampD[cslide]).astype(np.float32)
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


###metrics: AUC (binary: prob of positive class; multi-class: OVR macro), accuracy,
###balanced accuracy, macro F1, confusion matrix
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


###2-layer MLP: Linear(in,hidden) -> ReLU -> Dropout -> Linear(hidden,C) logits
def build_mlp(n_in,n_classes,hidden_dim,dropout_p):
    net=nn.Sequential(
        nn.Linear(n_in,hidden_dim,bias=True),
        nn.ReLU(),
        nn.Dropout(p=dropout_p),
        nn.Linear(hidden_dim,n_classes,bias=True),
    )
    return net.to(dtype=torch.float32)


###train the MLP on CPU (Adam + cross-entropy + L1/L2 on all parameters), return test-set probs
def mlp_predict_proba(Xtrain,Ytrain,Xtest,Ytest,mlp_params):
    device=torch.device('cpu')
    n_classes=int(max(Ytrain.max(),Ytest.max())+1)
    n_in=int(Xtrain.shape[1])

    torch.manual_seed(random_state)
    net=build_mlp(n_in,n_classes,mlp_params['hidden_dim'],mlp_params['dropout'])

    Xtr=torch.from_numpy(np.ascontiguousarray(Xtrain)).to(device=device,dtype=torch.float32)
    Ytr=torch.from_numpy(np.ascontiguousarray(Ytrain.astype(np.int64))).to(device=device,dtype=torch.long)
    Xte=torch.from_numpy(np.ascontiguousarray(Xtest)).to(device=device,dtype=torch.float32)

    bs=min(mlp_params['batch_size'],max(1,Xtr.shape[0]))
    loader=DataLoader(TensorDataset(Xtr,Ytr),batch_size=bs,shuffle=True,drop_last=False)
    opt=torch.optim.Adam(net.parameters(),lr=mlp_params['learning_rate'])
    ce_loss=nn.CrossEntropyLoss()

    net.train()
    for _ep in range(mlp_params['n_epochs']):
        for xb,yb in loader:
            opt.zero_grad()
            logits=net(xb)
            loss=ce_loss(logits,yb)
            l1=torch.zeros((),dtype=torch.float32)
            l2=torch.zeros((),dtype=torch.float32)
            for p in net.parameters():
                l1=l1+p.abs().sum()
                l2=l2+(p*p).sum()
            loss=loss+mlp_params['reg_l1']*l1+mlp_params['reg_l2']*l2
            loss.backward()
            opt.step()

    net.eval()
    with torch.no_grad():
        logits_te=net(Xte)
        probs=torch.softmax(logits_te,dim=1).cpu().numpy().astype(np.float64)

    return probs


###fit one classifier on (Xtrain,Ytrain), score on (Xtest,Ytest)
###on failure (e.g. QDA covariance not full rank for a small class), return nan metrics + reason
###instead of raising, so one bad (classifier,hyperparam) combo does not stop the whole run
def fit_and_score(kind,params,Xtrain,Ytrain,Xtest,Ytest):
    try:
        if kind=='lda':
            clf=LinearDiscriminantAnalysis(solver='svd')
            clf.fit(Xtrain,Ytrain)
            probs=clf.predict_proba(Xtest)
        elif kind=='qda':
            clf=QuadraticDiscriminantAnalysis(reg_param=params['reg_param'])
            clf.fit(Xtrain,Ytrain)
            probs=clf.predict_proba(Xtest)
        elif kind=='rf':
            clf=RandomForestClassifier(n_estimators=params['n_estimators'],max_depth=params['max_depth'],random_state=random_state,n_jobs=1)
            clf.fit(Xtrain,Ytrain)
            probs=clf.predict_proba(Xtest)
        elif kind=='lr':
            clf=LogisticRegression(C=params['C'],max_iter=1000,class_weight='balanced',random_state=random_state)
            clf.fit(Xtrain,Ytrain)
            probs=clf.predict_proba(Xtest)
        elif kind=='mlp':
            probs=mlp_predict_proba(Xtrain,Ytrain,Xtest,Ytest,params)
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
    P['Cvec']=Cvec
    P['num_trees']=num_trees
    P['rf_max_depths']=rf_max_depths
    P['QDA_reg_vals']=QDA_reg_vals
    P['classifiers']=classifiers
    P['mlp_hidden_dim']=mlp_hidden_dim
    P['mlp_dropout']=mlp_dropout
    P['mlp_batch_size']=mlp_batch_size
    P['mlp_n_epochs']=mlp_n_epochs
    P['mlp_reg_l1']=mlp_reg_l1
    P['mlp_reg_l2']=mlp_reg_l2
    P['mlp_learning_rate']=mlp_learning_rate

    if params is not None:
        for k in params.keys():
            P[k]=params[k]

    return P


###main entry point: classify sampler reps given a train csv + a test csv
def classify_wsi_sampler(train_csv,test_csv,sampler_reps_file,params=None):
    P=resolve_params(params)

    sampD=load_sampler_reps(sampler_reps_file)
    train_slides,train_labs=load_split_csv(train_csv)
    test_slides,test_labs=load_split_csv(test_csv)

    Xtrain,Ytrain=build_XY(sampD,train_slides,train_labs)
    Xtest,Ytest=build_XY(sampD,test_slides,test_labs)

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

        if 'qda' in P['classifiers']:
            for rp in P['QDA_reg_vals']:
                M=fit_and_score('qda',{'reg_param':rp},pXtrain,Ytrain,pXtest,Ytest)
                R={'classifier':'qda','PC_dim':cdim,'reg_param':rp}
                R.update(M)
                results.append(R)

        if 'lr' in P['classifiers']:
            for cval in P['Cvec']:
                M=fit_and_score('lr',{'C':cval},pXtrain,Ytrain,pXtest,Ytest)
                R={'classifier':'lr','PC_dim':cdim,'C':cval}
                R.update(M)
                results.append(R)

        if 'rf' in P['classifiers']:
            for ntree in P['num_trees']:
                for depth in P['rf_max_depths']:
                    M=fit_and_score('rf',{'n_estimators':ntree,'max_depth':depth},pXtrain,Ytrain,pXtest,Ytest)
                    R={'classifier':'rf','PC_dim':cdim,'n_estimators':ntree,'max_depth':depth}
                    R.update(M)
                    results.append(R)

        if 'mlp' in P['classifiers']:
            mlp_params={
                'hidden_dim':P['mlp_hidden_dim'],
                'dropout':P['mlp_dropout'],
                'batch_size':P['mlp_batch_size'],
                'n_epochs':P['mlp_n_epochs'],
                'reg_l1':P['mlp_reg_l1'],
                'reg_l2':P['mlp_reg_l2'],
                'learning_rate':P['mlp_learning_rate'],
            }
            M=fit_and_score('mlp',mlp_params,pXtrain,Ytrain,pXtest,Ytest)
            R={'classifier':'mlp','PC_dim':cdim}
            R.update(M)
            results.append(R)

    return results
