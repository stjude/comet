### Core integrative (WSI + methylation) tumor-classification code: PCA + {LDA, QDA,
### Random Forest, Logistic Regression, 2-layer MLP} on concat([WSI sampler feats, top-N
### methylation probes]). Same models as the WSI core; LR/MLP use the same torch training
### as the methylation core (single Linear layer for LR, Linear->ReLU->Dropout->Linear for
### MLP), matching the real integrative pipeline (separate LDA / torch LR / torch MLP /
### sklearn QDA+RF scripts, same reference).
###
### Two PCA modes (module hyperparameter 'pca_modes'), matching the real
### 'concat_pca' / 'double_pca' preprocessing options:
###   'double_pca' — PCA(WSI block, PC_dim comps) || PCA(meth block, PC_dim comps),
###                  each fit on train only, then concatenated.
###   'concat_pca' — PCA on the full raw [WSI || meth] concatenated vector, in ONE PCA fit
###                  on train. To keep the total feature count comparable to double_pca
###                  (which has 2*PC_dim total components), concat_pca uses 2*PC_dim
###                  components.
###
### Takes a train CSV, a test CSV (columns: slide,label), a path to a sampler
### representations pickle (slide -> (10,fdim) array), and a path to a methylation
### representations pickle (see github_repo/data). The WSI/methylation join is per slide:
### WSI features come straight from the sampler pickle; methylation features are resolved
### via slide -> SentrixID using the methylation pickle's own metadata (same as the
### methylation core).

import os
import pickle

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, accuracy_score, balanced_accuracy_score, f1_score, confusion_matrix

### cap torch intra-op threads: on shared many-core nodes, unbounded threading across
### concurrently running fits causes contention/thrashing instead of speedup
torch.set_num_threads(int(os.environ.get('OMP_NUM_THREADS','8')))

###############################################################################
# Default hyperparameters (integrative WSI + methylation classification)
###############################################################################
### per-modality PCA target for 'double_pca'; 'concat_pca' uses 2*PC_dim (see module docstring)
PC_dims=[500,1000]
pca_modes=['concat_pca','double_pca']
num_trees=[50,100,200]
rf_max_depths=[2,3,5,10]
QDA_reg_vals=[0,0.25,0.5]
classifiers=['lda','qda','rf','lr','mlp']
random_state=0

### torch training (shared by 'lr' and 'mlp'; matches the WSI/methylation torch hyperparameters)
torch_batch_size=64
torch_n_epochs=50
torch_reg_l1=1e-5
torch_reg_l2=1e-5
torch_learning_rate=0.05

### 'mlp' architecture only ('lr' is a single Linear(in,C) layer, no hidden layer)
mlp_hidden_dim=512
mlp_dropout=0.2


###load sampler reps pickle: {slide: (10,fdim) array}
def load_sampler_reps(sampler_reps_file):
    sampD=pickle.load(open(sampler_reps_file,'rb'))
    return sampD


###load methylation reps pickle: {'probe_ids','sentrix_ids','beta','metadata'}
def load_methylation_reps(methylation_reps_file):
    methD=pickle.load(open(methylation_reps_file,'rb'))
    return methD


###load a train/test split csv (columns: slide,label) — same schema as the WSI/methylation cores
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


###build X,Y from sampler dict + methylation dict + list of slides/labels.
###each row is concat([WSI sampN, methylation beta]); a slide is kept only if it has both a
###WSI sampler entry and a methylation match (via slide->Sentrix) and label>-1.
###also returns n_wsi (WSI feature count), needed to split X back apart for 'double_pca'.
def build_XY(sampD,methD,slides,labs):
    sampDslides=list(sampD.keys())
    n_wsi=len(np.ravel(sampD[sampDslides[0]]))

    beta=methD['beta']
    sentrix_ids=methD['sentrix_ids']
    sentrix_index={sx:i for i,sx in enumerate(sentrix_ids)}
    slide_to_sentrix=build_slide_to_sentrix(methD['metadata'])
    n_meth=beta.shape[1]

    snum=len(slides)
    fdim=n_wsi+n_meth
    X=np.zeros((snum,fdim),dtype=np.float32)
    Y=np.zeros((snum,),dtype=np.int32)

    cnt=0
    for i in range(snum):
        cslide=slides[i]
        clab=labs[i]
        sx=slide_to_sentrix.get(cslide)
        if (cslide in sampDslides) and (sx is not None) and (clab>-1):
            wsi_feat=np.ravel(sampD[cslide]).astype(np.float32)
            meth_feat=beta[sentrix_index[sx],:].astype(np.float32)
            X[cnt,:]=np.concatenate([wsi_feat,meth_feat],axis=0)
            Y[cnt]=clab
            cnt=cnt+1

    X=X[:cnt,:]
    Y=Y[:cnt]
    return X,Y,n_wsi


###train-fold-only PCA (fit on train, applied to train+test)
def pca_reduce(Xtrain,Xtest,n_components):
    ncomp=min(n_components,Xtrain.shape[0],Xtrain.shape[1])
    pca=PCA(n_components=ncomp,svd_solver='full',random_state=random_state).fit(Xtrain)
    pXtrain=pca.transform(Xtrain)
    pXtest=pca.transform(Xtest)
    return pXtrain,pXtest


###'concat_pca': one PCA fit on the full raw [WSI||meth] concat vector (2*PC_dim components)
def concat_pca_reduce(Xtrain,Xtest,cdim):
    return pca_reduce(Xtrain,Xtest,2*cdim)


###'double_pca': separate train-fold PCAs on the WSI block and the methylation block, then concat
def double_pca_reduce(Xtrain,Xtest,n_wsi,cdim):
    Xw_train,Xm_train=Xtrain[:,:n_wsi],Xtrain[:,n_wsi:]
    Xw_test,Xm_test=Xtest[:,:n_wsi],Xtest[:,n_wsi:]

    pXw_train,pXw_test=pca_reduce(Xw_train,Xw_test,cdim)
    pXm_train,pXm_test=pca_reduce(Xm_train,Xm_test,cdim)

    pXtrain=np.hstack([pXw_train,pXm_train])
    pXtest=np.hstack([pXw_test,pXm_test])
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
        elif kind=='qda':
            clf=QuadraticDiscriminantAnalysis(reg_param=params['reg_param'])
            clf.fit(Xtrain,Ytrain)
            probs=clf.predict_proba(Xtest)
        elif kind=='rf':
            clf=RandomForestClassifier(n_estimators=params['n_estimators'],max_depth=params['max_depth'],random_state=random_state,n_jobs=1)
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
    P['pca_modes']=pca_modes
    P['num_trees']=num_trees
    P['rf_max_depths']=rf_max_depths
    P['QDA_reg_vals']=QDA_reg_vals
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


###run every configured classifier on one already-PCA-reduced (pXtrain,pXtest), appending to results
def run_classifiers(P,pXtrain,Ytrain,pXtest,Ytest,cdim,pca_mode,results):
    if 'lda' in P['classifiers']:
        M=fit_and_score('lda',{},pXtrain,Ytrain,pXtest,Ytest)
        R={'classifier':'lda','PC_dim':cdim,'pca_mode':pca_mode}
        R.update(M)
        results.append(R)

    if 'qda' in P['classifiers']:
        for rp in P['QDA_reg_vals']:
            M=fit_and_score('qda',{'reg_param':rp},pXtrain,Ytrain,pXtest,Ytest)
            R={'classifier':'qda','PC_dim':cdim,'pca_mode':pca_mode,'reg_param':rp}
            R.update(M)
            results.append(R)

    if 'rf' in P['classifiers']:
        for ntree in P['num_trees']:
            for depth in P['rf_max_depths']:
                M=fit_and_score('rf',{'n_estimators':ntree,'max_depth':depth},pXtrain,Ytrain,pXtest,Ytest)
                R={'classifier':'rf','PC_dim':cdim,'pca_mode':pca_mode,'n_estimators':ntree,'max_depth':depth}
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
        R={'classifier':'lr','PC_dim':cdim,'pca_mode':pca_mode}
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
        R={'classifier':'mlp','PC_dim':cdim,'pca_mode':pca_mode}
        R.update(M)
        results.append(R)

    return results


###main entry point: classify concat([WSI,methylation]) reps given a train csv + a test csv
def classify_integrative(train_csv,test_csv,sampler_reps_file,methylation_reps_file,params=None):
    P=resolve_params(params)

    sampD=load_sampler_reps(sampler_reps_file)
    methD=load_methylation_reps(methylation_reps_file)
    train_slides,train_labs=load_split_csv(train_csv)
    test_slides,test_labs=load_split_csv(test_csv)

    Xtrain,Ytrain,n_wsi=build_XY(sampD,methD,train_slides,train_labs)
    Xtest,Ytest,_n_wsi=build_XY(sampD,methD,test_slides,test_labs)

    results=[]

    if Xtrain.shape[0]<2 or Xtest.shape[0]<1 or len(np.unique(Ytrain))<2:
        return results

    for cdim in P['PC_dims']:
        if 'concat_pca' in P['pca_modes']:
            pXtrain,pXtest=concat_pca_reduce(Xtrain,Xtest,cdim)
            results=run_classifiers(P,pXtrain,Ytrain,pXtest,Ytest,cdim,'concat_pca',results)

        if 'double_pca' in P['pca_modes']:
            pXtrain,pXtest=double_pca_reduce(Xtrain,Xtest,n_wsi,cdim)
            results=run_classifiers(P,pXtrain,Ytrain,pXtest,Ytest,cdim,'double_pca',results)

    return results
