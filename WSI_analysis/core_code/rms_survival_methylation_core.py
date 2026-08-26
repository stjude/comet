### Core RMS survival code (methylation only): t-test feature screening + PCA + a
### penalized Cox proportional-hazards model on top-N-variable-probe methylation
### representations. Same Cox pipeline as core_code/rms_survival_sampler_core.py — see
### that file's docstring for the general design (also core_code/rms_survival_integrative_core.py).
###
### Takes a train CSV, a test CSV (columns: slide,sentrix,time,event,binary_3yr — see
### github_repo/data/gen_rms_survival_metadata.py), and a path to a methylation
### representations pickle (see github_repo/data/gen_meth_top5000_sample_reps.py):
###   {'probe_ids','sentrix_ids','beta':(N,fdim) float16,'metadata':{...}}
### Samples are joined to the beta matrix by the CSV's own 'sentrix' column (already
### resolved once in gen_rms_survival_metadata.py), not by re-deriving slide->Sentrix here.

import pickle

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score
from sksurv.metrics import concordance_index_censored
from lifelines import CoxPHFitter

###############################################################################
# Default hyperparameters (RMS survival, methylation)
###############################################################################
PCA_N_COMPONENTS=50
TTEST_P_THRESHOLD=0.05
cox_penalizer=0.1
cox_l1_ratio=0.5
random_state=0


###load methylation reps pickle: {'probe_ids','sentrix_ids','beta','metadata'}
def load_methylation_reps(methylation_reps_file):
    return pickle.load(open(methylation_reps_file,'rb'))


###load a train/test split csv (columns: slide,sentrix,time,event,binary_3yr)
def load_split_csv(csv_file):
    return pd.read_csv(csv_file)


###build X (methylation beta) + time/event/binary_3yr from a split dataframe
###(skip rows with no Sentrix match in the beta matrix or an invalid — negative — binary_3yr)
def build_XY(methD,df):
    beta=methD['beta']
    sentrix_index={sx:i for i,sx in enumerate(methD['sentrix_ids'])}

    rows,times,events,bys=[],[],[],[]
    for sx,t,e,by in zip(df['sentrix'],df['time'],df['event'],df['binary_3yr']):
        if sx not in sentrix_index or by<0:
            continue
        rows.append(beta[sentrix_index[sx],:].astype(np.float32))
        times.append(t)
        events.append(e)
        bys.append(by)

    fdim=beta.shape[1]
    if not rows:
        return (np.zeros((0,fdim),dtype=np.float32),np.zeros((0,),dtype=np.float32),
                np.zeros((0,),dtype=np.float32),np.zeros((0,),dtype=np.float32))

    X=np.vstack(rows).astype(np.float32)
    T=np.asarray(times,dtype=np.float32)
    E=np.asarray(events,dtype=np.float32)
    BY=np.asarray(bys,dtype=np.float32)
    return X,T,E,BY


###fit a penalized Cox model on the training PCA features + (T,E), score risk on test
###mirrors benchmark/utils.py:mycox — coef-weighted linear risk score, min-max normalized
def cox_predict_risk(Xtrain,Ttrain,Etrain,Xtest,penalizer,l1_ratio):
    fnum=Xtrain.shape[1]
    colnames=['f%d'%i for i in range(fnum)]+['T','E']
    df=pd.DataFrame(np.concatenate((Xtrain,Ttrain[:,None],Etrain[:,None]),axis=1),columns=colnames)

    cph=CoxPHFitter(penalizer=penalizer,l1_ratio=l1_ratio).fit(df,'T','E')
    coef=np.asarray(list(cph.summary['coef']))[np.newaxis,:]

    h_test=np.sum(coef*Xtest,axis=1)
    probs=(h_test-np.min(h_test))/(np.max(h_test)-np.min(h_test)+1e-8)
    return probs


###t-test feature screening (train-only) + PCA (train-fit, applied to train+test)
def ttest_pca_reduce(Xtrain,Xtest,BYtrain,pdim,pthresh):
    _,pvals=stats.ttest_ind(Xtrain[BYtrain==0],Xtrain[BYtrain==1],equal_var=False,nan_policy='omit')
    pvals=np.asarray(pvals,dtype=np.float64)
    sel=np.where(pvals<pthresh)[0]
    if len(sel)<pdim:
        return None,None,len(sel)

    ncomp=min(pdim,Xtrain.shape[0])
    pca=PCA(n_components=ncomp,svd_solver='full',random_state=random_state).fit(Xtrain[:,sel])
    pXtrain=pca.transform(Xtrain[:,sel]).astype(np.float32)
    pXtest=pca.transform(Xtest[:,sel]).astype(np.float32)
    return pXtrain,pXtest,len(sel)


###C-index + AUC(binary_3yr) given a continuous risk score on the test set
def survival_metrics(Ttest,Etest,BYtest,risk):
    cindex,_,_,_,_=concordance_index_censored(Etest>0,Ttest,risk)
    try:
        auc=float(roc_auc_score(BYtest,risk))
    except (ValueError,TypeError):
        auc=float('nan')
    return {'c_index':float(cindex),'auc':auc}


###fit+score one (PC_dim) setting; on failure return nan metrics + reason instead of
###raising, so one bad setting does not stop the whole run
def fit_and_score(pdim,pthresh,penalizer,l1_ratio,Xtrain,Ttrain,Etrain,BYtrain,Xtest,Ttest,Etest,BYtest):
    try:
        pXtrain,pXtest,n_sel=ttest_pca_reduce(Xtrain,Xtest,BYtrain,pdim,pthresh)
        if pXtrain is None:
            return {'c_index':float('nan'),'auc':float('nan'),'ok':False,
                    'reason':'ttest_features_lt_pca (n_sel=%d)'%n_sel}

        risk=cox_predict_risk(pXtrain,Ttrain,Etrain,pXtest,penalizer,l1_ratio)
        M=survival_metrics(Ttest,Etest,BYtest,risk)
        M['ok']=True
        M['reason']=None
    except Exception as exc:
        M={'c_index':float('nan'),'auc':float('nan'),'ok':False,'reason':str(exc)}

    return M


###merge user-supplied params dict with the module defaults (only overrides matching keys)
def resolve_params(params):
    P={}
    P['PCA_N_COMPONENTS']=PCA_N_COMPONENTS
    P['TTEST_P_THRESHOLD']=TTEST_P_THRESHOLD
    P['cox_penalizer']=cox_penalizer
    P['cox_l1_ratio']=cox_l1_ratio

    if params is not None:
        for k in params.keys():
            P[k]=params[k]

    return P


###main entry point: Cox survival on methylation reps given a train csv + a test csv
def classify_rms_survival_methylation(train_csv,test_csv,methylation_reps_file,params=None):
    P=resolve_params(params)

    methD=load_methylation_reps(methylation_reps_file)
    train_df=load_split_csv(train_csv)
    test_df=load_split_csv(test_csv)

    Xtrain,Ttrain,Etrain,BYtrain=build_XY(methD,train_df)
    Xtest,Ttest,Etest,BYtest=build_XY(methD,test_df)

    results=[]
    if Xtrain.shape[0]<10 or Xtest.shape[0]<5 or len(np.unique(BYtrain))<2:
        return results

    M=fit_and_score(P['PCA_N_COMPONENTS'],P['TTEST_P_THRESHOLD'],P['cox_penalizer'],P['cox_l1_ratio'],
                     Xtrain,Ttrain,Etrain,BYtrain,Xtest,Ttest,Etest,BYtest)
    R={'modality':'methylation','PC_dim':P['PCA_N_COMPONENTS'],
       'n_train':int(Xtrain.shape[0]),'n_test':int(Xtest.shape[0])}
    R.update(M)
    results.append(R)

    return results
