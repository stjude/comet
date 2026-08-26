### Wrapper for core RMS survival, methylation only (see
### core_code/rms_survival_methylation_core.py). Same design as run_rms_survival_sampler.py
### — see that file's docstring — just swaps in the methylation core + representations pickle.
###
### Usage:
###   python run_rms_survival_methylation.py --metadata ../data/RMS_survival_metadata.p \
###       --splits rms_survival_cv_splits.p --methylation-reps ../data/meth_top5000_sample_reps.p \
###       --out rms_survival_methylation_results.p
###   python run_rms_survival_methylation.py ... --n-splits 5   # quick test

import argparse
import os
import pickle
import sys

import pandas as pd

SCRIPT_DIR=os.path.dirname(os.path.abspath(__file__))
REPO_DIR=os.path.dirname(SCRIPT_DIR)
CORE_DIR=os.path.join(REPO_DIR,'core_code')
if CORE_DIR not in sys.path:
    sys.path.insert(0,CORE_DIR)

from rms_survival_methylation_core import classify_rms_survival_methylation  # noqa: E402
from safe_io import resolve_path, safe_pickle_load, safe_pickle_dump  # noqa: E402

TEMP_DIR=os.path.join(SCRIPT_DIR,'temp_rms_survival_methylation_csv')


def safedir(masir):
    if not os.path.isdir(masir):
        os.mkdir(masir)
    return


###how many splits are present in the splits pickle
def count_splits(splitsD):
    N=0
    while ('split_'+str(N)) in splitsD:
        N=N+1
    return N


###write one split's train/test rows to CSV (columns: slide,sentrix,time,event,binary_3yr)
def write_split_csvs(cohort,idx,ni,tag,temp_dir):
    Y=cohort['Y']
    lab_ind=cohort['lab_ind']
    slides=cohort['slides']
    sentrix=cohort['sentrix']

    csv_path=os.path.join(temp_dir,'split_%d_%s.csv'%(ni,tag))
    D={
        'slide':[slides[i] for i in idx],
        'sentrix':[sentrix[i] for i in idx],
        'time':[Y[i,0] for i in idx],
        'event':[Y[i,1] for i in idx],
        'binary_3yr':[Y[i,lab_ind] for i in idx],
    }
    pd.DataFrame.from_dict(D).to_csv(csv_path,index=False)
    return csv_path


###run the core Cox pipeline over every split, book-keeping results
def run_all_splits(metadata_file,splits_file,methylation_reps_file,n_splits=None,params=None):
    cohort=safe_pickle_load(metadata_file)
    splitsD=safe_pickle_load(splits_file)

    N=n_splits
    if N is None:
        N=count_splits(splitsD)

    safedir(TEMP_DIR)

    all_results=[]
    for ni in range(N):
        print('split '+str(ni))
        split=splitsD['split_%d'%ni]
        train_csv=write_split_csvs(cohort,split['train_idx'],ni,'train',TEMP_DIR)
        test_csv=write_split_csvs(cohort,split['test_idx'],ni,'test',TEMP_DIR)

        results=classify_rms_survival_methylation(train_csv,test_csv,methylation_reps_file,params=params)
        for R in results:
            R['split']=ni
            all_results.append(R)

    return all_results


###aggregate mean/std of each metric across splits, grouped by modality + hyperparameters
def summarize_results(all_results):
    df=pd.DataFrame(all_results)
    if len(df)==0:
        return df

    metric_cols=['c_index','auc']
    non_key_cols=metric_cols+['split','ok','reason','n_train','n_test']
    key_cols=[c for c in df.columns if c not in non_key_cols]

    summary=df.groupby(key_cols,dropna=False)[metric_cols].agg(['mean','std']).reset_index()
    return summary


def main():
    ap=argparse.ArgumentParser(description='RMS survival wrapper (methylation only)')
    ap.add_argument('--metadata',required=True,help='cleaned cohort metadata pickle')
    ap.add_argument('--splits',required=True,help='train/test splits pickle')
    ap.add_argument('--methylation-reps',required=True,help='methylation representations pickle')
    ap.add_argument('--out',default='rms_survival_methylation_results.p')
    ap.add_argument('--n-splits',type=int,default=None,help='limit number of splits (e.g. for a quick test)')
    args=ap.parse_args()

    all_results=run_all_splits(args.metadata,args.splits,args.methylation_reps,n_splits=args.n_splits)
    summary=summarize_results(all_results)

    D={}
    D['all_results']=all_results
    D['summary']=summary
    safe_pickle_dump(D,args.out)

    print(summary)
    print('wrote '+args.out)


if __name__=='__main__':
    main()
