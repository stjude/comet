### Wrapper for core WSI tumor classification (see core_code/wsi_tumor_classification_core.py).
### Reads a train/test SPLITS pickle (schema: split_i -> slides_train, slides_lab_train,
### slides_test, slides_lab_test), writes per-split train/test CSVs, hands them to the
### core classifier, and books-keeps the metrics across all splits/hyperparameters.
###
### Usage:
###   python run_wsi_tumor_classification.py --splits splits.p --sampler-reps ../data/univ2_20x_sample_reps.p --out results.p
###   python run_wsi_tumor_classification.py --splits splits.p --sampler-reps ../data/univ2_20x_sample_reps.p --n-splits 3   # quick test

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

from wsi_tumor_classification_core import classify_wsi_sampler  # noqa: E402

TEMP_DIR=os.path.join(SCRIPT_DIR,'temp_wsi_tumor_class_csv')


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


###write one split's train/test slide+label lists to CSV files (columns: slide,label)
def write_split_csvs(D,ni,temp_dir):
    train_csv=os.path.join(temp_dir,'split_'+str(ni)+'_train.csv')
    test_csv=os.path.join(temp_dir,'split_'+str(ni)+'_test.csv')

    Dtr={}
    Dtr['slide']=D['slides_train']
    Dtr['label']=D['slides_lab_train']
    pd.DataFrame.from_dict(Dtr).to_csv(train_csv,index=False)

    Dte={}
    Dte['slide']=D['slides_test']
    Dte['label']=D['slides_lab_test']
    pd.DataFrame.from_dict(Dte).to_csv(test_csv,index=False)

    return train_csv,test_csv


###run the core classifier over every split, book-keeping all (split,classifier,hyperparam) results
def run_all_splits(splits_file,sampler_reps_file,n_splits=None,params=None):
    splitsD=pickle.load(open(splits_file,'rb'))

    N=n_splits
    if N is None:
        N=count_splits(splitsD)

    safedir(TEMP_DIR)

    all_results=[]
    for ni in range(N):
        print('split '+str(ni))
        D=splitsD['split_'+str(ni)]
        train_csv,test_csv=write_split_csvs(D,ni,TEMP_DIR)

        results=classify_wsi_sampler(train_csv,test_csv,sampler_reps_file,params=params)
        for R in results:
            R['split']=ni
            all_results.append(R)

    return all_results


###aggregate mean/std of each metric across splits, grouped by classifier + hyperparameters
def summarize_results(all_results):
    df=pd.DataFrame(all_results)
    if len(df)==0:
        return df

    metric_cols=['ovr_auc','accuracy','balanced_accuracy','f1_macro']
    non_key_cols=metric_cols+['confusion_matrix','split','ok','reason']
    key_cols=[c for c in df.columns if c not in non_key_cols]

    summary=df.groupby(key_cols,dropna=False)[metric_cols].agg(['mean','std']).reset_index()
    return summary


def main():
    ap=argparse.ArgumentParser(description='WSI sampler-reps tumor classification wrapper')
    ap.add_argument('--splits',required=True,help='train/test splits pickle')
    ap.add_argument('--sampler-reps',required=True,help='sampler representations pickle')
    ap.add_argument('--out',default='wsi_tumor_class_results.p')
    ap.add_argument('--n-splits',type=int,default=None,help='limit number of splits (e.g. for a quick test)')
    args=ap.parse_args()

    all_results=run_all_splits(args.splits,args.sampler_reps,n_splits=args.n_splits)
    summary=summarize_results(all_results)

    D={}
    D['all_results']=all_results
    D['summary']=summary
    pickle.dump(D,open(args.out,'wb'))

    print(summary)
    print('wrote '+args.out)


if __name__=='__main__':
    main()
