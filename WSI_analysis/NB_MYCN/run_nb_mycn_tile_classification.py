### NB MYCN amplification (binary) tile-level classification — thin wrapper around
### core_code/tile_level_classification_core.py (per-tile logistic regression / 2-layer
### MLP, trained on tumor tiles only, slide prediction = mean tile-softmax; see
### github_repo/data/gen_nb_tiles.py). Tumor-tile filtering uses
### github_repo/data/NB_tumor_indicator.p (see data/gen_nb_tumor_indicator.py) instead of
### the core's RMS default. Same design as NB_MYCN/run_nb_mycn_mil.py: reads a train/test
### SPLITS pickle (schema: split_i -> slides_train, slides_lab_train, slides_test,
### slides_lab_test), writes per-split train/test CSVs (columns: slide,label), hands them
### to the core classifier, and books-keeps the metrics across all splits/hyperparameters.
###
### Usage:
###   python run_nb_mycn_tile_classification.py --tiles-dir ../data/NB_tiles --out results.p
###   python run_nb_mycn_tile_classification.py --tiles-dir ../data/NB_tiles --n-splits 5   # quick test
###   python run_nb_mycn_tile_classification.py --tiles-dir ../data/NB_tiles --device cuda  # force GPU

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

from tile_level_classification_core import classify_tile_level  # noqa: E402

DEFAULT_SPLITS=os.path.join(SCRIPT_DIR,'nb_mycn_cv_splits.p')
DEFAULT_TUM_INDICATOR=os.path.join(REPO_DIR,'data','NB_tumor_indicator.p')
TEMP_DIR=os.path.join(SCRIPT_DIR,'temp_nb_mycn_tile_csv')


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


###run the core tile-level classifier over every split, book-keeping all
###(split,classifier,hyperparam) results
def run_all_splits(splits_file,tiles_dir,n_splits=None,params=None):
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

        results=classify_tile_level(train_csv,test_csv,tiles_dir,params=params)
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
    ap=argparse.ArgumentParser(description='NB MYCN amplification tile-level classification wrapper')
    ap.add_argument('--splits',default=DEFAULT_SPLITS,help='NB MYCN splits pickle')
    ap.add_argument('--tiles-dir',required=True,help='tile-features directory (one pickle per slide)')
    ap.add_argument('--tum-indicator',default=DEFAULT_TUM_INDICATOR,help='tumor-tile indicator pickle')
    ap.add_argument('--device',default='auto',help='cpu | cuda | auto (default: auto-detect)')
    ap.add_argument('--out',default='nb_mycn_tile_results.p')
    ap.add_argument('--n-splits',type=int,default=None,help='limit number of splits (e.g. for a quick test)')
    args=ap.parse_args()

    params={'device':args.device,'tum_indicator_path':args.tum_indicator}
    all_results=run_all_splits(args.splits,args.tiles_dir,n_splits=args.n_splits,params=params)
    summary=summarize_results(all_results)

    D={}
    D['all_results']=all_results
    D['summary']=summary
    pickle.dump(D,open(args.out,'wb'))

    print(summary)
    print('wrote '+args.out)


if __name__=='__main__':
    main()
