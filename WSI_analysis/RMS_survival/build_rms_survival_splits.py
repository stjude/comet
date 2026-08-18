#!/usr/bin/env python3
### Build stratified train/test CV splits for RMS survival, from the already-cleaned
### cohort metadata (see data/gen_rms_survival_metadata.py). Same design as
### RMSsurvival/build_rms_survival_splits.py:build_splits_dict (test_size=0.2, stratified
### on binary_3yr, fixed seed) — reimplemented here directly against our own metadata
### pickle rather than copied, because the shared minimal-example H&E sampler pickle
### (data/univ2_20x_sample_reps.p) covers 3 more slides than the production run behind the
### real published splits, so the two cohorts aren't index-for-index identical.
###
### Output: rms_survival_cv_splits.p
###   {'n_splits','test_size','random_seed_base','cohort_n_samples',
###    'split_i': {'train_idx','test_idx'}}  # indices into data/RMS_survival_metadata.p's
###                                          # 'Y'/'slides'/'sentrix'/'patient_ids' arrays

import os
import pickle

import numpy as np
from sklearn.model_selection import train_test_split

SCRIPT_DIR=os.path.dirname(os.path.abspath(__file__))
REPO_DIR=os.path.dirname(SCRIPT_DIR)
METADATA_FILE=os.path.join(REPO_DIR,'data','RMS_survival_metadata.p')
OUT_PATH=os.path.join(SCRIPT_DIR,'rms_survival_cv_splits.p')

N_SPLITS=200
TEST_SIZE=0.2
SPLIT_RANDOM_SEED=42


def main():
    cohort=pickle.load(open(METADATA_FILE,'rb'))
    Y=cohort['Y']
    lab_ind=cohort['lab_ind']
    snum=len(Y)
    inds=np.arange(snum)

    splitsD={
        'n_splits':int(N_SPLITS),
        'test_size':float(TEST_SIZE),
        'random_seed_base':int(SPLIT_RANDOM_SEED),
        'cohort_n_samples':int(snum),
    }

    rng=np.random.RandomState(SPLIT_RANDOM_SEED)
    n_fail=0
    for ni in range(N_SPLITS):
        fold_seed=int(rng.randint(0,2**31-1))
        try:
            inds_train,inds_test=train_test_split(
                inds,test_size=TEST_SIZE,stratify=Y[:,lab_ind],random_state=fold_seed,
            )
            splitsD['split_%d'%ni]={
                'train_idx':np.asarray(inds_train,dtype=np.int64),
                'test_idx':np.asarray(inds_test,dtype=np.int64),
            }
        except Exception:
            n_fail=n_fail+1

    pickle.dump(splitsD,open(OUT_PATH,'wb'),protocol=4)

    print('Wrote %s'%OUT_PATH)
    print('built %d/%d splits (failed=%d)'%(N_SPLITS-n_fail,N_SPLITS,n_fail))
    d0=splitsD['split_0']
    print('split_0: train=%d test=%d'%(len(d0['train_idx']),len(d0['test_idx'])))


if __name__=='__main__':
    main()
