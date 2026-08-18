### RMS subtype (ERMS / ARMS / Spindle) classification — thin wrapper around the same
### core code used for the (disease-level) tumor classification tasks
### (core_code/wsi_tumor_classification_core.py, methylation_tumor_classification_core.py,
### integrative_tumor_classification_core.py). Nothing RMS-specific lives in the core code:
### this script just points it at the RMS subtype splits (same schema: split_i ->
### slides_train, slides_lab_train, slides_test, slides_lab_test — 3 classes, ERMS=0/ARMS=1/
### Spindle=2 — see rms_subtype_cv_splits.p) and a modality ('wsi' | 'methylation' |
### 'integrative').
###
### Unlike the disease-level tasks, RMS subtyping does not use QDA or Random Forest (matches
### the real RMS subtype pipeline, which only sweeps LDA / logistic regression / MLP) —
### enforced here by restricting 'classifiers' to ['lda','lr','mlp'] by default.
###
### Usage:
###   python run_rms_subtype_classification.py --modality wsi --sampler-reps ../data/univ2_20x_sample_reps.p
###   python run_rms_subtype_classification.py --modality methylation --methylation-reps ../data/meth_top5000_sample_reps.p
###   python run_rms_subtype_classification.py --modality integrative --sampler-reps ../data/univ2_20x_sample_reps.p --methylation-reps ../data/meth_top5000_sample_reps.p
###   (add --n-splits 3 for a quick test; default splits file: rms_subtype_cv_splits.p next to this script)

import argparse
import os
import pickle
import sys

SCRIPT_DIR=os.path.dirname(os.path.abspath(__file__))
REPO_DIR=os.path.dirname(SCRIPT_DIR)
TC_DIR=os.path.join(REPO_DIR,'tumor_classification')
if TC_DIR not in sys.path:
    sys.path.insert(0,TC_DIR)

import run_wsi_tumor_classification as wsi_wrap  # noqa: E402
import run_methylation_tumor_classification as meth_wrap  # noqa: E402
import run_integrative_tumor_classification as integ_wrap  # noqa: E402

DEFAULT_SPLITS=os.path.join(SCRIPT_DIR,'rms_subtype_cv_splits.p')

### no QDA / RF for RMS subtyping (matches the real RMS subtype pipeline: mlda/logreg/mlp only)
RMS_CLASSIFIERS=['lda','lr','mlp']


def run_rms_subtype(modality,splits_file=DEFAULT_SPLITS,sampler_reps_file=None,methylation_reps_file=None,n_splits=None,params=None):
    P={'classifiers':RMS_CLASSIFIERS}
    if params is not None:
        for k in params.keys():
            P[k]=params[k]

    if modality=='wsi':
        all_results=wsi_wrap.run_all_splits(splits_file,sampler_reps_file,n_splits=n_splits,params=P)
        summarize=wsi_wrap.summarize_results
    elif modality=='methylation':
        all_results=meth_wrap.run_all_splits(splits_file,methylation_reps_file,n_splits=n_splits,params=P)
        summarize=meth_wrap.summarize_results
    elif modality=='integrative':
        all_results=integ_wrap.run_all_splits(splits_file,sampler_reps_file,methylation_reps_file,n_splits=n_splits,params=P)
        summarize=integ_wrap.summarize_results
    else:
        raise ValueError('unknown modality: '+modality)

    summary=summarize(all_results)
    return all_results,summary


def main():
    ap=argparse.ArgumentParser(description='RMS subtype (ERMS/ARMS/Spindle) classification wrapper')
    ap.add_argument('--modality',required=True,choices=['wsi','methylation','integrative'])
    ap.add_argument('--splits',default=DEFAULT_SPLITS,help='RMS subtype splits pickle')
    ap.add_argument('--sampler-reps',default=None,help='WSI sampler representations pickle (wsi/integrative)')
    ap.add_argument('--methylation-reps',default=None,help='methylation representations pickle (methylation/integrative)')
    ap.add_argument('--out',default='rms_subtype_class_results.p')
    ap.add_argument('--n-splits',type=int,default=None,help='limit number of splits (e.g. for a quick test)')
    args=ap.parse_args()

    if args.modality in ('wsi','integrative') and args.sampler_reps is None:
        raise ValueError('--sampler-reps is required for modality='+args.modality)
    if args.modality in ('methylation','integrative') and args.methylation_reps is None:
        raise ValueError('--methylation-reps is required for modality='+args.modality)

    all_results,summary=run_rms_subtype(
        args.modality,
        splits_file=args.splits,
        sampler_reps_file=args.sampler_reps,
        methylation_reps_file=args.methylation_reps,
        n_splits=args.n_splits,
    )

    D={}
    D['all_results']=all_results
    D['summary']=summary
    pickle.dump(D,open(args.out,'wb'))

    print(summary)
    print('wrote '+args.out)


if __name__=='__main__':
    main()
