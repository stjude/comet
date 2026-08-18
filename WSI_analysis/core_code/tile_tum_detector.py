### Tile-level tumor detector: a pre-trained LDA classifier on Virchow-v2 class-token
### features (base3 scale index 0, 20x, 1280-dim) that scores each tile as tumor/non-tumor.
### Used upstream of tile-level classification (see RMSsubtype_tile/build_rmst_tile_cache.py
### for the production usage this mirrors) to restrict training/eval to tumor tiles only.
###
### The classifier artifact itself (github_repo/data/tum_mlda_clf.p) is a small pre-fit
### sklearn LinearDiscriminantAnalysis (weights only, no patient data).

import os
import pickle

import numpy as np

SCRIPT_DIR=os.path.dirname(os.path.abspath(__file__))
REPO_DIR=os.path.dirname(SCRIPT_DIR)
DEFAULT_CLF_PATH=os.path.join(REPO_DIR,'data','tum_mlda_clf.p')

tum_prob_threshold=0.5


def load_tum_clf(clf_path=DEFAULT_CLF_PATH):
    return pickle.load(open(clf_path,'rb'))


###vclass_feats: (n_tiles,1280) float32 Virchow-v2 class-token features (scale index 0)
###returns tumor probability per tile
def predict_tumor_probs(vclass_feats,clf):
    X=np.ascontiguousarray(vclass_feats).astype(np.float32)
    probs=clf.predict_proba(X)[:,1]
    return probs


###returns a boolean tumor mask (probs>threshold)
def predict_tumor_mask(vclass_feats,clf,threshold=None):
    if threshold is None:
        threshold=tum_prob_threshold
    probs=predict_tumor_probs(vclass_feats,clf)
    return probs>threshold
