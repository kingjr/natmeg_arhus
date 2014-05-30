
print("######################################################################")
print("# Parallel n-split k-stratified-fold continuous SVM Scikitlearn MVPA #")
print("# (c) Jean-Remi King 2012, jeanremi.king [at] gmail [dot] com        #")
print("######################################################################")

# Implementation of a multivariate pattern analysis based on the scikit-learn
# toolbox (http://scikit-learn.org/stable/). It reads two .mat files
# (filenameX, filenamey) created by 'jr_classify.m'
#
# Function:
#     skl_king_parallel.py filenameX filenamey [number_of_cores]
#
# Inputs:
#     in filenameX:
#       Xm:     samples x features x classification matrix (e.g. trials x
#               chans x time)
#     in filenamey:
#       y:      vector indicating the class of each sample. Negative values
#               will be used for generalization only. 0 indicates to-be-
#               ignored samples.
#       y2:     cost/weights applied on each sample
#       path:   export directory
#       nameX:  export filename X
#       namey:  export filename y
#       folding:type of folding(e.g. stratified)
#       n_splits:number of splits
#       n_folds: number of folds
#       C:       SVM penalization parameter
#       compute_probas:     compute logit fit
#       compute_predict:    compute traditional SVM
#       fs_n:   number of univariate features selected for classification
#       dims:   classification performed on dims dimensions
#       dims_tg:classification generalized on dims_tg dimensions
#
# Ouputs:
#       predict: prediction matrix (split x samples x dims x dimsg)
#       predictg:same as predict for generalized samples
#       probas:  probas matrix (split x samples x dims x dimsg x class)
#       probasg: same as probas for generalized samples
#       coef:    weight hyperplan vector
#       all_folds:folding report (split x fold x samples)
#       y_all:   original y
#       y:       training y
#       yg:      generalized y
#       filenameX:
#       filenamey:
#
# Results are reported in: path + nameX + '_' + namey + "_results.mat"
###############################################################################
# (c) Jean-Remi King: jeanremi.king [at] gmail [dot] com
###############################################################################

# update    2013 01 03: input binary format
# update    2012 12 20: remove np.copy, add compute distance
# update    2012 11 29: fix 3rd dimension issue
# update    2012 11 13: fix bug str output on some python versions
# update    2012 11 02: change stratified kfolding y by y2
# update    2012 11 02: add np.copy to Xtrain and Xtest
# update    2012 11 01: correct feature selection coef bug when at 100 %
# update    2012 10 23: correct leaveoneout bug
# update    2012 10 23: correct major n_split new_order error
# update    2012 10 18: correct python/matlab dim incompatibility
# update    2012 10 18: correct error fs between 99 and 100 && remove Kbest
# update    2012 10 17: correct error n_features shape and add nice
# update    2012 10 01: correct prediction error+change loading results option
# update    2012 09 14: handle fs float error
# update    2012 09 14: pass n_cores to sys.arg
# version   2012 09 13: implementation of parallelization

###############################################################################
print("LIBRARY")
import sys as sys
import numpy as np
from scipy import stats
from sklearn import svm
from sklearn.cross_validation import StratifiedKFold, LeaveOneOut, KFold
from sklearn.feature_selection import SelectPercentile, SelectKBest, f_classif
from sklearn.externals.joblib import Parallel, delayed
import scipy.io as sio
from sklearn.preprocessing import Scaler
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import precision_score

###############################################################################
print("INPUT DATA")
#-- get argument to load specific file
filenameX = str(sys.argv[1])
filenamey = str(sys.argv[2])
if len(sys.argv) <= 3:
    n_cores = -1
else:
    n_cores = int(sys.argv[3])

print("cores: " + str(n_cores))
print(filenameX)
print(filenamey)

#-- load data dimension
mat = sio.loadmat(filenamey)
Xdims = mat["Xdims"].reshape(3)

#-- Load binary data into python
Xm_all = np.fromfile(filenameX, dtype=np.float64)
Xm_all = Xm_all.reshape(Xdims[2], Xdims[1], Xdims[0]).transpose([2, 1, 0])

#-- load classification parameters
mat = sio.loadmat(filenamey)
dims = mat["dims"]  # select time windows to compute
dims = np.reshape(dims, dims.size) - 1  # reshape for skl compatibility
dims_tg = mat["dims_tg"] - 1

mat = sio.loadmat(filenamey, squeeze_me=True)
path = mat["path"]
nameX = mat["nameX"]
namey = mat["namey"]
folding = mat["folding"]
n_splits = mat["n_splits"]
n_folds = mat["n_folds"]  # fold number
svm_C = mat["C"]  # svm penalization parameter
compute_probas = mat["compute_probas"]
compute_predict = mat["compute_predict"]
compute_distance = mat["compute_distance"]
fs_n = mat["fs"]  # feature selection
y_all = mat["y"]  # class used for train and test
print(Xm_all.shape)
print(y_all.shape)
y2_all = mat["y2"]  # class used for sample weights

#-- build training and generalizing classes
Xm = Xm_all[y_all > 0, :, :]  # training categories
Xmg = Xm_all[y_all < 0, :, :]  # generalization categories
y = y_all[y_all > 0]
yg = y_all[y_all < 0]
y2 = y2_all[y_all > 0]

n_samples, n_features, unused = Xm.shape
n_samplesg, unused, unused = Xmg.shape
n_featuresg = n_features
n_dims = dims.shape[0]
n_dimsg = n_dims

n_dims_tg = dims_tg.shape[1]
n_dimsg_tg = dims_tg.shape[1]

n_classes = np.unique(y).shape[0]
#deal with sample_weight
sample_weight = np.ones(y.shape[0])
classes = np.unique(y2)
for c in range(classes.shape[0]):
    sample_weight[y2 == classes[c]] = 1. / (np.sum(y2 == classes[c]))

###############################################################################
print("PREPARE CLASSIFICATION")


#-- classifier

clf = GridSearchCV(svm.SVC(kernel='linear', probability=True),
    {'C': svm_C}, score_func=precision_score)

#-- normalizer
scaler = Scaler()

#-- feature selection
if fs_n > 1:
    fs = SelectKBest(f_classif, k=fs_n)
elif fs_n == -1:
    fs = SelectKBest(f_classif, k=1)
else:
    fs = SelectPercentile(f_classif, percentile=fs_n * 100)


#-- results initialization
if compute_predict:
    predict = np.zeros([n_splits, n_samples, n_dims, n_dims_tg]) ** np.nan
    predictg = np.zeros([n_splits, n_samplesg, n_dimsg, n_dimsg_tg, n_folds]) ** np.nan
else:
    predict = []
    predictg = []

if compute_probas:
    probas = np.zeros([n_splits, n_samples, n_dims, n_dims_tg, n_classes]) ** np.nan
    probasg = np.zeros([n_splits, n_samplesg, n_dimsg, n_dimsg_tg, n_classes, n_folds]) ** np.nan
else:
    probas = []
    probasg = []


if compute_distance:
    distance = np.zeros([n_splits, n_samples, n_dims, n_dims_tg, n_classes]) ** np.nan
    distanceg = np.zeros([n_splits, n_samplesg, n_dimsg, n_dimsg_tg, n_classes, n_folds]) ** np.nan
else:
    distance = []
    distanceg = []

coef = np.empty([n_splits, n_folds, n_dims, n_classes * (n_classes - 1) / 2, n_features]) ** 0
all_folds = np.zeros([n_splits, n_folds, n_samples]) ** np.nan

###############################################################################
#-- Define parallel cross validation


def my_pipeline(train, test,
    Xm_shfl, y_shfl, sw_shfl, Xmg,
    dims, fs, scaler, clf,
    n_samples, n_dims, n_dims_tg, n_classes):
    # indicate opened fold
    sys.stdout.write("<")
    sys.stdout.flush()
    # initialize results within a given fold
    if compute_predict:
        predict = np.zeros([n_samples, n_dims, n_dims_tg]) ** np.nan
        predictg = np.zeros([n_samplesg, n_dimsg, n_dimsg_tg]) ** np.nan
    else:
        predict = []
        predictg = []
    if compute_probas:
        probas = np.zeros([n_samples, n_dims, n_dims_tg, n_classes]) ** np.nan
        probasg = np.zeros([n_samplesg, n_dimsg, n_dimsg_tg, n_classes]) ** np.nan
    else:
        probas = []
        probasg = []
    if compute_distance:
        distance = np.zeros([n_samples, n_dims, n_dims_tg, n_classes]) ** np.nan
        distanceg = np.zeros([n_samplesg, n_dimsg, n_dimsg_tg, n_classes]) ** np.nan
    else:
        distance = []
        distanceg = []
    coef = np.empty([n_dims, n_classes * (n_classes - 1) / 2, n_features]) ** 0
    # apply different classification along dimension 0
    for d in range(0, dims.shape[0]):
        Xtrain = Xm_shfl[train, :, dims[d]]
        ytrain = y_shfl[train]
        sw_train = sw_shfl[train]
        # (deal with NaN samples in training)
        ytrain = ytrain[~np.isnan(np.nansum(Xtrain, axis=1))]
        sw_train = sw_train[~np.isnan(np.nansum(Xtrain, axis=1))]
        Xtrain = Xtrain[~np.isnan(np.nansum(Xtrain, axis=1)), :]
        if np.unique(ytrain).shape[0] > 1:
            # feature selection
            fs.fit(Xtrain, ytrain)
            Xtrain = fs.transform(Xtrain)
            # normalization
            scaler.fit(Xtrain)
            Xtrain = scaler.transform(Xtrain)
            # Grid search
            clf.fit(Xtrain, ytrain, sample_weight=sw_train, cv=5)
            # select best classifier
            _clf = clf.best_estimator_
            # print(clf.best_params_)
            # retrieve features selected during univariate selection
            if fs_n == -1 or (fs_n > 1):
                # uni_features = sorted(range(len(fs.pvalues_)),key=lambda x:fs.pvalues_[x])
                uni_features = range(0, _clf.coef_.shape[1])
            else:
                uni_features = fs.pvalues_ <= stats.scoreatpercentile(fs.pvalues_, fs.percentile)
            # retrieve hyperplan (unselected features as 0)
            coef[d, :, uni_features] = scaler.inverse_transform(_clf.coef_).T
            # generalize across all time points
            for d_tg in range(0, n_dims_tg):
                # select data
                Xtest = Xm_shfl[test, :, dims_tg[d, d_tg]]
                # handles NaNs
                test_nan = np.isnan(np.nansum(Xtest, axis=1))
                Xtest = Xtest[~test_nan, :]
                # feature selection from training
                Xtest = fs.transform(Xtest)
                # normalize from training
                Xtest = scaler.transform(Xtest)
                # generalize test samples
                if (Xtest.shape[0] - np.sum(test_nan)) > 0:
                    if compute_predict:
                        predict[test[~test_nan], d, d_tg] = _clf.predict(Xtest)
                    if compute_probas:
                        probas[test[~test_nan], d, d_tg, :] = _clf.predict_proba(Xtest)
                    if compute_distance:
                        distance[test[~test_nan], d, d_tg, :] = _clf.decision_function(Xtest)  # correct!
                # predict on generalization sample
                # select data
                Xtestg = Xmg[:, :, dims_tg[d, d_tg]]
                # handles NaNs
                test_nan = np.isnan(np.nansum(Xtestg, axis=1))
                if (Xtestg.shape[0] - np.sum(test_nan)) > 0:
                    Xtestg = Xtestg[~test_nan, :]
                    # preproc feature selection and normalization
                    Xtestg = fs.transform(Xtestg)
                    Xtestg = scaler.transform(Xtestg)
                    # compute prediction
                    if compute_predict:
                        predictg[~test_nan, d, d_tg] = clf.predict(Xtestg)
                    if compute_probas:
                        probasg[~test_nan, d, d_tg, :] = clf.predict_proba(Xtestg)
                    if compute_distance:
                        distanceg[~test_nan, d, d_tg, :] = clf.decision_function(Xtestg)  # correct!
    # summarize fold results
    out = {
        'coef': coef,
        'predict': predict,
        'predictg': predictg,
        'probas': probas,
        'probasg': probasg,
        'distance': distance,
        'distanceg': distanceg}
    # indicate end of fold
    sys.stdout.write(">")
    sys.stdout.flush()
    return out

###############################################################################
print("CLASSIFY")
#-- Shuffle split
for split in range(n_splits):
    print("split " + str(split))
    #-- shuffle order in case this is not the first split
    new_order = np.array(range(y.shape[0]))
    if split > 0:
        np.random.shuffle(new_order)
        y_shfl = y
        y_shfl = y_shfl[new_order]
        y2_shfl = y2
        y2_shfl = y2_shfl[new_order]
        Xm_shfl = Xm
        Xm_shfl = Xm_shfl[new_order, :, :]
        sw_shfl = sample_weight
        sw_shfl = sw_shfl[new_order]
    else:
        y_shfl = y
        y2_shfl = y2
        Xm_shfl = Xm
        sw_shfl = sample_weight
    #-- define crossvalidation
    if folding == 'stratified':
        cv = StratifiedKFold(y2_shfl, k=n_folds)
    elif folding == 'kfolding':
        cv = KFold(n=y2_shfl.shape[0], k=n_folds)
    elif folding == 'leaveoneout':
        n_folds = y_shfl.shape[0]
        cv = LeaveOneOut(n=y_shfl.shape[0])
    else:
        print("unknown crossvalidation method!")
    # Cross-validation computed in parallel
    out = Parallel(n_jobs=n_cores)(delayed(my_pipeline)(
        train=train,
        test=test,
        Xm_shfl=Xm_shfl,
        y_shfl=y_shfl,
        sw_shfl=sw_shfl,
        Xmg=Xmg,
        dims=dims,
        fs=fs,
        scaler=scaler,
        clf=clf,
        n_samples=n_samples,
        n_dims=n_dims,
        n_dims_tg=n_dims_tg,
        n_classes=n_classes) for fold, (train, test) in enumerate(cv))
    # reorder results folds and splits
    for fold, (train, test) in enumerate(cv):
        all_folds[split, fold, train] = 1
        all_folds[split, fold, test] = 0
        coef[split, fold, :, :, :] = out[fold]['coef']
        if compute_predict:
            predict[split, new_order[test], :, :] = out[fold]['predict'][test, :, :]
            predictg[split, :, :, :, fold] = out[fold]['predictg']
        if compute_probas:
            probas[split, new_order[test], :, :, :] = out[fold]['probas'][test, :, :, :]
            probasg[split, :, :, :, :, fold] = out[fold]['probasg']
        if compute_distance:
            distance[split, new_order[test], :, :, :] = out[fold]['distance'][test, :, :, :]
            distanceg[split, :, :, :, :, fold] = out[fold]['distanceg']
    all_folds[split, :, new_order] = all_folds[split, :, :].T

###############################################################################
print("EXPORT DATA")
mat['predict'] = predict
mat['predictg'] = predictg
mat['probas'] = probas
mat['probasg'] = probasg
mat['distance'] = distance
mat['distanceg'] = distanceg
mat['coef'] = coef
mat['all_folds'] = all_folds
mat['y_all'] = y_all
mat['y'] = y
mat['yg'] = yg
mat['filenameX'] = filenameX
mat['filenamey'] = filenamey

print nameX
print namey
print path

output = str(path) + str(nameX) + '_' + str(namey) + "_results.mat"

print(output)
sio.savemat(output, mat)
