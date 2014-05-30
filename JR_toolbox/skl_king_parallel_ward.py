
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
from sklearn.cluster import WardAgglomeration
import scipy.io as sio
from sklearn.preprocessing import Scaler

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

#-- Load data into python
mat = sio.loadmat(filenameX)
Xm_all = mat["Xm"]  # data

#-- load classification parameters
mat = sio.loadmat(filenamey)
path = mat["path"][0]
nameX = mat["nameX"][0]
namey = mat["namey"][0]
folding = mat["folding"][0]
n_splits = mat["n_splits"]  # svm penalization parameter
n_splits = np.reshape(n_splits, n_splits.size)
n_folds = mat["n_folds"]  # fold number
n_folds = np.reshape(n_folds, n_folds.size)
svm_C = mat["C"]  # svm penalization parameter
svm_C = np.reshape(svm_C, svm_C.size)
compute_probas = mat["compute_probas"]  # svm penalization parameter
compute_probas = np.reshape(compute_probas, compute_probas.size)
compute_predict = mat["compute_predict"]  # svm penalization parameter
compute_predict = np.reshape(compute_predict, compute_predict.size)
fs_n = mat["fs"]  # feature selection
fs_n = np.reshape(fs_n, fs_n.size)
dims = mat["dims"]  # select time windows to compute
dims = np.reshape(dims, dims.size) - 1  # reshape for skl compatibility
dims_tg = mat["dims_tg"] - 1  # svm penalization parameter
y_all = mat["y"]  # class used for train and test
y_all = np.reshape(y_all, y_all.size)  # reshape for skl compatibility
y2_all = mat["y2"]  # class used for sample weights
y2_all = np.reshape(y2_all, y2_all.size)  # reshape for skl compatibility

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
n_features = 50
#--crossvalidation
if folding == 'stratified':
    cv = StratifiedKFold(y, k=n_folds)
elif folding == 'kfolding':
    cv = KFold(n=y.shape[0], k=n_folds)
elif folding == 'leaveoneout':
    n_folds[0] = y.shape[0]
    cv = LeaveOneOut(n=y.shape[0])
else:
    print("unknown crossvalidation method!")


#-- classifier
clf = svm.SVC(kernel='linear', probability=True, C=svm_C)

#-- normalizer
scaler = Scaler()

#-- Clustering
n_clusters = 100
cluster = WardAgglomeration(n_clusters=n_clusters, connectivity=None,
    compute_full_tree='auto')

#-- feature selection
fs_n = round(n_features * fs_n) / n_features
if fs_n == 100:
    fs = SelectKBest(f_classif, k=n_features)
else:
    fs = SelectPercentile(f_classif, percentile=fs_n)


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

coef = np.empty([n_splits, n_folds, n_dims, n_classes * (n_classes - 1) / 2, n_features]) ** 0
all_folds = np.zeros([n_splits, n_folds, n_samples]) ** np.nan
y_shfl = np.copy(y)
Xm_shfl = np.copy(Xm)
sw_shfl = np.copy(sample_weight)


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
            # ward clustering
            cluster.fit(np.mean(Xtrain, axis=0))
            Xtrain = cluster.transform(Xtrain)
            # feature selection
            fs.fit(Xtrain, ytrain)
            Xtrain = fs.transform(Xtrain)
            # normalization
            scaler.fit(Xtrain)
            Xtrain = scaler.transform(Xtrain)
            # SVM fit
            clf.fit(Xtrain, ytrain, sample_weight=sw_train)
            # retrieve features selected during univariate selection
            uni_features = fs.pvalues_ <= stats.scoreatpercentile(fs.pvalues_, fs.percentile)
            #uni_features = cluster.inverse_transform(uni_features)
            # retrieve hyperplan (unselected features as 0)
            # coef[d, :, uni_features] = scaler.inverse_transform(clf.coef_).T
            # generalize across all time points
            for d_tg in range(0, n_dims_tg):
                # select data
                Xtest = Xm_shfl[test, :, dims_tg[d, d_tg]]
                # handles NaNs
                test_nan = np.isnan(np.nansum(Xtest, axis=1))
                Xtest = Xtest[~test_nan, :]
                # clustering
                Xtest = cluster.transform(Xtest)
                # feature selection from training
                Xtest = fs.transform(Xtest)
                # normalize from training
                Xtest = scaler.transform(Xtest)
                # generalize test samples
                if (Xtest.shape[0] - np.sum(test_nan)) > 0:
                    if compute_predict:
                        predict[test[~test_nan], d, d_tg] = clf.predict(Xtest)
                    if compute_probas:
                        probas[test[~test_nan], d, d_tg, :] = clf.predict_proba(Xtest)
                # predict on generalization sample
                # select data
                Xtestg = Xmg[:, :, dims_tg[d, d_tg]]
                # handles NaNs
                test_nan = np.isnan(np.nansum(Xtestg, axis=1))
                if (Xtestg.shape[0] - np.sum(test_nan)) > 0:
                    Xtestg = Xtestg[~test_nan, :]
                    # preproc feature selection and normalization
                    Xtestg = cluster.transform(Xtestg)
                    Xtestg = fs.transform(Xtestg)
                    Xtestg = scaler.transform(Xtestg)
                    # compute prediction
                    if compute_predict:
                        predictg[~test_nan, d, d_tg] = clf.predict(Xtestg)
                    if compute_probas:
                        probasg[~test_nan, d, d_tg, :] = clf.predict_proba(Xtestg)
    # summarize fold results
    out = {
        'coef': coef,
        'predict': predict,
        'predictg': predictg,
        'probas': probas,
        'probasg': probasg}
    # indicate end of fold
    sys.stdout.write(">")
    sys.stdout.flush()
    return out

###############################################################################
print("CLASSIFY")
#-- Shuffle split
for split in range(n_splits):
    print("split " + str(split))
    # shuffle order in case this is not the first split
    new_order = np.array(range(y.shape[0]))
    if split > 0:
        np.random.shuffle(new_order)
        y_shfl[new_order] = np.copy(y)
        Xm_shfl[new_order, :, :] = np.copy(Xm)
        sw_shfl[new_order] = np.copy(sample_weight)
        cv = StratifiedKFold(y_shfl, k=n_folds)
    # Cross-validation computed in parallel
    # run parallel computation
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
        n_classes=n_classes) for train, test in cv)
    # reorder results folds and splits
    for fold, (train, test) in enumerate(cv):
        all_folds[split, fold, train] = 1
        all_folds[split, fold, test] = 0
        coef[split, fold, :, :, :] = out[fold]['coef']
        if compute_predict:
            predict[split, test, :, :] = out[fold]['predict'][new_order[test], :, :]
            if out[fold]['predictg']:
                predictg[split, :, :, :, fold] = out[fold]['predictg']
        if compute_probas:
            probas[split, test, :, :, :] = out[fold]['probas'][new_order[test], :, :, :]
            if out[fold]['probasg']:
                probasg[split, :, :, :, :, fold] = out[fold]['probasg']
    all_folds[split, :, :] = all_folds[split, :, new_order].T


###############################################################################
print("EXPORT DATA")
mat['predict'] = predict
mat['predictg'] = predictg
mat['probas'] = probas
mat['probasg'] = probasg
mat['coef'] = coef
mat['all_folds'] = all_folds
mat['y_all'] = y_all
mat['y'] = y
mat['yg'] = yg
mat['filenameX'] = filenameX
mat['filenamey'] = filenamey

output = path + nameX + '_' + namey + "_ward_results.mat"
print(output)
sio.savemat(output, mat)
