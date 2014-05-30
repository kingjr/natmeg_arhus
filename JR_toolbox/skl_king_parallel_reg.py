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


# update    2013 10 18: correction wsize dims_tg default values
# update    2013 09 20: save_allg: compress gen data for memory 
# update    2013 09 20: StandardScaler update
# update    2013 06 18: windows compatibility, dims_tg=0  bug correction
# update    2013 03 16: optional save coef; need to finish save all g
# update    2013 01 26: correct TOI time generalization bug for single samples
# update    2013 01 24: add wsize (not adapted for time gen)+ default parameters
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
from sklearn.preprocessing import StandardScaler
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
 

def default(var, value, dict):
    if var in dict:
        out = dict[var]
    else:
        out = value
    return out
 
#-- classification parameters
mat = sio.loadmat(filenamey, squeeze_me=True)
y_all = mat["y"]  # class used for train and test
path = default("path", "", mat)
nameX = default("nameX", "default", mat)
namey = default("namey", "default", mat)
folding = default("folding", "stratified", mat)
n_splits = default("n_splits", 1, mat)
n_folds = default("n_folds", 5, mat)
svm_C = default("C", 1, mat)
compute_probas = default("compute_probas", True, mat)
compute_predict = default("compute_predict", False, mat)
compute_distance = default("compute_distance", False, mat)
fs_n = default("fs", 1, mat)
y2_all = default("y2", y_all, mat)
wsize = default("wsize", 1, mat)
fb_clf = default("fb_clf", False, mat)
fb_split = default("fb_split", True, mat)
save_coef = default("save_coef", True, mat)
save_allg = default("save_allg", True, mat)
fb_fold = default("fb_fold", True, mat)
trans = default("transpose", [0, 1, 2], mat)
Xdim = mat["Xdim"]
 
#-- Load X data
Xm_all = np.fromfile(filenameX, dtype=np.float32)
Xm_all = Xm_all.reshape(Xdim[2], Xdim[1], Xdim[0]).transpose([2, 1, 0])
Xm_all = Xm_all.transpose(trans) # specific sample feature clf dimensions
mat2 = sio.loadmat(filenamey)
 
#-- X-specific parameters
if "dims" in mat2:
    dims = mat2["dims"]  # select time windows to compute
    dims = np.reshape(dims, dims.size) - 1  # reshape for skl compatibility
else:
    dims = np.array(range(0,Xm_all.shape[2]-wsize+1))
 
if "dims_tg" in mat2:
    dims_tg = mat2["dims_tg"]
    if np.prod(dims_tg.shape)==1:
        if dims_tg[0,0] == 0: # all time gen
            dims_tg = dims[np.newaxis]
            dims_tg = np.tile(dims_tg,[dims_tg.shape[1],1]);
    else:
        dims_tg = dims_tg - 1
        if np.shape(dims_tg.shape)[0] == 1:  # /!\ Matlab squeeze array
            dims_tg = dims_tg.reshape([1, dims_tg.shape[0]])
else:
    # dims_tg = np.transpose([range(0,Xm_all.shape[2])])
    dims_tg = dims[np.newaxis]
    dims_tg = dims_tg.T

mat["dims"] = dims
mat["dims_tg"] = dims_tg
 
#-- pre format X dimensions
features = default("features", range(0,Xm_all.shape[1]), mat)
Xm_all = Xm_all[:,features,:]  # a priori feature selection;test
 
###############################################################################
#-- build training and generalizing classes
print(Xm_all.shape)
print(y_all.shape)
Xm = Xm_all[y_all > -1000, :, :]  # training categories
Xmg = Xm_all[y_all < -1000, :, :]  # generalization categories
y = y_all[y_all > -1000]
yg = y_all[y_all < -1000]
y2 = y2_all[y_all > -1000]
 
n_samples, n_features, unused = Xm.shape
n_samplesg, unused, unused = Xmg.shape
n_featuresg = n_features
n_dims = dims.shape[0]
n_dimsg = n_dims
print(n_dims)

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
if np.size(svm_C)>1:
    clf = GridSearchCV(svm.SVC(kernel='linear', probability=compute_probas),
    {'C': svm_C}, score_func=precision_score)
else:
    #clf = svm.SVC(kernel='linear', probability=True, C=svm_C)
    clf = svm.SVR(kernel='linear', C=svm_C)
 
#-- normalizer
scaler = StandardScaler()
 
#-- feature selection
if fs_n > 1:
    fs = SelectKBest(f_classif, k=fs_n)
elif fs_n == -1:
    fs = SelectKBest(f_classif, k=1)
else:
    fs = SelectPercentile(f_classif, percentile=fs_n * 100)
 
#-- results initialization
(predict, predictg, probas, probasg, distance, distanceg) = ([], [], [], [], [], [])
if compute_predict:
    predict = np.zeros([n_splits, n_samples, n_dims, n_dims_tg]) ** np.nan
    predictg = np.zeros([n_splits, n_samplesg, n_dimsg, n_dimsg_tg, n_folds]) ** np.nan
if compute_probas:
    probas = np.zeros([n_splits, n_samples, n_dims, n_dims_tg, n_classes]) ** np.nan
    if save_allg:
        probasg = np.zeros([n_splits, n_samplesg, n_dimsg, n_dimsg_tg, n_classes, n_folds]) ** np.nan
    else:
        probasg = np.zeros([n_splits, n_samplesg, n_dims, n_dims_tg, n_classes]) ** np.nan
if compute_distance:
    distance = np.zeros([n_splits, n_samples, n_dims, n_dims_tg, n_classes]) ** np.nan
    distanceg = np.zeros([n_splits, n_samplesg, n_dimsg, n_dimsg_tg, n_classes, n_folds]) ** np.nan
if save_coef:
    coef = np.empty([n_splits, n_folds, n_dims, n_classes * (n_classes - 1) / 2, n_features * wsize]) ** 0
all_folds = np.zeros([n_splits, n_folds, n_samples]) ** np.nan
 
###############################################################################
#-- Define parallel cross validation
def my_pipeline(train, test,
    Xm_shfl, y_shfl, sw_shfl, Xmg,
    dims, fs, scaler, _clf, svm_C,
    n_samples, n_samplesg,
    n_dims, n_dims_tg, n_classes, wsize,
    fb_fold,fb_clf,save_coef):
    # indicate opened fold
    if fb_fold:
        sys.stdout.write("<")
        sys.stdout.flush()
    # initialize results within a given fold
    (predict, predictg, probas, probasg, distance, distanceg) = ([], [], [], [], [], [])
    if compute_predict:
        predict = np.zeros([n_samples, n_dims, n_dims_tg]) ** np.nan
        predictg = np.zeros([n_samplesg, n_dimsg, n_dimsg_tg]) ** np.nan
    if compute_probas:
        probas = np.zeros([n_samples, n_dims, n_dims_tg, n_classes]) ** np.nan
        probasg = np.zeros([n_samplesg, n_dimsg, n_dimsg_tg, n_classes]) ** np.nan
    if compute_distance:
        distance = np.zeros([n_samples, n_dims, n_dims_tg, n_classes]) ** np.nan
        distanceg = np.zeros([n_samplesg, n_dimsg, n_dimsg_tg, n_classes]) ** np.nan
    if save_coef:
        coef = np.empty([n_dims, n_classes * (n_classes - 1) / 2, n_features * wsize]) ** 0
    else:
        coef = []
    # apply different classification along dimension 0
    for d in range(0, dims.shape[0]):
        if fb_clf:
            sys.stdout.write(".")
            sys.stdout.flush()
        Xtrain = Xm_shfl[train, :, dims[d]:(dims[d] + wsize)]
        if np.shape(Xtrain.shape)[0] == 3:
            Xtrain = Xtrain.reshape(Xtrain.shape[0], Xtrain.shape[1] * Xtrain.shape[2])
        ytrain = y_shfl[train]
        sw_train = sw_shfl[train]
        # (deal with NaN samples in training)
        # ytrain = ytrain[~np.isnan(np.nansum(Xtrain, axis=1))]
        # sw_train = sw_train[~np.isnan(np.nansum(Xtrain, axis=1))]
        # Xtrain = Xtrain[~np.isnan(np.nansum(Xtrain, axis=1)), :]
        ytrain = ytrain[~np.isnan(np.sum(Xtrain, axis=1))]
        sw_train = sw_train[~np.isnan(np.sum(Xtrain, axis=1))]
        Xtrain = Xtrain[~np.isnan(np.sum(Xtrain, axis=1)), :]
        if np.unique(ytrain).shape[0] > 1:
            # feature selection
            fs.fit(Xtrain, ytrain)
            Xtrain = fs.transform(Xtrain)
            # normalization
            scaler.fit(Xtrain)
            Xtrain = scaler.transform(Xtrain)
            # SVM fit
            # Grid search
            if np.size(svm_C)>1:
                _clf.fit(Xtrain, ytrain, sample_weight=sw_train, cv=4)
                clf = _clf.best_estimator_
            else:
                clf = _clf
            clf.fit(Xtrain, ytrain, sample_weight=sw_train)
            if save_coef:
                # retrieve features selected during univariate selection
                uni_features = np.argsort(fs.pvalues_)[0:scaler.inverse_transform(clf.coef_).shape[1]]
                coef[d, :, uni_features] = scaler.inverse_transform(clf.coef_).T
            # generalize across all time points
            for d_tg in range(0, n_dims_tg):
                # select data
                Xtest = Xm_shfl[test, :, dims_tg[d, d_tg]:(dims_tg[d, d_tg] + wsize)]
                if np.shape(Xtest.shape)[0] == 3:
                    Xtest = Xtest.reshape(Xtest.shape[0], Xtest.shape[1] * Xtest.shape[2])  # adapt wsize for time generalization
                # handles NaNs
                # test_nan = np.isnan(np.nansum(Xtest, axis=1))
                test_nan = np.isnan(np.sum(Xtest, axis=1))
                Xtest = Xtest[~test_nan, :]
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
                    if compute_distance:
                        distance[test[~test_nan], d, d_tg, :] = clf.decision_function(Xtest)  # correct!
                # predict on generalization sample
                # select data
                Xtestg = Xmg[:, :, dims_tg[d, d_tg]:(dims_tg[d, d_tg] + wsize)]
                if np.shape(Xtestg.shape)[0] == 3:
                    Xtestg = Xtestg.reshape(Xtestg.shape[0], Xtestg.shape[1] * Xtestg.shape[2])  # adapt wsize for time generalization
                # handles NaNs
                # test_nan = np.isnan(np.nansum(Xtestg, axis=1))
                test_nan = np.isnan(np.sum(Xtestg, axis=1))
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
    if fb_fold:
        sys.stdout.write(">")
        sys.stdout.flush()
    return out
 
###############################################################################
print("CLASSIFY")
#-- Shuffle split
for split in range(n_splits):
    if fb_split:
        sys.stdout.write("*")
        sys.stdout.flush()
    #-- shuffle order in case this is not the first split
    new_order = np.array(range(y.shape[0]))
    if split >= 0:
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
        cv = StratifiedKFold(y2_shfl, n_folds=n_folds)
    elif folding == 'kfolding':
        cv = KFold(n=y2_shfl.shape[0], n_folds=n_folds)
    elif folding == 'leaveoneout':
        n_folds = y_shfl.shape[0]
        cv = LeaveOneOut(n=y_shfl.shape[0])
    else:
        print("unknown crossvalidation method!")
    # Cross-validation computed in parallel
    print(n_cores)
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
        _clf=clf,
        svm_C=svm_C,
        n_samples=n_samples,
        n_samplesg=n_samplesg,
        n_dims=n_dims,
        n_dims_tg=n_dims_tg,
        n_classes=n_classes,
        wsize=wsize,
        fb_fold=fb_fold,
        fb_clf=fb_clf,
        save_coef=save_coef) for fold, (train, test) in enumerate(cv))
    # reorder results folds and splits
    for fold, (train, test) in enumerate(cv):
        all_folds[split, fold, train] = 1
        all_folds[split, fold, test] = 0
        if save_coef:
            coef[split, fold, :, :, :] = out[fold]['coef']
        if compute_predict:
            predict[split, new_order[test], :, :] = out[fold]['predict'][test, :, :]
            predictg[split, :, :, :, fold] = out[fold]['predictg']
        if compute_probas:
            probas[split, new_order[test], :, :, :] = out[fold]['probas'][test, :, :, :]
            if save_allg:
                probasg[split, :, :, :, :, fold] = out[fold]['probasg']
            elif yg.shape[0]>0:
                if fold == 0:
                    probasg[split, :, :, :, :] = out[fold]['probasg']/n_folds
                else:
                    probasg[split, :, :, :, :] += out[fold]['probasg']/n_folds
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
if save_coef:
    mat['coef'] = coef
mat['all_folds'] = all_folds
mat['y_all'] = y_all
mat['y'] = y
mat['yg'] = yg
mat['filenameX'] = filenameX
mat['filenamey'] = filenamey
 
print(nameX)
print(namey)
print(path)
 
output = str(path) + str(nameX) + '_' + str(namey) + "_results.mat"
 
print(output)
sio.savemat(output, mat)
 
if False:  # for debugging
 split = 0
 if fb_split:
     sys.stdout.write("*")
     sys.stdout.flush()
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
     cv = StratifiedKFold(y2_shfl, n_folds=n_folds)
 elif folding == 'kfolding':
     cv = KFold(n=y2_shfl.shape[0], n_folds=n_folds)
 elif folding == 'leaveoneout':
     n_folds = y_shfl.shape[0]
     cv = LeaveOneOut(n=y_shfl.shape[0])
 else:
     print("unknown crossvalidation method!")
 # Cross-validation computed in parallel
 print(n_cores)

 for fold, (train, test) in enumerate(cv):
     out = my_pipeline(
        train=train,
        test=test,
        Xm_shfl=Xm_shfl,
        y_shfl=y_shfl,
        sw_shfl=sw_shfl,
        Xmg=Xmg,
        dims=dims,
        fs=fs,
        scaler=scaler,
        _clf=clf,
        svm_C=svm_C,
        n_samples=n_samples,
        n_samplesg=n_samplesg,
        n_dims=n_dims,
        n_dims_tg=n_dims_tg,
        n_classes=n_classes,
        wsize=wsize,
        fb_fold=fb_fold,
        fb_clf=fb_clf,
        save_coef=save_coef)