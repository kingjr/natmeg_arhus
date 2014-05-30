print("################################################################################")
print("# Implementation of a multivariate pattern analysis based on  the scikitlearn   ")
print("# toolbox (http://scikit-learn.org/stable/). It reads a matlab file containing  ")
print("#     Xm:      a matrix of trials x chans x timepoint.                          ")
print("#     y:       a vector indicating the class of each trial                      ")
print("# The classification algorithm is based on a support vector machine.            ")
print("# (c) Jean-Remi King 2012, jeanremi.king [at] gmail.com                         ")
print("################################################################################")

################################################################################
print("LIBRARY")
import sys as sys
import numpy as np
from scipy import stats
from sklearn import svm
from sklearn.cross_validation import StratifiedKFold, LeaveOneOut, KFold
from sklearn.feature_selection import SelectPercentile, f_classif
import scipy.io as sio
from sklearn.preprocessing import Scaler

################################################################################
print("INPUT DATA")
#-- get argument to load specific file
filenameX = str(sys.argv[1])
filenamey = str(sys.argv[2])

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
generalize_time = mat["generalize_time"]  # svm penalization parameter
generalize_time = np.reshape(generalize_time, generalize_time.size)
fs_n = mat["fs"]  # feature selection
fs_n = np.reshape(fs_n, fs_n.size)
dims = mat["dims"]  # select time windows to compute
dims = np.reshape(dims, dims.size) - 1  # reshape for skl compatibility
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

n_samples, n_features, n_dims = Xm.shape
n_samplesg, n_featuresg, n_dimsg = Xmg.shape
n_classes = np.unique(y).shape[0]
#deal with sample_weight
sample_weight = np.ones(y.shape[0])
classes = np.unique(y2)
for c in range(classes.shape[0]):
    sample_weight[y2 == classes[c]] = 1. / (np.sum(y2 == classes[c]))

################################################################################
print("PREPARE CLASSIFICATION")

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

#-- feature selection
fs = SelectPercentile(f_classif, percentile=fs_n)
#-- grid search
#parameters = {'svm__C': (1e-6,1e-3, 1e-1, .4)}
#clf       = GridSearchCV(svm, parameters,n_jobs=1)

#-- initialize results
predict = np.zeros([n_splits, n_samples, n_dims]) ** np.nan
probas = np.zeros([n_splits, n_samples, n_dims, n_classes]) ** np.nan
predictg = np.zeros([n_splits, n_samplesg, n_dimsg, n_folds]) ** np.nan
probasg = np.zeros([n_splits, n_samplesg, n_dimsg, n_classes, n_folds]) ** np.nan
coef = np.empty([n_splits, n_folds, n_dims, n_classes * (n_classes - 1) / 2, n_features]) ** 0
all_folds = np.zeros([n_splits, n_folds, n_samples]) ** np.nan
y_shfl = np.copy(y)
Xm_shfl = np.copy(Xm)
sw_shfl = np.copy(sample_weight)


################################################################################
print("CLASSIFY...")
#-- shufflesplit
# repeat stratified kfolding for getting rid off the folding artefacts
for split in range(n_splits):
    print(split)
    # shuffle order
    new_order = np.array(range(y.shape[0]))
    if split > 0:
        np.random.shuffle(new_order)
        y_shfl[new_order] = np.copy(y)
        Xm_shfl[new_order, :, :] = np.copy(Xm)
        sw_shfl[new_order] = np.copy(sample_weight)
        cv = StratifiedKFold(y_shfl, k=n_folds)
    # Stratified crossvalidation
    for fold, (train, test) in enumerate(cv):
        print(fold)
        all_folds[split, fold, train] = 1
        all_folds[split, fold, test] = 0
        for d in range(0, dims.shape[0]):
            Xtrain = Xm_shfl[train, :, dims[d]]
            ytrain = y_shfl[train]
            sw_train = sw_shfl[train]
            # (deal with NaN in training)
            ytrain = ytrain[~np.isnan(np.nansum(Xtrain, axis=1))]
            sw_train = sw_train[~np.isnan(np.nansum(Xtrain, axis=1))]
            Xtrain = Xtrain[~np.isnan(np.nansum(Xtrain, axis=1)), :]
            if np.unique(ytrain).shape[0] > 1:
                # feature selection (find the 50% most discriminative channels)
                fs.fit(Xtrain, ytrain)         # find
                Xtrain = fs.transform(Xtrain)  # remove unnecessary channels
                # normalization
                scaler.fit(Xtrain)            # find
                Xtrain = scaler.transform(Xtrain)  # apply zscore
                # SVM fit
                clf.fit(Xtrain, ytrain, sample_weight=sw_train)
                # retrieve hyperplan feature identification
                coef[split, fold, dims[d], :, :] = 0  # initialize
                #--- univariate
                uni_features = fs.pvalues_ <= stats.scoreatpercentile(fs.pvalues_, fs.percentile)
                #--- multivariate
                coef[split, fold, dims[d], :, uni_features] = clf.coef_.T
                # predict cross val (deal with NaN in testing)
                Xtest = Xm_shfl[test, :, dims[d]]
                test_nan = np.isnan(np.nansum(Xtest, axis=1))
                Xtest = fs.transform(Xtest)
                Xtest = scaler.transform(Xtest)
                if (Xtest.shape[0] - np.sum(test_nan)) > 0:
                    if compute_predict:
                        predict[split, test[~test_nan], dims[d]] = clf.predict(Xtest[~test_nan, :])
                    if compute_probas:
                        probas[split, test[~test_nan], dims[d], :] = clf.predict_proba(Xtest[~test_nan, :])
                        if np.sum(test_nan) > 0:
                            probas[split, test[test_nan], dims[d], :] = np.nan
                # predict cross val on generalization sample (deal with NaN in testing)
                Xtestg = Xmg[:, :, dims[d]]
                test_nan = np.isnan(np.nansum(Xtestg, axis=1))
                Xtestg = fs.transform(Xtestg)
                Xtestg = scaler.transform(Xtestg)
                if (Xtestg.shape[0] - np.sum(test_nan)) > 0:
                    if compute_predict:
                        predictg[split, ~test_nan, dims[d], fold] = clf.predict(Xtestg[~test_nan, :])
                    if compute_probas:
                        probasg[split, ~test_nan, dims[d], :, fold] = clf.predict_proba(Xtestg[~test_nan, :])
                        if np.sum(test_nan) > 0:
                            probasg[split, test_nan, dims[d], :, fold] = np.nan
    #-- reorder results
    predict[split, :, :] = predict[split, new_order, :]
    probas[split, :, :, :] = probas[split, new_order, :, :]
    all_folds[split, :, :] = all_folds[split, :, new_order].T


################################################################################
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

output = path + nameX + '_' + namey + "_results.mat"
print(output)
sio.savemat(output, mat)
