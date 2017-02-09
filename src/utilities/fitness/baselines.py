#!/usr/bin/env python
import sys
from collections import Counter

import numpy as np
from sklearn.linear_model import LinearRegression, ElasticNet

from utilities.fitness.get_data import get_data


def pprint(a, format_string ='{0:.2f}'):
    # adapted from http://stackoverflow.com/a/18287838
    return "[" + ", ".join(format_string.format(v,i) for i,v in enumerate(a)) + "]"


def fit_maj_class(train_X, train_y, test_X):
    """Use the majority class, for a binary problem..."""
    train_y = train_y.astype(int)
    classes = set(train_y) # often just {0, 1} or {-1, 1}
    maj = Counter(train_y).most_common(1)[0][0] # majority class
    model = "Majority class %d" % maj
    yhat_train = maj * np.ones(len(train_y))
    yhat_test = maj * np.ones(len(test_y))
    return model, yhat_train, yhat_test
    

def fit_const(train_X, train_y, test_X):
    """Use the mean of the y training values as a predictor."""
    mn = np.mean(train_y)
    yhat_train = np.ones(len(train_y)) * mn
    yhat_test = np.ones(len(test_y)) * mn
    model = "Const %.2f" % mn
    return model, yhat_train, yhat_test


def fit_lr(train_X, train_y, test_X):
    """Use linear regression to predict."""
    lr = LinearRegression()
    lr.fit(train_X, train_y)
    yhat_train = lr.predict(train_X)
    yhat_test = lr.predict(test_X)
    model = "LR int %.2f coefs %s" % (lr.intercept_, pprint(lr.coef_))
    return model, yhat_train, yhat_test

    
def fit_enet(train_X, train_y, test_X):
    """Use linear regression to predict -- elastic net is LR with L1
    and L2 regularisation."""
    enet = ElasticNet()
    enet.fit(train_X, train_y)
    model = "ElasticNet int %.2f coefs %s" % (enet.intercept_, pprint(enet.coef_))
    yhat_train = enet.predict(train_X)
    yhat_test = enet.predict(test_X)
    return model, yhat_train, yhat_test


    
if __name__ == "__main__":

    dataset_name = sys.argv[1]
    if len(sys.argv) > 2:
        metric = sys.argv[2]
    else:
        metric = "rmse"

    s = "from .error_metric import " + metric + " as metric"
    exec(s)

    train_X, train_y, test_X, test_y = get_data(dataset_name)
    train_X = train_X.T
    test_X = test_X.T

    methods = [fit_maj_class, fit_const, fit_lr, fit_enet]
    for fit in methods:
        model, train_yhat, test_yhat = fit(train_X, train_y, test_X)
        error_train = metric(train_y, train_yhat)
        error_test = metric(test_y, test_yhat)
        print("%s %s %s train error %.2f test error %.2f" %
              (metric.__name__, fit.__name__, model, error_train, error_test))
