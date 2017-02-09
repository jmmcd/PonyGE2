import numpy as np
from sklearn.metrics.classification import f1_score as sklearn_f1_score


def mae(y, yhat):
    """Calculate mean absolute error between inputs."""
    return np.mean(np.abs(y - yhat))
mae.maximise = False


def rmse(y, yhat):
    """Calculate root mean square error between inputs."""
    return np.sqrt(np.mean(np.square(y - yhat)))
rmse.maximise = False


def mse(y, yhat):
    """Calculate mean square error between inputs."""
    return np.mean(np.square(y - yhat))
mse.maximise = False


def hinge(y, yhat):
    """Hinge loss is a suitable loss function for classification.  Here y is
    the true values (-1 and 1) and yhat is the "raw" output of the individual,
    ie a real value. The classifier will use sign(yhat) as its prediction."""
    # NB not np.max. maximum does element-wise max
    return np.sum(np.maximum(0, 1 - y * yhat))
hinge.maximise = False


def f1_score(y, yhat):
    """The F_1 score is a metric for classification which tries to balance
    precision and recall, ie both true positives and true negatives.
    For F_1 score higher is better."""

    # convert real values to boolean with a zero threshold
    yhat = (yhat > 0)
    # this can give a runtime warning because if we predict the same
    # value for all samples (trivial individuals will do so), f-score
    # is undefined (see sklearn implementation). However it's a
    # warning, we can ignore.
    return sklearn_f1_score(y, yhat)
f1_score.maximise = True
