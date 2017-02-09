import numpy as np
import urllib.request
import subprocess

"""Well-known Boston housing dataset. We just add a 1-line header and do a
deterministic shuffle and split.

"""

url = "http://www.dcc.fc.up.pt/~ltorgo/Regression/housing.tar.gz"
filename = "housing.tar.gz"
urllib.request.urlretrieve(url, filename)

subprocess.call(["gunzip", filename])
subprocess.call(["tar", "xf", filename[:11]])

d = np.genfromtxt("Housing/housing.data", delimiter=",")
np.random.seed(0)
np.random.shuffle(d)

idx = int(0.7 * len(d))
d_train = d[:idx]
d_test = d[idx:]

colnames = "#CRIM ZN INDUS CHAS NOX RM AGE DIS RAD TAX PTRATIO B LSTAT CLASS"
np.savetxt("Housing-Train.csv", d_train, delimiter=",", header=colnames)
np.savetxt("Housing-Test.csv", d_test, delimiter=",", header=colnames)
