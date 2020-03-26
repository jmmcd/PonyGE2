import numpy as np
import json

def arc_from_json(filename):

    with open(filename, 'r') as f:
        obj = json.loads(f.readline())
        
    n_train = len(obj['train'])
    n_test  = len(obj['test'])
    
    x_train = []
    y_train = []
    x_test  = []
    y_test  = []
    
    for d in obj['train']:
        x_train.append(np.array(d['input']))
        y_train.append(np.array(d['output']))
    for d in obj['test']:
        x_test.append(np.array(d['input']))
        y_test.append(np.array(d['output']))
        
    return x_train, y_train, x_test, y_test