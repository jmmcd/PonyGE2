import numpy as np

def mode(m, sel):
    color_dist = np.zeros(10)
    for x, y in sel:
        color = m[x, y]
        if color is not None:
            color_dist[color] += 1

    return np.argmax(color_dist) 

def first(m, sel):
    return m[sel[0][0], sel[0][1]]

def rectangle(point1, point2):
    pass

def neighborhood(point, radius):
    pass

def selection_equals(m, sel1, sel2):
    c_1 = np.array(list(map(lambda x, y: m[x, y], sel1)))
    c_2 = np.array(list(map(lambda x, y: m[x, y], sel2)))
    return np.all(c_1 == c_2)

def any_hidden(m, sel):
    for x, y in sel:
        if m[x, y] is None:
            return True

    return False

def all_hidden(m, sel):
    for x, y in sel:
        if m[x, y] is not None:
            return False

    return True