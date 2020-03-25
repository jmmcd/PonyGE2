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
    x1, y1 = point1
    x2, y2 = point2
    x_incr = 1 if x2 >= x1 else -1
    y_incr = 1 if y2 >= y1 else -1
    return [(x, y) for x in range(x1, x2+x_incr, x_incr)
                   for y in range(y1, y2+y_incr, y_incr)]

def neighborhood(point, radius):
    x, y = point
    return [(x+dx, y+dy) for dx in range(x-radius, x+radius+1)
                         for dy in range(y-radius, y+radius+1)]

def selection_equals(m, sel1, sel2):
    colors1 = [ m[x, y] for x, y in sel1 ]
    colors2 = [ m[x, y] for x, y in sel2 ] 
    return colors1 == colors2

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
