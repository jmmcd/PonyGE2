import numpy as np

def get(m, point):
    sizex, sizey = m.shape
    x, y = point
    if (0 <= x < sizex) and (0 <= y < sizey):
        return int(m[x, y])
    else:
        return None

def get_sel(m, sel):
    sizex, sizey = m.shape
    cropped_sel = [ (x, y) for x, y in sel if (0 <= x < sizex) and (0 <= y < sizey) ]
    return [ int(m[x, y]) for x, y in cropped_sel ]

def mode(m, sel):
    color_dist = np.zeros(10)
    for col in get_sel(m, sel):
        color_dist[col] += 1
    return np.argmax(color_dist)

def first(m, sel):
    return get(m, sel[0])

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
    colors1 = [ get(m, point) for point in sel1 ]
    colors2 = [ get(m, point) for point in sel2 ] 
    return colors1 == colors2

def contains(m, sel, col):
    return col in get_sel(m, sel)

def any_hidden(m, sel):
    return len(get_sel(m, sel)) < len(sel)

def all_hidden(m, sel):
    return len(get_sel(m, sel)) == 0
