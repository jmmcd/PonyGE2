def mode(m, sel):
    pass

def first(m, sel):
    pass

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
    pass

def any_hidden(m, sel):
    pass

def all_hidden(m, sel):
    pass
