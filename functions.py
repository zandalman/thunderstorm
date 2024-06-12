import numpy as np

def mag(vec): 
    return np.sqrt(np.sum(vec**2))

def exponential(x, a, b):
    return a * b**x

def interp(x0, x, y, logx=False, logy=False):
    ''' 
    Interpolate data to a given point. 
    
    Args
    x0: point to interpolate
    x, y: data points
    logx, logy: interpolate in logspace
    '''
    if logx: x0, x = np.log(x0), np.log(x)
    if logy: y = np.log(y)
    if x0 < x[0]:
        y0 = x[0] + (y[1]-y[0])/(x[1]-x[0]) * (x0 - x[0])
    elif x0 > x[-1]:
        y0 = x[-1] + (y[-1]-y[-2])/(x[-1]-x[-2]) * (x0 - x[-1])
    else:
        y0 = np.interp(x0, x, y)
    if logy: y0 = np.exp(y0)
    return y0