import numpy as np

def mag(vec): 
    return np.sqrt(np.sum(vec**2))

def exponential(x, a, b):
    return a * b**x

def interp(x0, x, y, logx=False, logy=False, llim=None, ulim=None):
    ''' 
    Interpolate data to a given point. 
    
    Args
    x0: point(s) to interpolate
    x, y: data points
    logx, logy: interpolate in logspace
    '''
    is_arr = True
    if type(x0) != np.ndarray: x0 = np.array([x0]); is_arr=False
    if logx: x0, x = np.log(x0), np.log(x)
    if logy: y = np.log(y)
    cond_low  = x0 < x[0]
    cond_high = x0 > x[-1]
    cond_med  = np.logical_and(~cond_low, ~cond_high) 
    y0 = np.zeros_like(x0, dtype=float)
    y0[cond_low]  = y[0] if x[1]==x[0] else y[0] + (y[1]-y[0])/(x[1]-x[0]) * (x0[cond_low] - x[0])
    y0[cond_high] = y[-1] if x[-1] == x[-2] else y[-1] + (y[-1]-y[-2])/(x[-1]-x[-2]) * (x0[cond_high] - x[-1])
    y0[cond_med]  = np.interp(x0[cond_med], x, y)
    if logy: y0 = np.exp(y0)
    if llim is not None: y0[y0<llim] = llim
    if ulim is not None: y0[y0>ulim] = ulim
    if not is_arr: y0 = y0[0]
    return y0