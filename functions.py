import numpy as np

def mag(vec): 
    return np.sqrt(np.sum(vec**2))

def exponential(x, a, b):
    return a * b**x