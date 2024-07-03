import numpy as np
from scipy.special import erf

def calc_conf_per(data, sig=1):
    frac = erf(sig/np.sqrt(2))/2
    mean = np.mean(data, axis=0)
    llim = np.percentile(data, 100*(1-frac)/2, axis=0)
    ulim = np.percentile(data, 100*(1+frac)/2, axis=0)
    return mean, llim, ulim

def calc_conf_boot(data, sig=1, num_sample=1000, estimator=np.mean):
    frac = erf(sig / np.sqrt(2)) / 2
    est = estimator(data, axis=0)
    num_part, num_data = data.shape
    idx_boot = np.random.randint(0, num_part, size=(num_sample, num_part))
    est_boot = estimator(data[idx_boot], axis=1)
    llim = np.percentile(est_boot, (1 - frac) / 2 * 100, axis=0)
    ulim = np.percentile(est_boot, (1 + frac) / 2 * 100, axis=0)
    return est, llim, ulim
