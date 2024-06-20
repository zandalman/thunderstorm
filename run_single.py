# import libraries for I/O
import sys, os
import argparse
from types import SimpleNamespace

# import libraries for computations
import numpy as np
from scipy.integrate import trapz as integrate
from scipy.integrate import cumtrapz as cumintegrate
import uuid

# import custom functions
import const
from functions import interp
from sim import Sim

# parse command line arguments
parser = argparse.ArgumentParser(prog='run.py', description='Run a Monte Carlo electron transport simulation.')
parser.add_argument('npac', type=int, help='number of packets')
parser.add_argument('rho', type=float, help='density [g/cc]')
parser.add_argument('ener', type=float, help='energy [eV]')
parser.add_argument('tmax', type=float, help='simulation duration [s]')
parser.add_argument('-t', '--time', type=float, default=10., help='time [day]')
parser.add_argument('-r', '--run', type=int, default=1, help='KNe simulation run')
parser.add_argument('-rp', '--rproc', type=bool, default=True, help='strong r-process')
parser.add_argument('-bt', '--bturb', type=float, default=0., help='turbulent magnetic field amplitude')
parser.add_argument('-bc', '--bco', type=float, default=0., help='coherent magnetic field amplitude')
parser.add_argument('-o', '--outfile', type=str, default='data.out', help='output file')

def append(filename, data):
    with open(filename, 'a') as f:
        f.write(data)

if __name__ == "__main__":

    # read command line arguments
    args = parser.parse_args()
    time = args.time * const.day
    rproc = 'strong' if args.rproc else 'weak'

    # remove output file if it exists
    if os.path.exists(args.outfile): os.remove(args.outfile)

    # start packet count
    npac = 0

    # read abundances and beta decay spectra from r-process simulations
    data_rproc = np.load(os.path.join('%s.npz' % rproc), allow_pickle=True)
    data_rproc = SimpleNamespace(**data_rproc)
    idx_rproc = np.searchsorted(data_rproc.time, time)
    spec_beta_cdf = cumintegrate(data_rproc.beta_spec[idx_rproc], data_rproc.ener, initial=0)/integrate(data_rproc.beta_spec[idx_rproc], data_rproc.ener)

    if time < data_rproc.time[0]:  print('Warning: time below minimum value in r-process simulations')
    if time > data_rproc.time[-1]: print('Warning: time above maximum value in r-process simulations')

    append(args.outfile, 'npac: %d\n' % args.npac)
    append(args.outfile, 'rho: %.3g g/cc\n' % args.rho)
    append(args.outfile, 'ener: %.3g keV\n' % (args.ener/1e3))
    append(args.outfile, '\n')
    append(args.outfile, 'packet id, step, event label, event category, time [s], x [cm], y [cm], z [cm], data1, data2\n')
    
    # create simulation object
    sim = Sim(ener=args.ener, rho=args.rho, do_moller=False)
    for Zelem in range(1, 99):
        ab = data_rproc.ab_elem[idx_rproc, Zelem]
        if ab == 0: continue
        A_avg = data_rproc.A_avg[idx_rproc, Zelem]
        sim.add_spec(Z=Zelem, A=A_avg, ab=ab)
    if args.bturb > 0: sim.add_Bturb(Bmag=args.bturb)
    if args.bco > 0: sim.add_Bco(Bmag=args.bco, Bhat=np.array([0, 1, 0]))

    while True:

        id_pac = npac
        
        # reset simulation
        sim.reset()

        # run the simulation
        nstep = 0
        while sim.time < args.tmax:
            sim.step()
            nstep += 1
        
        # write data
        data = ''
        for i, ev in enumerate(sim.ev_list):
            data += '%d,%d,%s,%s,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n' % (id_pac, nstep, ev.label, ev.cat, ev.time, ev.x, ev.y, ev.z, ev.data1, ev.data2)
        append(args.outfile, data)

        # increment packet counts
        npac += 1
        if npac >= args.npac: break
