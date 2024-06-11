# system modules
import os
from types import SimpleNamespace

# math modules
import numpy as np
from scipy.integrate import trapz as integrate
from scipy.integrate import cumtrapz as cumintegrate
from scipy.interpolate import interp1d

# custom modules
from functions import mag
from helper import *
import const

# periodic table
elements = SimpleNamespace(**np.load('elements.npz'))

# beta particle data
data_EEDL = SimpleNamespace(**np.load('EEDL.npz', allow_pickle=True))

class Sim(object):
    '''
    Monte Carlo charged particle transport simulation.
    
    Args
    
    ener: partlce energy [eV]
    elem: background element
    rho:  background density [g/cm^3]

    m_i:   particle mass [g]
    q_i:   particle charge [esu]
    
    do_Bco:   use a coherent magnetic field
    do_Bturb: use a turbulent magnetic field
    do_iso:   use isotropic scattering
    seed:     random number generator seed
    num:      number of interpolation points

    do_ion:   include ionization of background ions
    do_exc:   include excitation of background ions
    do_brem:  include Bremsstrahlung with background ions
    do_scat:  include elastic scattering with background ions
    '''
    def __init__(self, ener=300000, rho=1e-16, m_i=const.m_e, q_i=const.e, do_iso=False, seed=None, num=int(1e12), do_ion=True, do_exc=True, do_brem=True, do_scat=True):
        
        self.ener  = ener
        self.rho   = rho
        
        self.m_i  = m_i
        self.q_i  = q_i

        self.do_Bco_1  = False
        self.do_Bco_2 = False
        self.do_Bturb = False
        self.do_iso   = do_iso

        self.do_ion  = do_ion
        self.do_exc  = do_exc
        self.do_brem = do_brem
        self.do_scat = do_scat
        self.spec_list = []
        
        self.seed = seed
        self.rng  = np.random.default_rng(self.seed)
        self.num = num

        self.idx_event_list = []
        self.coord_list     = []
        self.time_list      = []
        
        self.time  = 0
        self.coord = np.zeros(3)
        self.vhat  = self.rand_dir()

    def reset(self):
        ''' Reset the simulation. '''
        self.idx_event_list = []
        self.coord_list     = []
        self.time_list      = []

        self.time  = 0
        self.coord = np.zeros(3)
        self.vhat  = self.rand_dir()  

    def add_Bturb(self, Bturb, q=5/3, fturb=1, Lmax=1e16):
        '''
        Add a turbulent magnetic field.

        Bturb: turbulent magnetic field amplitude [G]
        q:     power law exponent of magnetic turbulence spectrum
        fturb: ratio of turbulent to magnetic energy
        Lmax:  largest scale of magnetic turbulence [cm]
        '''
        self.Bturb = Bturb
        self.q     = q
        self.fturb = fturb
        self.Lmax  = Lmax
        self.do_Bturb = True
    
    def add_Bco(self, Bmag, Bhat, method=1):
        ''' 
        Add a coherent magnetic field. 
        
        Args
        Bmag: coherent magnetic field amplitude [G]
        Bhat: coherent magnetic field direction
        method:
            1. transport the particle in the direction of the coherent magnetic field 
            2. use the coherent magnetic field to add anisotropy to the turbulent field scattering
        '''
        self.Bmag = Bmag
        self.Bhat = Bhat
        if method == 1:
            self.do_Bco_1 = True
        else:
            assert self.do_Bturb == True, "Error: method 2 requires a turbulent magnetic field, but no such field was found"
            self.do_Bco_2 = True

            chi = self.Bmag/self.Bturb
            self.sin_be = np.linspace(-1, 1, 1024)
            cos_be = np.sqrt(1 - self.sin_be**2)
            self.sin_be_cdf = np.zeros_like(self.sin_be)

            cond = np.abs(self.sin_be) < 1/max(chi, 1)
            cos_al_p = -chi*self.sin_be[cond]**2 + cos_be[cond] * np.sqrt(1 - chi**2 * self.sin_be[cond]**2)
            cos_al_m = -chi*self.sin_be[cond]**2 - cos_be[cond] * np.sqrt(1 - chi**2 * self.sin_be[cond]**2)

            self.sin_be_cdf[cond] = 1/2 * (1 + np.sign(self.sin_be)[cond] * (1 + 1/2 * (cos_al_m - cos_al_p)))
            self.sin_be_cdf[~cond] = 1/2 * (1 + np.sign(self.sin_be)[~cond])

    def add_spec(self, Z, A, ab):
        ''' 
        Add a species to the background plasma. 
        
        Args
        Z:  atomic number
        A:  atomic mass
        ab: abundance (mass fraction)
        '''
        spec_data = SimpleNamespace(Z=Z, A=A, ab=ab)
        spec_data.name = elements.name[Z-1]
        spec_data.symbol = elements.symbol[Z-1]
        
        spec_data.sig_ion_name_list = []
        spec_data.sig_ion_func_list = []

        for i, dtype_name in enumerate(data_EEDL.dtype_list):

            data = data_EEDL.data[Z-1, i]
            if np.all(data) == None: continue

            if dtype_name[:3] == 'sig':

                ener, sig = data
                sig = sig*ab/A # effective cross section
                sig_func = interp1d(np.log(ener), np.log(sig), fill_value='extrapolate')
            
                if dtype_name == 'sig_ion' and len(dtype_name) > 7:

                    spec_data.sig_ion_name_list.append(dtype_name)
                    spec_data.sig_ion_func_list.append(sig_func)
                
                else:

                    setattr(spec_data, '%s_func' % dtype_name, sig_func)

            elif dtype_name == 'th_scat':
        
                spec_data.cos_th_cdf = np.logspace(-9, 0, self.num)
                cos_th_grid = np.zeros((len(data), self.num))

                ener_list = []
                for i, data_ener in enumerate(data):
                    
                    ener, cos_th, cos_th_pdf = data_ener
                    ener_list.append(ener)

                    cos_th_cdf     = cumintegrate(cos_th_pdf, cos_th, initial=0) / integrate(cos_th_pdf, cos_th)
                    cos_th_func    = interp1d(cos_th_cdf, cos_th, fill_value='extrapolate')
                    cos_th_grid[i] = cos_th_func(spec_data.cos_th_cdf)

                ener_list = np.array(ener_list)
                spec_data.cos_th_scat_func = interp1d(np.log(ener_list), cos_th_grid.T, fill_value='extrapolate')

        self.spec_list.append(spec_data)
    
    def calc_sig_Bturb(self):
        ''' Compute the effective cross section of the particle to a turbulent magnetic field. '''
        beta = self.vmag/const.c
        beta_A = self.Bturb / np.sqrt(4*np.pi*self.rho*const.c**2) # Alfven speed relative to the speed of light
        func_beta_A = (1-(beta_A/beta)**(2-self.q))/(2-self.q) - (1-(beta_A/beta)**(4-self.q))/(4-self.q)
        
        Om = self.q_i*self.Bturb/(self.m_i*const.c) # particle gyro-frequency
        kmin = const.c/(Om*self.Lmax) # minimum wavenumber of the magnetic turbulence spectrum
        
        lam_Bturb = const.c*beta**(2-self.q)*self.gam**(2-self.q)*func_beta_A * 2/(np.pi*(self.q-1)*self.fturb*const.c*kmin) * (const.c*kmin/Om)**(2-self.q)
        self.sig_Bturb = 1/(lam_Bturb*self.rho*const.N_A)

    def rand_dir(self):
        ''' Compute a direction by uniformally sampling the sphere. '''
        uhat = np.zeros(3)

        phi    = 2.0*np.pi*self.rng.random()
        cos_th = 2*self.rng.random()-1
        sin_th = np.sqrt(1-cos_th**2)

        uhat[X] = np.sin(phi)*sin_th
        uhat[Y] = np.cos(phi)*sin_th
        uhat[Z] = cos_th

        return uhat
    
    def get_event_list(self):
        ''' Create a list of possible interactions. '''
        self.sig_list   = []
        self.event_list = []

        for spec_data in self.spec_list:
            
            if self.do_ion:

                for i, sig_ion_name in enumerate(spec_data.sig_ion_name_list):
                
                    sig_ion = np.exp(spec_data.sig_ion_func_list[i](np.log(self.ener)))
                    self.sig_list.append(sig_ion)
                    self.event_list.append(SimpleNamespace(
                        func = self.dep_ener,
                        args = dict(),
                        spec = spec_data.Z,
                        cat = 'ionization',
                        label = '%s ionization (%s)' % (spec_data.symbol, sig_ion_name[8:])
                    ))

            if self.do_exc:

                sig_exc = np.exp(spec_data.sig_exc_func(np.log(self.ener)))
                self.sig_list.append(sig_exc)
                self.event_list.append(SimpleNamespace(
                    func = self.dep_ener,
                    args = dict(),
                    spec = spec_data.Z,
                    cat = 'excitation',
                    label = '%s excitation' % (spec_data.symbol)
                ))

            if self.do_brem:

                sig_brem = np.exp(spec_data.sig_brem_func(np.log(self.ener)))
                self.sig_list.append(sig_brem)
                self.event_list.append(SimpleNamespace(
                    func = self.dep_ener,
                    args = dict(),
                    spec = spec_data.Z,
                    cat = 'Bremsstrahlung',
                    label = '%s Bremsstrahlung' % (spec_data.symbol)
                ))

            if self.do_scat:

                cos_th_scat = spec_data.cos_th_scat_func(np.log(self.ener))
                sig_scat_la = np.exp(spec_data.sig_scat_la_func(np.log(self.ener)))
                self.sig_list.append(sig_scat_la)
                self.event_list.append(SimpleNamespace(
                    func = self.scat,
                    args = dict(do_iso=self.do_iso, trig_th_scat=cos_th_scat, trig_th_scat_cdf=spec_data.cos_th_cdf),
                    spec = spec_data.Z,
                    cat = 'elastic scatter',
                    label = '%s elastic scatter (large angle)' % (spec_data.symbol)
                ))

                sig_scat_sa = np.exp(spec_data.sig_scat_func(np.log(self.ener))) - sig_scat_la
                self.sig_list.append(sig_scat_sa)
                self.event_list.append(SimpleNamespace(
                    func = self.dep_ener,
                    args = dict(),
                    spec = spec_data.Z,
                    cat = 'elastic scatter',
                    label = '%s elastic scatter (small angle)' % (spec_data.symbol)
                ))

        if self.do_Bturb:
            
            self.calc_sig_Bturb()
            self.sig_list.append(self.sig_Bturb)
            if self.do_Bco_2:
                args = dict(do_iso=False, trig_th_scat=self.sin_be, trig_th_scat_cdf=self.sin_be_cdf, scat_dir=self.Bhat, do_sin=True)
            else:
                args = dict(do_iso=True)
            
            self.event_list.append(SimpleNamespace(
                func = self.scat,
                args = args,
                spec = -1,
                cat = 'Bturb scatter',
                label = 'Bturb scatter'
            ))

        self.sig_tot = np.sum(self.sig_list)

    def dep_ener(self):

        pass

    def scat(self, do_iso, trig_th_scat=None, trig_th_scat_cdf=None, scat_dir=None, do_sin=False):
        '''
        Scatter the particle packet:
            1. Sample the CDF of the scattering kernel to determine the scattering angle
            2. Compute a unit vector perpendicular to the velocity
            3. Rotate the velocity about the perpendicular vector by the scattering angle using the Rodrigues formula
            4. Rotate the velocity about the original velocity by a random angle

        do_iso:           use isotropic scattering
        trig_th_scat:     sin/cos of scattering angle
        trig_th_scat_cdf: cumulative distribution function of the sin/cos of scattering angle
        scat_dir:         unit vector representing the direction off which to scatter
        do_sin:           use sine instead of cosine
        '''
        if do_iso:
            
            self.vhat = self.rand_dir()
        
        else:
            
            if np.all(scat_dir) == None: scat_dir = self.vhat

            phi_scat = 2.0*np.pi*self.rng.random()
            if do_sin:
                sin_th_scat = trig_th_scat[np.searchsorted(trig_th_scat_cdf, self.rng.random())]
                cos_th_scat = np.sqrt(1-sin_th_scat**2)
            else:
                cos_th_scat = trig_th_scat[np.searchsorted(trig_th_scat_cdf, self.rng.random())]
                sin_th_scat = np.sqrt(1-cos_th_scat**2)
            
            vperp = np.zeros(3)
            vperp[X], vperp[Y] = scat_dir[Y], -scat_dir[X]
            vperp /= np.sqrt(scat_dir[X]**2 + scat_dir[Y]**2)
            
            vhatnew = scat_dir * cos_th_scat      + np.cross(vperp, scat_dir)   * sin_th_scat      + vperp    * np.sum(vperp*scat_dir)   * (1-cos_th_scat)
            vhatnew = vhatnew  * np.cos(phi_scat) + np.cross(scat_dir, vhatnew) * np.sin(phi_scat) + scat_dir * np.sum(scat_dir*vhatnew) * (1-np.cos(phi_scat))
            self.vhat = vhatnew/mag(vhatnew)
    
    def move(self):
        ''' Move the packet until the next interaction. '''
        dis = -np.log(1.0-self.rng.random()) / (self.rho*const.N_A*self.sig_tot) # sample distance to collision by inversion technique
        self.time += dis/self.vmag

        if self.do_Bco_1:
            self.coord += dis*np.sum(self.vhat*self.Bhat)*self.Bhat
        else:
            self.coord += dis*self.vhat

    def interaction(self):
        ''' Chose which interaction takes place. '''
        xi = self.rng.random()
        sigfrac_cum = np.cumsum(self.sig_list)/self.sig_tot
        self.idx_event = np.searchsorted(sigfrac_cum, xi)
        self.event_list[self.idx_event].func(**self.event_list[self.idx_event].args)

        self.idx_event_list.append(self.idx_event)
        self.coord_list.append([self.coord[X], self.coord[Y], self.coord[Z]])
        self.time_list.append(self.time)
        
    def step(self):
        ''' Evolve the particle packet one step. '''
        self.get_event_list()
        self.move()
        self.interaction()

    @property
    def gam(self):
        ''' Lorentz factor '''
        return 1+self.ener*const.eV/(self.m_i*const.c**2)
    
    @property
    def vmag(self):
        ''' Particle velocity magnitide '''
        return np.sqrt(1 - 1/self.gam**2)*const.c  
    
    @property
    def n(self):
        ''' Total number density '''

