# system modules
import os
from types import SimpleNamespace

# math modules
import numpy as np
from scipy.integrate import trapz as integrate
from scipy.integrate import cumtrapz as cumintegrate
from scipy.interpolate import interp1d

# custom modules
from functions import mag, interp
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

    do_ion:   include ionization of background ions
    do_exc:   include excitation of background ions
    do_brem:  include Bremsstrahlung with background ions
    do_scat:  include elastic scattering with background ions
    '''
    def __init__(self, ener=300000, rho=1e-16, m_i=const.m_e, q_i=const.e, do_iso=False, seed=None, do_ion=True, do_exc=True, do_brem=True, do_scat=True):
        
        self.ener_init = ener
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
        
        self.reset()

    def reset(self):
        ''' Reset the simulation. '''
        self.ener = self.ener_init
        self.time  = 0
        self.coord = np.zeros(3)
        self.vhat  = self.rand_dir()  

        self.ev_label_list = []
        self.ev_cat_list   = []
        self.ev_Z_list     = []
        self.time_list     = []
        self.coord_list    = []
        self.ener_list     = []

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
        
        spec_data.sig_ion_name_list   = []
        spec_data.sig_ion_x_list      = []
        spec_data.sig_ion_y_list      = []
        spec_data.ener_bind_list      = []
        spec_data.data_ener_loss_list = []

        for i, dtype_name in enumerate(data_EEDL.dtype_list):

            data = data_EEDL.data[Z-1, i]
            if type(data) == type(None): continue
            
            if dtype_name[:7] == 'sig_ion' and len(dtype_name) > 7:

                    ener, sig, ener_bind = data
                    sig = sig*ab/A # effective cross section
                    
                    spec_data.sig_ion_name_list.append(dtype_name)
                    spec_data.sig_ion_x_list.append(ener)
                    spec_data.sig_ion_y_list.append(sig)
                    spec_data.ener_bind_list.append(ener_bind)
                
            elif dtype_name[:3] == 'sig':

                    ener, sig = data
                    sig = sig*ab/A # effective cross section
                    
                    setattr(spec_data, '%s_x' % dtype_name, ener)
                    setattr(spec_data, '%s_y' % dtype_name, sig)

            elif dtype_name == 'th_scat':
        
                spec_data.data_th_scat = data

            elif dtype_name[:8] == 'spec_ion' and len(dtype_name) > 7:

                spec_data.data_ener_loss_list.append(data)

            elif dtype_name == 'spec_brem':

                spec_data.ener_loss_brem_x, spec_data.ener_loss_brem_y = data

            elif dtype_name == 'spec_exc':

                spec_data.ener_loss_exc_x, spec_data.ener_loss_exc_y = data

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
                
                    if self.ener < spec_data.ener_bind_list[i]: continue
                    ener_loss = self.calc_ener_loss(spec_data.data_ener_loss_list[i])
                    sig_ion = interp(self.ener, spec_data.sig_ion_x_list[i], spec_data.sig_ion_y_list[i], logx=True, logy=True)
                    self.sig_list.append(sig_ion)
                    self.event_list.append(SimpleNamespace(
                        func = self.dep_ener,
                        args = dict(ener_loss=ener_loss),
                        Z = spec_data.Z,
                        cat = 'ionization',
                        label = '%s ionization (%s)' % (spec_data.symbol, sig_ion_name[8:])
                    ))

            if self.do_exc:

                ener_loss_exc = interp(self.ener, spec_data.ener_loss_exc_x, spec_data.ener_loss_exc_y, logx=True, logy=True)
                sig_exc = interp(self.ener, spec_data.sig_exc_x, spec_data.sig_exc_y, logx=True, logy=True)
                self.sig_list.append(sig_exc)
                self.event_list.append(SimpleNamespace(
                    func = self.dep_ener,
                    args = dict(ener_loss=ener_loss_exc),
                    Z = spec_data.Z,
                    cat = 'excitation',
                    label = '%s excitation' % (spec_data.symbol)
                ))

            if self.do_brem:

                ener_loss_brem = interp(self.ener, spec_data.ener_loss_brem_x, spec_data.ener_loss_brem_y, logx=True, logy=True)
                sig_brem = interp(self.ener, spec_data.sig_brem_x, spec_data.sig_brem_y, logx=True, logy=True)
                self.sig_list.append(sig_brem)
                self.event_list.append(SimpleNamespace(
                    func = self.dep_ener,
                    args = dict(ener_loss=ener_loss_brem),
                    Z = spec_data.Z,
                    cat = 'Bremsstrahlung',
                    label = '%s Bremsstrahlung' % (spec_data.symbol)
                ))

            if self.do_scat:

                cos_th_scat = self.calc_cos_th_scat(spec_data.data_th_scat)
                sig_scat_la = interp(self.ener, spec_data.sig_scat_la_x, spec_data.sig_scat_la_y, logx=True, logy=True)
                self.sig_list.append(sig_scat_la)
                self.event_list.append(SimpleNamespace(
                    func = self.scat,
                    args = dict(cos_th=cos_th_scat),
                    Z = spec_data.Z,
                    cat = 'elastic scatter',
                    label = '%s elastic scatter (large angle)' % (spec_data.symbol)
                ))

                sig_scat_sa = interp(self.ener, spec_data.sig_scat_x, spec_data.sig_scat_y, logx=True, logy=True) - sig_scat_la
                self.sig_list.append(sig_scat_sa)
                self.event_list.append(SimpleNamespace(
                    func = self.dep_ener,
                    args = dict(ener_loss=0),
                    Z = spec_data.Z,
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
                Z = None,
                cat = 'Bturb scatter',
                label = 'Bturb scatter'
            ))

        self.sig_tot = np.sum(self.sig_list)

    def calc_ener_loss(self, data):
        '''
        Compute the energy lost in an ionization interaction.

        Args
        data: A list of tuples (ener, cos_th, cos_th_pdf)
            ener:          energy [eV]
            ener_loss:     energy lost [eV]
            ener_loss_pdf: PDF of the energy lost
        '''
        xi = np.random.random() # compute a random number

        nener = len(data)
        ener_list, ener_loss_list = np.zeros(nener), np.zeros(nener)

        for i, data_ener in enumerate(data):

            ener_list[i], ener_loss, ener_loss_pdf = data_ener
            ener_loss_cdf = cumintegrate(ener_loss_pdf, ener_loss, initial=0) / integrate(ener_loss_pdf, ener_loss)
            ener_loss_list[i] = interp(xi, ener_loss_cdf, ener_loss, logx=True, logy=True)

        ener_loss = interp(self.ener, ener_list, ener_loss_list, logx=True, logy=True)
        return ener_loss

    def dep_ener(self, ener_loss):
        ''' Deposit energy. '''
        self.ener = self.ener - ener_loss

    def calc_cos_th_scat(self, data):
        '''
        Compute the cosine of the scattering angle.

        Args
        data: A list of tuples (ener, cos_th, cos_th_pdf)
            ener:       energy [eV]
            cos_th:     cosine of the scattering angle
            cos_th_pdf: PDF of the cosine of the scattering angle
        '''
        xi = self.rng.random() # compute a random number

        nener = len(data)
        ener_list, cos_th_list = np.zeros(nener), np.zeros(nener)
        
        for i, data_ener in enumerate(data):
        
            ener_list[i], cos_th, cos_th_pdf = data_ener
            cos_th_cdf = cumintegrate(cos_th_pdf, cos_th, initial=0) / integrate(cos_th_pdf, cos_th)
            cos_th_list[i] = interp(xi, cos_th_cdf, cos_th)

        cos_th = interp(self.ener, ener_list, cos_th_list, logx=True)
        cos_th = max(min(cos_th, 1), -1)
        return cos_th
    
    def scat(self, cos_th=None, shat=None):
        '''
        Scatter the particle packet:
            1. Compute a unit vector perpendicular to the scattering direction
            2. Rotate the velocity about the perpendicular vector by the scattering angle using the Rodrigues formula
            3. Rotate the velocity about the scattering direction by a random angle using the Rodrigues formula

        cos_th: scattering angle
        shat:   unit vector representing the direction off which to scatter
        '''
        if cos_th == None:
            self.vhat = self.rand_dir()
        else:
            if np.all(shat) == None: shat = self.vhat

            phi = 2.0*np.pi*self.rng.random()
            sin_th = np.sqrt(1-cos_th**2)
            
            vperp = np.zeros(3)
            vperp[X], vperp[Y] = shat[Y], -shat[X]
            vperp /= np.sqrt(shat[X]**2 + shat[Y]**2)
            
            vhatnew = shat * cos_th + np.cross(vperp, shat) * sin_th + vperp * np.sum(vperp*shat) * (1-cos_th)
            vhatnew = vhatnew * np.cos(phi) + np.cross(shat, vhatnew) * np.sin(phi) + shat * np.sum(shat*vhatnew) * (1-np.cos(phi))
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
        event = self.event_list[self.idx_event]
        event.func(**self.event_list[self.idx_event].args)

        self.ev_label_list.append(event.label)
        self.ev_cat_list.append(event.cat)
        self.ev_Z_list.append(event.Z)
        self.coord_list.append([self.coord[X], self.coord[Y], self.coord[Z]])
        self.time_list.append(self.time)
        self.ener_list.append(self.ener)
        
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

