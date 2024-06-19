# system modules
import os
from types import SimpleNamespace

# math modules
import numpy as np
from scipy.integrate import trapz as integrate
from scipy.integrate import cumtrapz as cumintegrate

# custom modules
from functions import mag, interp
from helper import *
import const

# periodic table
elements = SimpleNamespace(**np.load('elements.npz'))

# beta particle data
data_EEDL = SimpleNamespace(**np.load('EEDL.npz', allow_pickle=True))

def calc_cos_th_scat(xi, ener, data):
    '''
    Compute the cosine of the scattering angle in an elastic scattering interaction.

    Args
    xi:   random number
    ener: kinetic energy [eV]
    data: a list of tuples (ener, cos_th, cos_th_pdf)
        ener:       energy [eV]
        cos_th:     cosine of the scattering angle
        cos_th_pdf: PDF of the cosine of the scattering angle

    Returns
    cos_th: scattering angle
    '''
    nener = len(data)
    ener_list, cos_th_list = np.zeros(nener), np.zeros(nener)
    
    for i, data_ener in enumerate(data):
    
        ener_list[i], cos_th, cos_th_pdf = data_ener
        cos_th_cdf = cumintegrate(cos_th_pdf, cos_th, initial=0) / integrate(cos_th_pdf, cos_th)
        cos_th_list[i] = interp(xi, cos_th_cdf, cos_th)

    cos_th = interp(ener, ener_list, cos_th_list, llim=-1, ulim=1)
    return cos_th

def calc_ener_loss(xi, ener, data):
    '''
    Compute the energy lost in an ionization interaction.

    Args
    xi:   random number
    ener: kinetic energy [eV]
    data: a list of tuples (ener, cos_th, cos_th_pdf)
        ener:          energy [eV]
        ener_loss:     energy lost [eV]
        ener_loss_pdf: PDF of the energy lost

    Returns
    ener_loss: energy lost [eV]
    '''
    nener = len(data)
    ener_list, ener_loss_list = np.zeros(nener), np.zeros(nener)

    for i, data_ener in enumerate(data):

        ener_list[i], ener_loss, ener_loss_pdf = data_ener
        ener_loss_cdf = cumintegrate(ener_loss_pdf, ener_loss, initial=0) / integrate(ener_loss_pdf, ener_loss)
        ener_loss_list[i] = interp(xi, ener_loss_cdf, ener_loss)

    ener_loss = interp(ener, ener_list, ener_loss_list)
    return ener_loss

def calc_moller(xi, ener, frac_lim=1e-3, num=1024):
    ''' 
    Compute the cross section, energy loss, and scattering angle in a Moller (electron-electron) scattering event.
    
    Args
    xi:       random number
    ener:     kinetic energy [eV]
    frac_lim: relative kinetic energy loss below which we ignore Moller scattering
    num:      number of points to sample the CDF

    Returns
    sig:       cross section [1/cm^2]
    cos_th:    scattering angle
    ener_loss: energy loss [eV]
    '''        
    sin_th_cm_lim = np.sqrt(1 - (1-2*frac_lim)**2)
    sin_th_cm = np.logspace(np.log10(sin_th_cm_lim), 0, num)[::-1]

    ener_tot_cm = np.sqrt(const.m_e*const.c**2/2 * (2*const.m_e*const.c**2 + ener*const.eV))
    fac1 = (ener_tot_cm**2-const.m_e**2*const.c**4)**2
    fac2 = (8*ener_tot_cm**4 - 4*const.m_e**2*const.c**4*ener_tot_cm**2 - const.m_e**4*const.c**8) / fac1
    fac3 = 4*(2*ener_tot_cm**2 - const.m_e**2*const.c**4)**2 / fac1
    fac4 = const.alpha**2*const.hbar**2*const.c**2/(4*ener_tot_cm**2)
    sig_cdf = 4*np.pi*fac4*(1-sin_th_cm + fac2*(1-1/sin_th_cm) - fac3/3*(1-1/sin_th_cm**3) )
    sig = sig_cdf[-1]
    sig_cdf /= sig

    sin_th_cm = interp(xi, sig_cdf, sin_th_cm)
    cos_th_cm = np.sqrt(1-sin_th_cm**2)
    ener_loss = ener/2*(1-cos_th_cm)
    cos_th = np.sqrt((2*const.m_e*const.c**2+ener*const.eV)*(1+cos_th_cm)/(4*const.m_e*const.c**2+ener*const.eV*(1+cos_th_cm)))

    return sig, cos_th, ener_loss

def calc_sig_Bturb(Bmag, fturb, ener, rho, m_i, q_i, q=5/3, Lmax=1e16):
    ''' 
    Compute the effective cross section of the particle to a turbulent magnetic field. 
    
    Args
    Bmag:  total magnetic field amplitude [G]
    fturb: ratio of turbulent to total magnetic energy
    ener:  kinetic energy [eV]
    rho:   density [g/cm^3]
    m_i:   particle mass [g]
    q_i:   particle charge [esu]
    q:     power law exponent of magnetic turbulence spectrum
    Lmax:  largest scale of magnetic turbulence [cm]

    Returns
    sig_Bturb: effective cross section [1/cm^2]
    '''
    gam = 1 + ener*const.eV/(m_i*const.c**2)
    beta = np.sqrt(1 - 1/gam**2)
    beta_A = Bmag / np.sqrt(4*np.pi*rho*const.c**2) # Alfven speed relative to the speed of light
    func_beta_A = (1-(beta_A/beta)**(2-q))/(2-q) - (1-(beta_A/beta)**(4-q))/(4-q)
    
    Om = q_i*Bmag/(m_i*const.c) # particle gyro-frequency
    kmin = const.c/(Om*Lmax)    # minimum wavenumber of magnetic turbulence spectrum
    lam_Bturb = const.c*beta**(2-q)*gam**(2-q)*func_beta_A * 2/(np.pi*(q-1)*fturb*const.c*kmin) * (const.c*kmin/Om)**(2-q)
    sig_Bturb = 1/(lam_Bturb*rho*const.N_A)

    return sig_Bturb

class Sim(object):
    '''
    Monte Carlo charged particle transport simulation.
    
    Args
    
    ener: partlce energy [eV]
    elem: background element
    rho:  background density [g/cm^3]

    m_i: particle mass [g]
    q_i: particle charge [esu]

    do_iso:    use isotropic scattering
    do_ion:    include ionization of background ions
    do_exc:    include excitation of background ions
    do_brem:   include Bremsstrahlung with background ions
    do_scat:   include elastic scattering with background ions
    do_moller: include moller scattering with background electrons

    seed: random number generator seed
    '''
    def __init__(self, ener=300000, rho=1e-16, m_i=const.m_e, q_i=const.e, do_iso=False, seed=None, do_ion=True, do_exc=True, do_brem=True, do_scat=True, do_moller=True, do_sync=True):
        
        self.ener_init = ener
        self.rho       = rho
        
        self.m_i  = m_i
        self.q_i  = q_i

        self.do_Bfield = False
        self.Bmag_co   = 0.
        self.Bhat_co   = np.array([0., 0., 0.])
        self.Bmag_turb = 0.
        self.do_iso    = do_iso

        self.do_ion    = do_ion
        self.do_exc    = do_exc
        self.do_brem   = do_brem
        self.do_scat   = do_scat
        self.do_moller = do_moller
        self.do_sync   = do_sync
        self.spec_list = []
        
        self.seed = seed
        self.rng  = np.random.default_rng(self.seed)
        
        self.reset()

    def reset(self, ener=None):
        ''' Reset the simulation. '''
        self.ener_init = self.ener_init if ener == None else ener
        self.ener = self.ener_init if ener == None else ener
        self.time  = 0
        self.coord = np.zeros(3)
        self.vhat  = self.rand_dir()  
        self.ev_list = []

    def add_Bturb(self, Bmag, q=5/3, Lmax=1e16):
        '''
        Add a turbulent magnetic field.

        Bturb: turbulent magnetic field amplitude [G]
        q:     power law exponent of magnetic turbulence spectrum
        Lmax:  largest scale of magnetic turbulence [cm]
        '''
        self.do_Bfield = True
        self.Bmag_turb = Bmag
        self.q         = q
        self.Lmax      = Lmax
        self.new_Bturb()
    
    def add_Bco(self, Bmag, Bhat):
        ''' 
        Add a coherent magnetic field. 
        
        Args
        Bmag: coherent magnetic field amplitude [G]
        Bhat: coherent magnetic field direction
        '''
        self.do_Bfield = True
        self.Bmag_co   = Bmag
        self.Bhat_co   = Bhat
        self.new_Bturb()

    def add_elec(self, x_e):
        '''
        Add electrons to the background plasma.

        Args
        x_e: elecron fraction
        '''
        n_ion = 0
        for spec_data in self.spec_list:
            n_ion += self.rho*const.N_A*spec_data.ab
        self.n_e = x_e*n_ion / (1-x_e)
    
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
                    sig = sig*ab # effective cross section
                    
                    spec_data.sig_ion_name_list.append(dtype_name)
                    spec_data.sig_ion_x_list.append(ener)
                    spec_data.sig_ion_y_list.append(sig)
                    spec_data.ener_bind_list.append(ener_bind)
                
            elif dtype_name[:3] == 'sig':

                    ener, sig = data
                    sig = sig*ab # effective cross section
                    
                    setattr(spec_data, '%s_x' % dtype_name, ener)
                    setattr(spec_data, '%s_y' % dtype_name, sig)

            elif dtype_name == 'th_scat':
        
                spec_data.data_th_scat = data

            elif dtype_name[:8] == 'spec_ion':

                spec_data.data_ener_loss_list.append(data)

            elif dtype_name == 'spec_brem':

                spec_data.ener_loss_brem_x, spec_data.ener_loss_brem_y = data

            elif dtype_name == 'spec_exc':

                spec_data.ener_loss_exc_x, spec_data.ener_loss_exc_y = data

        self.spec_list.append(spec_data)

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
                    ener_loss_ion = calc_ener_loss(self.rng.random(), self.ener, spec_data.data_ener_loss_list[i])
                    sig_ion = interp(self.ener, spec_data.sig_ion_x_list[i], spec_data.sig_ion_y_list[i])
                    self.sig_list.append(sig_ion)
                    self.event_list.append(SimpleNamespace(
                        func = self.dep_ener,
                        args = dict(ener_loss=ener_loss_ion),
                        Z = spec_data.Z,
                        cat = 'ionization',
                        label = '%s ionization (%s)' % (spec_data.symbol, sig_ion_name[8:]),
                        data = [ener_loss_ion, ener_loss_ion-spec_data.ener_bind_list[i]]
                    ))

            if self.do_exc:
                
                ener_loss_exc = interp(self.ener, spec_data.ener_loss_exc_x, spec_data.ener_loss_exc_y)
                sig_exc = interp(self.ener, spec_data.sig_exc_x, spec_data.sig_exc_y)
                self.sig_list.append(sig_exc)
                self.event_list.append(SimpleNamespace(
                    func = self.dep_ener,
                    args = dict(ener_loss=ener_loss_exc),
                    Z = spec_data.Z,
                    cat = 'excitation',
                    label = '%s excitation' % (spec_data.symbol),
                    data = [ener_loss_exc, 0]
                ))

            if self.do_brem:

                ener_loss_brem = interp(self.ener, spec_data.ener_loss_brem_x, spec_data.ener_loss_brem_y)
                sig_brem = interp(self.ener, spec_data.sig_brem_x, spec_data.sig_brem_y)
                self.sig_list.append(sig_brem)
                self.event_list.append(SimpleNamespace(
                    func = self.dep_ener,
                    args = dict(ener_loss=ener_loss_brem),
                    Z = spec_data.Z,
                    cat = 'Bremsstrahlung',
                    label = '%s Bremsstrahlung' % (spec_data.symbol),
                    data = [ener_loss_brem, 0]
                ))

            if self.do_scat:

                cos_th_scat = calc_cos_th_scat(self.rng.random(), self.ener, spec_data.data_th_scat)
                sig_scat_la = interp(self.ener, spec_data.sig_scat_la_x, spec_data.sig_scat_la_y)
                self.sig_list.append(sig_scat_la)
                self.event_list.append(SimpleNamespace(
                    func = self.scat,
                    args = dict(cos_th=cos_th_scat),
                    Z = spec_data.Z,
                    cat = 'elastic scatter',
                    label = '%s elastic scatter (large angle)' % (spec_data.symbol),
                    data = [cos_th_scat, 0]
                ))

                sig_scat_sa = interp(self.ener, spec_data.sig_scat_x, spec_data.sig_scat_y) - sig_scat_la
                self.sig_list.append(sig_scat_sa)
                self.event_list.append(SimpleNamespace(
                    func = self.dep_ener,
                    args = dict(ener_loss=0),
                    Z = spec_data.Z,
                    cat = 'elastic scatter',
                    label = '%s elastic scatter (small angle)' % (spec_data.symbol),
                    data = [0, 0]
                ))

        if self.do_moller:

            sig_moller, cos_th_scat, ener_loss = calc_moller(self.rng.random(), self.ener)
            sig_moller *= self.n_e/(self.rho*const.N_A) # effective cross section
            self.sig_list.append(sig_moller)
            self.event_list.append(SimpleNamespace(
                func = [self.scat, self.dep_ener],
                args = [dict(cos_th=cos_th_scat), dict(ener_loss=ener_loss)],
                Z = None,
                cat = 'Moller scatter',
                label = 'Moller scatter',
                data = [cos_th_scat, ener_loss]
            ))

        if self.Bmag_turb>0.:
            
            fturb = (self.Bmag_turb/self.Bmag)**2
            sig_Bturb = calc_sig_Bturb(self.Bmag, fturb, self.ener, self.rho, self.m_i, self.q_i, q=self.q, Lmax=self.Lmax)
            self.sig_list.append(sig_Bturb)
            self.event_list.append(SimpleNamespace(
                func = self.new_Bturb,
                args = dict(),
                Z = None,
                cat = 'Bturb scatter',
                label = 'Bturb scatter',
                data = [0, 0]
            ))

        self.sig_tot = np.sum(self.sig_list)

    def new_Bturb(self):
        ''' Resample the direction of the turbulent magnetic field. '''
        Bvec = self.rand_dir()*self.Bmag_turb + self.Bhat_co*self.Bmag_co
        self.Bmag = mag(Bvec)
        self.Bhat = Bvec/self.Bmag
    
    def dep_ener(self, ener_loss):
        ''' Deposit energy. '''
        self.ener = self.ener - ener_loss
    
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
            
            vhatnew = shat*cos_th + np.cross(vperp, shat)*sin_th + vperp*np.sum(vperp*shat)*(1-cos_th)
            vhatnew = vhatnew*np.cos(phi) + np.cross(shat, vhatnew)*np.sin(phi) + shat*np.sum(shat*vhatnew)*(1-np.cos(phi))
            self.vhat = vhatnew/mag(vhatnew)
    
    def move(self):
        ''' Move the packet until the next interaction. '''
        dis = -np.log(1.0-self.rng.random()) / (self.rho*const.N_A*self.sig_tot) # sample distance to collision by inversion technique
        self.time += dis/self.vmag

        if self.do_Bfield:
            
            self.coord += dis*np.sum(self.vhat*self.Bhat)*self.Bhat
            phi = 2.0*np.pi*self.rng.random()
            vhatnew = self.vhat*np.cos(phi) + np.cross(self.Bhat, self.vhat)*np.sin(phi) + self.Bhat*np.sum(self.Bhat*self.vhat)*(1-np.cos(phi))
            self.vhat = vhatnew/mag(vhatnew)

            if self.do_sync:
                sin_alpha = np.sqrt(1 - np.sum(self.vhat*self.Bhat)**2)
                beta = self.vmag/const.c
                ener_loss = 2/3 * self.q_i**4*self.gam**2*beta**2*self.Bmag**2*sin_alpha**2 / (self.m_i**2*const.c**5) * dis/self.vmag
                self.dep_ener(ener_loss)

        else:

            self.coord += dis*self.vhat

    def interaction(self):
        ''' Chose which interaction takes place. '''
        xi = self.rng.random()
        sigfrac_cum = np.cumsum(self.sig_list)/self.sig_tot
        self.idx_event = np.searchsorted(sigfrac_cum, xi)
        event = self.event_list[self.idx_event]
        if type(event.func) == list:
            for func, args in zip(event.func, event.args):
                func(**args)
        else: event.func(**event.args)

        self.ev_list.append(SimpleNamespace(label=event.label, cat=event.cat, Z=event.Z, time=self.time, x=self.coord[X], y=self.coord[Y], z=self.coord[Z], ener=self.ener, data1=event.data[0], data2=event.data[1]))
        
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

