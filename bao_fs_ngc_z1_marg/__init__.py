"""A likelihood function using the FastPMRunner. Modified from the code at
https://github.com/Michalychforever/lss_montepython using this paper:
https://arxiv.org/abs/1909.05277

This code compares a measured 2-point fourier-space correlation function
(Multipoles of a redshift space power spectrum) from the BOSS galaxy survey to a theory model.
For the equations see Section 3.2 of https://arxiv.org/abs/1909.05277

I have removed all the 'theory model' parameters discussed at length
as the point of this is to explore numerical approximations.

I have also simplified this by:
1. Removing the Alcock-Paczynski test which is a bit complicated to implement, and does not dominate the constraints.
2. Removing the shot noise power which is close to zero
3. Only compute the monopole terms
"""

import os
import numpy as np
from numpy.fft import fft, ifft , rfft, irfft , fftfreq
from numpy import exp, log, log10, cos, sin, pi, cosh, sinh , sqrt
from scipy import interpolate
from scipy import special
from FastPMRunner.simulationic import SimulationICs

class Likelihood:
    """Likelihood function for the BOSS power spectrum."""
    def __init__(self):
        """Initialize the function, loading data and other useful functions that can be precomputed"""

        self.measurements_file = 'ngc_z1/pk.dat'
        self.covmat_file = 'ngc_z1/ngc_z1_unrec_mon+quad+alphas_cov_48bins.dat'
        self.window_file = 'ngc_z1/window.dat'

        # Specify k-range, redshift and nuisance parameters
        self.ksize = 48
        self.z      = 0.38

        ## LOAD IN DATA
        self.k = np.zeros(self.ksize,'float64')
        self.Pk0 = np.zeros(self.ksize,'float64')
        self.Pk2 = np.zeros(self.ksize,'float64')

        self.cov = np.zeros(
            (2*self.ksize+2, 2*self.ksize+2), 'float64')

        # Load covariance matrix
        datafile = open(self.covmat_file, 'r')
        for i in range(2*self.ksize+2):
            line = datafile.readline()
            while line.find('#') != -1:
                line = datafile.readline()
            for j in range(2*self.ksize+2):
                self.cov[i,j] = float(line.split()[j])
        datafile.close()

        # Load unreconstructed power spectrum
        datafile = open(self.measurements_file, 'r')
        for i in range(self.ksize):
            line = datafile.readline()
            while line.find('#') != -1:
                line = datafile.readline()
            self.k[i] = float(line.split()[0])
            self.Pk0[i] = float(line.split()[1])
            self.Pk2[i] = float(line.split()[2])
        datafile.close()

        ## Window function: only W0 is used.
        self.Nmax=128
        self.W0 = np.zeros((self.Nmax,1))
        self.W2 = np.zeros((self.Nmax,1))
        self.W4 = np.zeros((self.Nmax,1))
        datafile = open(self.window_file, 'r')
        for i in range(self.Nmax):
            line = datafile.readline()
            while line.find('#') != -1:
                line = datafile.readline()
            self.W0[i] = float(line.split()[0])
            self.W2[i] = float(line.split()[1])
            self.W4[i] = float(line.split()[2])
        datafile.close()

        # Precompute useful window function things
        kmax = 100.
        self.k0 = 5.e-4

        self.rmin = 0.01
        rmax = 1000.
        b = -1.1001
        bR = -2.001

        Delta = log(kmax/self.k0) / (self.Nmax - 1)
        Delta_r = log(rmax/self.rmin) / (self.Nmax - 1)
        i_arr = np.arange(self.Nmax)
        rtab = self.rmin * exp(Delta_r * i_arr)

        self.kbins3 = self.k0 * exp(Delta * i_arr)
        self.tmp_factor = exp(-1.*b*i_arr*Delta)
        self.tmp_factor2 = exp(-1.*bR*i_arr*Delta_r)[:,np.newaxis]

        jsNm = np.arange(-self.Nmax//2,self.Nmax//2+1,1)
        self.etam = b + 2*1j*pi*(jsNm)/self.Nmax/Delta

        def J_func(r,nu):
            gam = special.gamma(2+nu)
            r_pow = r**(-3.-1.*nu)
            sin_nu = np.sin(pi*nu/2.)
            J0 = -1.*sin_nu*r_pow*gam/(2.*pi**2.)
            return J0

        j0 = J_func(rtab.reshape(-1,1),self.etam.reshape(1,-1))
        self.J0_arr = j0[:,:,np.newaxis]

        self.etamR = bR + 2*1j*pi*(jsNm)/self.Nmax/Delta_r

        def Jk_func(k,nu):
            gam = special.gamma(2+nu)
            k_pow = k**(-3.-1.*nu)
            sin_nu = np.sin(pi*nu/2.)
            J0k = -1.*k_pow*gam*sin_nu*(4.*pi)
            return J0k

        j0k = Jk_func(self.kbins3.reshape(-1,1),self.etamR.reshape(1,-1))
        self.J0k_arr = j0k[:,:,np.newaxis]

        #Precompute inverse covariance matrix
        self.invcov = np.linalg.inv(self.cov)

        # Compute window response matrix
        self.response_matrix = np.zeros((self.ksize,self.Nmax))
        for i in range(self.Nmax):
            tmp_resp0 = self.window_response(i)
            self.response_matrix[:,i] = tmp_resp0[0]

        # Now convolve with window function
        self.invcovW = np.matmul(self.response_matrix.T,self.invcov)
        self.invcovWW = np.matmul(self.response_matrix.T,np.matmul(self.invcov,self.response_matrix))

    def window_response(self,k_index):
        """Window function. This is from 1607.03150 and includes only the monopole.
        Accurate for k > 0.0015."""
        Nmax = self.Nmax
        k0 = self.k0

        Pdiscrin0 = np.zeros(Nmax)
        Pdiscrin0[k_index] = 1

        cm0 = np.fft.fft(Pdiscrin0)/ Nmax
        cmsym0 = np.zeros(Nmax+1,dtype=np.complex_)

        all_i = np.arange(Nmax+1)
        f = (all_i+2-Nmax//2) < 1
        cmsym0[f] = k0**(-self.etam[f])*np.conjugate(cm0[-all_i[f]+Nmax//2])
        cmsym0[~f] = k0**(-self.etam[~f])*cm0[all_i[~f]-Nmax//2]

        cmsym0[-1] = cmsym0[-1] / 2
        cmsym0[0] = cmsym0[0] / 2

        xi0 = np.real(cmsym0*self.J0_arr).sum(axis=1)
        Xidiscrin0 = (xi0*self.W0)*self.tmp_factor2

        cmr0 = np.fft.fft(Xidiscrin0)/ Nmax

        cmsymr0 = np.zeros(Nmax+1,dtype=np.complex_)

        arr_i = np.arange(Nmax+1)
        f = (arr_i+2-Nmax//2)<1

        cmsymr0[f] = self.rmin**(-self.etamR[f])*np.conjugate(cmr0[-arr_i[f] + Nmax//2])
        cmsymr0[~f] = self.rmin**(-self.etamR[~f])* cmr0[arr_i[~f] - Nmax//2]

        cmsymr0[-1] = cmsymr0[-1] / 2
        cmsymr0[0] = cmsymr0[0] / 2

        P0t = np.real(cmsymr0*self.J0k_arr).sum(axis=1)
        P0int = interpolate.InterpolatedUnivariateSpline(self.kbins3,P0t)(self.k)
        return P0int

    def loglkl(self, params, box=384, npart=128, timesteps=10):
        """Compute the log-likelihood of the model, given the data and covariance loaded in __init__.
        The cosmology and model parameters are put in the params list.
        The fastpm accuracy parameters box, npart and timesteps are kw args.
        Model parameters:
        0 - h the hubble parameter
        1 - Omega_0 the total matter density
        2 - scalar_amp the primordial amplitude of fluctuations at k = 0.05 Mpc^{-1}
        3 - The linear bias, a variable which encodes how much galaxies differ from cold dark matter.
        """
        #Cosmology parameters
        hubble = params[0]
        omega0 = params[1]
        scalar_amp = params[2]
        #ns is not measured by BOSS very well, so fixed
        #ns = params[3]
        linear_bias = params[3]

        ## COMPUTE SPECTRA
        # Run fastPM for power spectrum with cosmological parameters.
        sim = SimulationICs(hubble = hubble,
                            omega0 = omega0,
                            scalar_amp = scalar_amp,
                            ns = 0.9649,
                            fastpm_bin = "fastpm",
                            box = box,
                            npart = npart,
                            timesteps = timesteps,
                            redend = self.z,
                            )
        sim.make_simulation()
        powerspec = sim.powerspecs[-1]
        kk = sim.kk[-1]
        #Galaxy power spectrum with multipoles, for reference.
        #Powergal = lambda mu, l: (linear_bias + ff * mu**2)**2 * sim.powerspecs[-1]
        #For higher moments we should do a multipole integral. Here just use l = 0 for simplicity.
        #Neglect the redshift space distortions because they are degenerate with linear bias
        #for the monopole.
        powergal = linear_bias**2 * powerspec
        powgalint = interpolate.InterpolatedUnivariateSpline(kk,powergal)

        # Now do a Fourier-transform to include the window function
        factor = np.exp(-1.*(self.kbins3*hubble/2.)**4.)*self.tmp_factor
        factint = interpolate.InterpolatedUnivariateSpline(self.kbins3*hubble, factor)
        Pdisc = powgalint(self.k)*factint(self.k)

        # Now compute chi^2: this is:
        #(Pth W - Pexp)^T C^-1 (Pth W - Pexp)
        chi2 = np.inner(Pdisc[:,0],np.inner(self.invcovWW,Pdisc[:,0]))
        chi2 += np.inner(self.Pk0,np.inner(self.invcov,self.Pk0))
        chi2 += -2.*np.inner(Pdisc[:,0],np.inner(self.invcovW,self.Pk0))

        # Compute log-likelihood
        loglkl = -0.5 * chi2
        return loglkl
