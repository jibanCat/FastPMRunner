"""
A class to automate the fastpm simulation
"""

class SimulationICs(object):
    """
    Class for creating the initial conditions for a fastpm simulation.

    There are a few things this class needs to do:
    - Generate linear theory input files (use python/make-pklin.py)
    - Generate fastpm .lua parameter file (ref: tests/nobodykit.lua)
    - Run fastpm simulation directly and generate power spectrum

    The class will store the parameters of the simulation.
    We also save a copy of the input and enough information to reproduce the
    results exactly in SimulationICs.json.
    Many things are left hard-coded.

    Init parameters:
    ----
    outdir     - Directory in which to save ICs
    box        - Box size in comoving Mpc/h
    npart      - Cube root of number of particles
    redshift   - redshift at which to generate ICs
    omegab     - baryon density.
    omegam     - Total matter density at z=0. (omega_m = omega_b + omega_cdm)
    hubble     - Hubble parameter, h, which is H0 / (100 km/s/Mpc)
    scalar_amp - A_s at k = 0.05, comparable to the Planck value.
    ns         - Scalar spectral index
    timesteps  - number of time steps for the simulation
    """
    def __init__(self, *,
            outdir: str = "nbodykit",   box: int = 384,  npart: int = 128,
            seed :         int   = 100,          redshift: float = 99,
            redend:        float = 0,            omega0:   float = 0.288, 
            omegab:        float = 0.0472,       hubble:   float = 0.7,
            scalar_amp:    float = 2.427e-9,     ns:       float = 0.97,
            fastpm_bin:    str = "fastpm",       timesteps: float = 10,
            python:        str = "python",
            cores:         int = 4) -> None:
        #Check that input is reasonable and set parameters
        #In Mpc/h
        assert box  < 20000
        self.box      = box

        #Cube root
        assert npart > 1 and npart < 16000
        self.npart    = int(npart)

        #Physically reasonable
        assert omega0 <= 1 and omega0 > 0
        self.omega0   = omega0

        assert omegab > 0 and omegab < 1
        self.omegab   = omegab

        assert redshift > 1 and redshift < 1100
        self.redshift = redshift


        assert redend >= 0 and redend < 1100
        self.redend = redend

        self.timesteps = timesteps

        assert hubble < 1 and hubble > 0
        self.hubble = hubble

        assert scalar_amp < 1e-7 and scalar_amp > 0
        self.scalar_amp = scalar_amp

        self.seed = seed