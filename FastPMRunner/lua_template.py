"""
A simple lua template to run fastpm with 5 cosmological parameters
"""
from typing import List


def simple_lua_string(
    box: int,
    npart: int,
    seed: int,
    omega0: float,
    omegab: float,
    hubble: float,
    scalar_amp: float,
    ns: float,
    time_start: float,
    time_end: float,
    timesteps: float,
    read_powerspectrum: str = "my_pk_linear.txt",
    write_powerspectrum: str = "nbodykit/powerspec",
    write_runpb_snapshot: bool = False,
    write_snapshot: bool = False,
    write_fof: bool = False,
):
    """
    Return a simple lua parameter file with 5 cosmologies
    """

    lua_template = '''-- parameter file
------ Size of the simulation -------- 

-- For Testing
nc = {nc}
boxsize = {boxsize}

-------- Time Sequence ----
-- linspace: Uniform time steps in a
time_step = linspace({time_start}, {time_end}, {timesteps})

output_redshifts= {{9.0, 2.0, 1.0, 0.0}}  -- redshifts of output

-- Cosmology --
omega_m = {omega_m}
h       = {h}
scalar_amp = {scalar_amp}
scalar_spectral_index = {scalar_spectral_index}

-- Start with a linear density field
-- Power spectrum of the linear density field: k P(k) in Mpc/h units
-- Must be compatible with the Cosmology parameter
read_powerspectrum= "{read_powerspectrum}"
linear_density_redshift = 0.0 -- the redshift of the linear density field.
random_seed= {seed}
particle_fraction = 1.0
--
-------- Approximation Method ---------------
force_mode = "fastpm"
-- kernel_type = "1_4"

growth_mode = "LCDM"

pm_nc_factor = 2
lpt_nc_factor = 1

np_alloc_factor= 4.0      -- Amount of memory allocated for particle

-------- Output ---------------
-- 1d power spectrum (raw), without shotnoise correction
write_powerspectrum = "{write_powerspectrum}"'''.format(
        nc=npart,
        boxsize=box,
        time_start=time_start,
        time_end=time_end,
        timesteps=timesteps,
        omega_m=omega0,
        h=hubble,
        scalar_amp=scalar_amp,
        scalar_spectral_index=ns,
        read_powerspectrum=read_powerspectrum,
        seed=seed,
        write_powerspectrum=write_powerspectrum,
    )

    if write_runpb_snapshot:
        lua_template += '''-- Dark matter particle outputs (all particles)
write_runpb_snapshot= "nbodykit/tpm" '''
    if write_snapshot:
        lua_template += '''write_snapshot = "nbodykit/fastpm"'''
    if write_fof:
        lua_template += '''write_fof = "nbodykit/fastpm"'''

    return lua_template
