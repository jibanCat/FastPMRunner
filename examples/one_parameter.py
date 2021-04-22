"""
An example code for varying only one parameter
and return power spectrum
"""

import os
import time
import numpy as np
from matplotlib import pyplot as plt

from FastPMRunner.simulationic import SimulationICs

def one_parameter_hubble(fastpm_bin: str = "fastpm"):
    """
    Vary one parameter: hubble parameter, h, which is H0 / (100 km/s/Mpc)

    Parameters:
    ----
    fastpm_bin: your fastpm binary
    """
    # hard coded output simulation folder
    folder = "simulation_files/"
    if not os.path.exists(folder):
        os.mkdir(folder)

    # simulation sub folder: containing files from each simulation
    subfolder_fn = lambda number : os.path.join(folder, "hubble_{:04d}".format(number))

    # hard coded sampling range, 0.65 ~ 0.75, 10 samples
    uniform_samples = np.linspace(0.65, 0.75, num=10)

    for i,hubble in enumerate(uniform_samples):
        # only vary hubble. others fixed
        sim = SimulationICs(outdir=subfolder_fn(i), hubble=hubble, fastpm_bin=fastpm_bin)

        # excute the simulation
        tic = time.time()
        print("Performing simulation {} ... ".format(i), end="")
        sim.make_simulation()

        print("took {:.3g} seconds".format(time.time() - tic))

        # plot a figure to show one parameter variation of the matter power spectrum
        # -1: the final redshift, z = 0, the current Universe
        plt.loglog(sim.kk[-1], sim.powerspecs[-1], label="hubble = {:.3g}".format(hubble))

    plt.xlabel("k (h/Mpc)")
    plt.ylabel("P(k)")
    plt.legend()
    plt.savefig(os.path.join(folder, "pk_one_param_hubble.pdf"), format="pdf", dpi=300)
