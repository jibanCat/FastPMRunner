from nbodykit.lab import *
from nbodykit.cosmology import WMAP9, LinearPower, Cosmology
import numpy
from matplotlib import pyplot as plt


def save_powerspec(
    omega0: float = 0.288,
    omegab: float = 0.0472,
    hubble: float = 0.7,
    scalar_amp: float = 2.427e-9,
    ns: float = 0.97,
    outfile: str = "powerspec.txt",
):
    """Generate linear powerspec and save into file"""
    omegacdm = (
        omega0 - omegab
    )  # dark matter density = total mass density - baryon density

    MYPlanck = Cosmology(
        m_ncdm=[],  # no massive neutrino
        Omega0_b=omegab,
        Omega0_cdm=omegacdm,
        h=hubble,
        n_s=ns,
        A_s=scalar_amp,
    )  # \
    # .match(sigma8=0.8159)
    pklin0 = LinearPower(MYPlanck, redshift=0.0)

    k = numpy.logspace(-3, 2, 10000, endpoint=True)

    numpy.savetxt(outfile, list(zip(k, pklin0(k))))


def plot_pk(outfile: str):
    kk, pk = numpy.loadtxt(outfile).T
    plt.loglog(kk, pk)
    plt.xlabel("k")
    plt.ylabel("P(k)")
