"""Functions from Physics to benchmark with aerospeed velocity """

import astropy.constants as const
import astropy.units as u
from numpy import pi
from plasmapy.physics import (
    cold_plasma_permittivity_SDP,
    cold_plasma_permittivity_LRP
)

class Dielectric:
    """
    Benchmark that times the performance of funcions from
    Physics/dielectric package
    """
    B = 2*u.T
    species = ['e', 'D+']
    n = [1e18*u.m**-3, 1e18*u.m**-3]
    omega = 3.7e9*(2*pi)*(u.rad/u.s)
    def setup(self):
        pass

    def time_cold_plasma_permittivity_SDP(self):
        permittivity = S, D, P = cold_plasma_permittivity_SDP(
                                                             Dielectric.B,
                                                             Dielectric.species,
                                                             Dielectric.n,
                                                             Dielectric.omega)
    def time_cold_plasma_permittivity_LRP(self):
        permittivity = S, D, P = cold_plasma_permittivity_LRP(
                                                             Dielectric.B,
                                                             Dielectric.species,
                                                             Dielectric.n,
                                                             Dielectric.omega)
