"""Functions from Atomic to benchmark with aerospeed velocity """

import astropy.constants as const
import astropy.units as u
from plasmapy.atomic import (
    known_isotopes,
    common_isotopes,
    stable_isotopes,
    isotopic_abundance,
    integer_charge,
    reduced_mass,
)


class Atomic:
    """
    Benchmark that times the performance of funcions from Atomic package
    """
    def setup(self):
        pass

    def time_known_isotopes(self):
        known_isotopes('H')
        known_isotopes('helium 1+')
        known_isotopes()[0:10]
        known_isotopes()[4:15]

    def time_common_isotopes(self):
        common_isotopes('H')
        common_isotopes('Fe')
        common_isotopes('Fe', most_common_only=True)
        common_isotopes()[0:7]

    def time_stable_isotopes(self):
        stable_isotopes('beryllium')
        stable_isotopes('Pb-209')
        stable_isotopes(118)
        stable_isotopes('U', unstable=True)[:5]

    def time_reduced_mass(self):
        reduced_mass('p+', 'e-')
        reduced_mass(5.4e-27 * u.kg, 8.6e-27 * u.kg)
        reduced_mass(6.4e-10 * u.kg, 8.6e-11 * u.kg)
        reduced_mass(1.4e-10 * u.kg, 2.6e-40 * u.kg)
