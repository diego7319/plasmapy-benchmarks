"""Functions from mathematics to benchmark with aerospeed velocity """

import astropy.units as u
from plasmapy.formulary import (
    plasma_dispersion_func,
    plasma_dispersion_func_deriv,
    Fermi_integral
)


class mathematics:
    """
    Benchmark that times the performance of funcions from
    mathematics/mathematics.py file
    """

    def setup(self):
        pass

    def time_plasma_dispersion_func(self):
        plasma_dispersion_func(-1.52+0.47j)

    def time_plasma_dispersion_func_deriv(self):
        plasma_dispersion_func_deriv(-1.52+0.47j)

    def time_Fermi_integral(self):
        Fermi_integral(1, 0)
