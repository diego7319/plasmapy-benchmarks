"""Functions from Physics to benchmark with aerospeed velocity """

import astropy.constants as const
import astropy.units as u
from numpy import pi
from plasmapy.formulary import parameters

#For dielectric class
from plasmapy.formulary import (
    cold_plasma_permittivity_SDP,
    cold_plasma_permittivity_LRP,
    permittivity_1D_Maxwellian
)

#For dimensionless class
from plasmapy.formulary.dimensionless import quantum_theta, beta


#For distribution class
from plasmapy.formulary import (
    Maxwellian_1D,
    Maxwellian_velocity_2D,
    Maxwellian_velocity_3D,
    Maxwellian_speed_1D,
    Maxwellian_speed_2D,
    Maxwellian_speed_3D,
    kappa_velocity_1D,
    kappa_velocity_3D
)

#For parameters class
from plasmapy.formulary.parameters import (
    mass_density,
    Alfven_speed,
    ion_sound_speed,
    thermal_speed,
    thermal_pressure,
    kappa_thermal_speed,
    Hall_parameter,
    gyrofrequency,
    gyroradius,
    plasma_frequency,
    Debye_length,
    Debye_number,
    inertial_length,
    magnetic_pressure,
    magnetic_energy_density,
    upper_hybrid_frequency,
    lower_hybrid_frequency,
)

#For quantum class
from plasmapy.formulary import(
    deBroglie_wavelength,
    thermal_deBroglie_wavelength,
    Fermi_energy,
    Thomas_Fermi_length,
    Wigner_Seitz_radius,
    chemical_potential,
)

#For relativity class
from plasmapy.formulary import Lorentz_factor


class dielectric:
    """
    Benchmark that times the performance of funcions from
    Physics/dielectric package
    """
    B = 2*u.T
    species = ['e', 'D+']
    n = [1e18*u.m**-3, 1e18*u.m**-3]
    omega = 3.7e9*(2*pi)*(u.rad/u.s)
    T = 30 * 11600 * u.K
    particle = 'Ne'
    z_mean = 8 * u.dimensionless_unscaled
    vTh = parameters.thermal_speed(T, particle, method="most_probable")
    kWave = omega / vTh
    n_permittivity_1D_Maxwellian = 1e18 * u.cm**-3


    def setup(self):
        pass

    def time_cold_plasma_permittivity_SDP(self):
        permittivity = S, D, P = cold_plasma_permittivity_SDP(
                                                             dielectric.B,
                                                             dielectric.species,
                                                             dielectric.n,
                                                             dielectric.omega)
    def time_cold_plasma_permittivity_LRP(self):
        permittivity = S, D, P = cold_plasma_permittivity_LRP(
                                                             dielectric.B,
                                                             dielectric.species,
                                                             dielectric.n,
                                                             dielectric.omega)

    def time_permittivity_1D_Maxwellian(self):
        permittivity_1D_Maxwellian(dielectric.omega,
                                   dielectric.kWave,
                                   dielectric.T,
                                   dielectric.n_permittivity_1D_Maxwellian,
                                   dielectric.particle,
                                   dielectric.z_mean)


class dimensionless:
    """
    Benchmark that times the performance of funcions from
    Physics.dimensionless
    """

    def setup(self):
        pass

    def time_quantum_theta(self):
        quantum_theta(1*u.eV, 1e20*u.m**-3)

    def time_beta(self):
        beta(1*u.eV, 1e20*u.m**-3, 1*u.T)


class distribution:
    """
    Benchmark that times the performance of functions from
    Physics.distribution
    """

    def setup(self):
        pass

    def time_Maxwellian_1D(self):
        Maxwellian_1D(v=1*u.m/u.s,
                     T=30000 * u.K,
                     particle='e',
                     v_drift=0 * u.m / u.s)

    def time_Maxwellian_velocity_2D(self):
        Maxwellian_velocity_2D(vx=1 * u.m / u.s,
                               vy=1 * u.m / u.s,
                               T=30000*u.K,
                               particle='e',
                               vx_drift=0 * u.m / u.s,
                               vy_drift=0 * u.m / u.s)

    def time_Maxwellian_velocity_3D(self):
        Maxwellian_velocity_3D(vx=1 * u.m / u.s,
                               vy=1 * u.m / u.s,
                               vz=1 * u.m / u.s,
                               T=30000 * u.K,
                               particle='e',
                               vx_drift=0 * u.m / u.s,
                               vy_drift=0 * u.m / u.s,
                               vz_drift=0 * u.m / u.s)

    def time_Maxwellian_speed_1D(self):
        Maxwellian_speed_1D(v=1 * u.m / u.s,
                            T=30000 * u.K,
                            particle='e',
                            v_drift=0 * u.m / u.s)

    def time_Maxwellian_speed_2D(self):
        Maxwellian_speed_2D(v=1 * u.m / u.s,
                            T=30000 * u.K,
                            particle='e',
                            v_drift=0 * u.m / u.s)

    def time_Maxwellian_speed_3D(self):
        Maxwellian_speed_3D(v=1 * u.m / u.s,
                            T=30000*u.K,
                            particle='e',
                            v_drift=0 * u.m / u.s)

    def time_kappa_velocity_1D(self):
        kappa_velocity_1D(v=1 * u.m / u.s,
                          T=30000*u.K,
                          kappa=4,
                          particle='e',
                          v_drift=0 * u.m / u.s)

    def time_kappa_velocity_3D(self):
        kappa_velocity_3D(vx=1 * u.m / u.s,
                          vy=1 * u.m / u.s,
                          vz=1 * u.m / u.s,
                          T=30000 * u.K,
                          kappa=4,
                          particle='e',
                          vx_drift=0 * u.m / u.s,
                          vy_drift=0 * u.m / u.s,
                          vz_drift=0 * u.m / u.s)


class parameters:
    """
    Benchmark that times the performance of functions from
    Physics.parameters.py
    """

    def setup(self):
        pass

    def time_mass_density(self):
        mass_density(1 * u.m ** -3,'p')


    def time_Alfven_speed(self):
        Alfven_speed(1 * u.T,
                    1e-8 * u.kg * u.m ** -3)

    def time_ion_sound_speed(self):
        ion_sound_speed(T_i=1.3232 * u.MK,
                        T_e=1.831 * u.MK,
                        ion='p',
                        gamma_e=1,
                        gamma_i=3)

    def time_thermal_speed(self):
        thermal_speed(1 * u.MK)

    def time_thermal_pressure(self):
        thermal_pressure(1*u.eV,
                         1e20/u.m**3)


    def time_kappa_thermal_speed(self):
        kappa_thermal_speed(5*u.eV,
                            4,
                            'p',
                            'mean_magnitude')

    def time_Hall_parameter(self):
        Hall_parameter(1e6 / u.m ** 3,5e6*u.K,0.1*u.T,'He-4 +1')

    def time_gyrofrequency(self):
        gyrofrequency(0.01*u.T,
                      particle='T+',
                      to_hz=True)

    def time_gyroradius(self):
        gyroradius(0.2*u.T,
                   particle='p+',
                   T_i=1e5*u.K)

    def time_plasma_frequency(self):
        plasma_frequency(1e19*u.m**-3,
                         particle='p',
                         to_hz=True)

    def time_Debye_length(self):
        Debye_length(5e6*u.K,
                     5e15*u.m**-3)

    def time_Debye_number(self):
        Debye_number(5e6*u.K,
                     5e9*u.cm**-3)

    def time_inertial_length(self):
        inertial_length(5 * u.m ** -3,
                        'He+')

    def time_magnetic_pressure(self):
        magnetic_pressure(0.1*u.T)

    def time_magnetic_energy_density(self):
        magnetic_energy_density(0.1*u.T)

    def time_upper_hybrid_frequency(self):
        upper_hybrid_frequency(0.2*u.T,
                               n_e=5e19*u.m**-3,
                               to_hz = True)

    def time_lower_hybrid_frequency(self):
        lower_hybrid_frequency(0.2*u.T,
                               n_i=5e19*u.m**-3,
                               ion='D+',
                               to_hz = True)


class quantum:
    """
    Benchmark that times the performance of functions from
    Physics.quantum
    """

    def setup(self):
        pass

    def deBroglie_wavelength(self):
        deBroglie_wavelength(V = 0 * u.m / u.s,
                             particle = 'D+')

    def time_thermal_deBroglie_wavelength(self):
        thermal_deBroglie_wavelength(1 * u.eV)

    def time_Fermi_energy(self):
        Fermi_energy(1e23 * u.cm**-3)

    def time_Thomas_Fermi_length(self):
        Thomas_Fermi_length(1e23 * u.cm**-3)

    def time_Wigner_Seitz_radius(self):
        Wigner_Seitz_radius(1e29 * u.m**-3)

    def time_chemical_potential(self):
        chemical_potential(n_e=1e21*u.cm**-3,
                           T=11000*u.K)


class relativity:
    """
    Benchmark that times the performance of functions from
    Physics.relativity
    """

    def setup(self):
        pass

    def time_Lorentz_factor(self):
        Lorentz_factor(1.4e8 * u.m / u.s)
