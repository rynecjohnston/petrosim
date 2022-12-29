"""
These are the helper functions, mostly differential equations used in the EC-AFC
model by Spera & Bohrson J. Petrol., 42, 999–1018, 2001.
"""

import functools

import numpy as np

from petrosim.models.ecafc import loader as ldr


@functools.lru_cache
def unnormalize_temp(Tnorm):
    """
    Unnormalize the temperature

    :param Tnorm: The normalized temperature
    :type Tnorm: float

    :return: The unnormalized temp
    :rtype: float
    """

    return ldr.params.Tlm * Tnorm - ldr.params.K


@functools.lru_cache
def tempK(Tnorm):
    """
    Find the temperature in Kelvin

    :param Tnorm: The normalized temperature
    :type Tnorm: float

    :return: The temperature in Kelvin
    :rtype: float
    """

    return ldr.params.Tlm * Tnorm


@functools.lru_cache
def Tm_slope(dT):
    """
    Scale the Tm slope by 1_000

    :param Tnorm: The normalized temperature
    :type Tnorm: float

    :return: The scaled Tm slope
    :rtype: float
    """

    return dT * 1e3


@functools.lru_cache
def phi(T, component):
    """
    Melt fraction 

    :param T: The component temperature
    :type T: float

    :param component: The component of interest
    :type component: str

    :return: The melt fraction
    :rtype: float
    """

    if component == ldr.ASSIMILANT:
        Tlx = ldr.params.Tla
    elif component == ldr.MAGMA:
        Tlx = ldr.params.Tlm
    return (T - ldr.params.Ts) / (Tlx - ldr.params.Ts)


@functools.lru_cache
def phi_prime(component):
    """
    Melt fraction T derivative

    :param component: The component of interest
    :type component: str

    :return: The melt fraction
    :rtype: float
    """

    if component == ldr.ASSIMILANT:
        Tlx = ldr.params.Tla
    elif component == ldr.MAGMA:
        Tlx = ldr.params.Tlm
    return ldr.params.Tlm / (Tlx - ldr.params.Ts)


@functools.lru_cache
def melt_productivity(T, component):
    """
    Melt productivity function (logistical form)

    :param component: The component of interest
    :type component: str

    :return: Melt productivity
    :rtype: float
    """

    loga_x, logb_x = 0, 0
    if component == ldr.ASSIMILANT:
        loga_x = ldr.params.loga_a
        logb_x = ldr.params.logb_a
    elif component == ldr.MAGMA:
        loga_x = ldr.params.loga_m
        logb_x = ldr.params.logb_m

    return 1 / (1 + (loga_x * np.exp(logb_x * phi(T, component))))


@functools.lru_cache
def melt_productivity_T_deriv(T, component):
    """
    Melt productivity function T derivative (logistical form)

    :param T: The component temperature
    :type T: float

    :param component: The component of interest
    :type component: str

    :return: Melt productivity T derivative
    :rtype: float
    """

    loga_x, logb_x = 0, 0
    if component == ldr.ASSIMILANT:
        loga_x = ldr.params.loga_a
        logb_x = ldr.params.logb_a
        Tlx = ldr.params.Tla

    elif component == ldr.MAGMA:
        loga_x = ldr.params.loga_m
        logb_x = ldr.params.logb_m
        Tlx = ldr.params.Tlm

    numer = -1 * loga_x * logb_x * np.exp(logb_x * phi(T, component))
    denom = (loga_x * np.exp(logb_x * phi(T, component)) + 1)**2
    fxprime = phi_prime(component) * numer / denom
    return fxprime


@functools.lru_cache
def mass_country_rock(T):
    """
    :param T: The temperature
    :type T: float

    :return: Mass of anatectic melt at given a temperature
    :rtype: float
    """

    temp = T / ldr.params.Tlm
    spec_nrg_m = ldr.params.cpm * (ldr.params.Tm0 - ldr.params.Tlm * temp)
    enth_cry_m = ldr.params.dhm * (1 - melt_productivity(T, ldr.MAGMA))
    spec_nrg_a = ldr.params.cpa * (ldr.params.Tlm * temp - ldr.params.Ta0)
    enth_fus_a = ldr.params.dha * melt_productivity(T, ldr.ASSIMILANT)
    return (spec_nrg_m + enth_cry_m) / (spec_nrg_a + enth_fus_a)


def mass_anatectic_melt(Ma0, fa):
    """
    :param Ma0: Mass of anatectic melt
    :type Ma0: float

    :param fa: Melt productivity of assimilant
    :type fa: float

    :return: Mass of the anatectic melt
    :rtype: float
    """

    return Ma0 * fa


def mass_cumulates(fm):
    """
    :param fa: Melt productivity of magma
    :type fa: float

    :return: Mass of cumulates
    :rtype: float
    """

    return 1 - fm


def mass_magma(Mastar, Mct):
    """
    :param Mastar: Mass of anatectic melt
    :type Mastar: float

    :param Mct: Mass of cumulates
    :type Mct: float

    :return: Mass of magma
    :rtype: float
    """

    return 1 + Mastar - Mct


def ratio_country_rock_to_cumulates(Ma0, Mct):
    """
    :param Ma0: Mass of country rock
    :type Ma0: float

    :param Mct: Mass of cumulates
    :type Mct: float

    :return: Ratio country rock to cumulates 
    :rtype: float
    """

    return Ma0 / Mct


def ratio_anatectic_to_cumulates(Mastar, Mct):
    """
    :param Mastar: Mass of anatectic melt
    :type Mastar: float

    :param Mct: Mass of cumulates
    :type Mct: float

    :return: Ratio of anatectic melt to cumulates
    :rtype: float
    """

    return Mastar / Mct


@functools.lru_cache
def conserv_energy_dTa_dTm(Tm, Ta, Ma0):
    """
    :param Tm: Normalized temperature of magma
    :type Tm: float

    :param Ta: Normalized temperature of assimilant
    :type Ta: float

    :param Ma0: Mass of country rock
    :type Ma0: float

    :return: Conservation of energy for restite-magma equilibrium
    :rtype: float
    """

    temp_m = tempK(Tm)
    temp_a = tempK(Ta)

    fm = melt_productivity(temp_m, ldr.MAGMA)
    fa = melt_productivity(temp_a, ldr.ASSIMILANT)
    fmprime = melt_productivity_T_deriv(temp_m, ldr.MAGMA)
    faprime = melt_productivity_T_deriv(temp_a, ldr.ASSIMILANT)

    term1 = -1 / Ma0
    term2_num = (ldr.params.Tlm * ldr.params.cpm) + (
        ldr.params.dhm * fmprime) + (Ma0 * ldr.params.cpa * ldr.params.Tlm *
                                     fa)
    term2_den = (
        ldr.params.Tlm * ldr.params.cpa *
        (1 - fa)) + (ldr.params.dha + ldr.params.cpa * ldr.params.Tlm *
                     (Tm - Ta)) * faprime
    return term1 * term2_num / term2_den


@functools.lru_cache
def conserv_mass_dMm_dTm(Tm, Ta, Ma0):
    """
    :param Tm: Normalized temperature of magma
    :type Tm: float

    :param Ta: Normalized temperature of assimilant
    :type Ta: float

    :param Ma0: Mass of country rock
    :type Ma0: float

    :return: Conservation of total mass for restite-magma equilibrium
    :rtype: float
    """

    temp_m = tempK(Tm)
    temp_a = tempK(Ta)

    fmprime = melt_productivity_T_deriv(temp_m, ldr.MAGMA)
    faprime = melt_productivity_T_deriv(temp_a, ldr.ASSIMILANT)
    dTa_dTm = conserv_energy_dTa_dTm(Tm, Ta, Ma0)
    return Ma0 * faprime * dTa_dTm + fmprime


@functools.lru_cache
def conc_trace_elem_dCm_dTm(Tm, Ta, Ma0, Mm, Cm, Ca0, Cm0, dTa_dTm):
    """
    :param Tm: Normalized temperature of magma
    :type Tm: float

    :param Ta: Normalized temperature of assimilant
    :type Ta: float

    :param Ma0: Mass of country rock
    :type Ma0: float

    :param Mm: Mass of magma
    :type Mm: float

    :param Cm: Concentration of trace element in magma at temperature Tm
    :type Cm: float

    :param Ca0: Concentration of trace element in assimilant initially
    :type Ca0: float

    :param Cm0: Concentration of trace element in magma initially
    :type Cm0: float

    :param dTa_dTm: Conservation of total mass for restite-magma equilibrium
    :type dTa_dTm: float

    :return: Variation of the concentration of trace element in the standing melt
    :rtype: float
    """

    temp_m = tempK(Tm)
    temp_a = tempK(Ta)

    fmprime = melt_productivity_T_deriv(temp_m, ldr.MAGMA)
    faprime = melt_productivity_T_deriv(temp_a, ldr.ASSIMILANT)
    Cabar = conc_trace_elem_in_anatectic_melt(Ta, Ca0) / Ca0
    Cmbar = Cm / Cm0
    s = Ca0 / Cm0
    Dm = distribution_coefficient(Tm, ldr.MAGMA)

    term1 = 1 / Mm
    term2 = Ma0 * (s * Cabar - Cmbar) * faprime * dTa_dTm
    term3 = Cmbar * (Dm - 1) * fmprime

    return term1 * (term2 + term3) * Cm0


@functools.lru_cache
def conc_trace_elem_in_anatectic_melt(Ta, Ca0):
    """
    :param Ta: Normalized temperature of assimilant
    :type Ta: float

    :param Ca0: Concentration of trace element in assimilant initially
    :type Ca0: float

    :return: Trace element conservation in the melt
    :rtype: float
    """

    temp_a = tempK(Ta)

    fa = melt_productivity(temp_a, ldr.ASSIMILANT)
    Da = distribution_coefficient(Ta, ldr.ASSIMILANT)
    return (Ca0 / Da) * ((1 - fa)**((1 - Da) / Da))


@functools.lru_cache
def isotop_ratio_dem_dTm(Tm, Ta, Ma0, Mm, Cm, Ca0, Cm0, dTa_dTm, ea0, em):
    """
    :param Tm: Normalized temperature of magma
    :type Tm: float

    :param Ta: Normalized temperature of assimilant
    :type Ta: float

    :param Ma0: Mass of country rock
    :type Ma0: float

    :param Mm: Mass of magma
    :type Mm: float

    :param Cm: Concentration of trace element in magma at temperature Tm
    :type Cm: float

    :param Ca0: Concentration of trace element in assimilant initially
    :type Ca0: float

    :param Cm0: Concentration of trace element in magma initially
    :type Cm0: float

    :param dTa_dTm: Conservation of total mass for restite-magma equilibrium
    :type dTa_dTm: float

    :param ea0: Isotopic ratio of trace element in assimilant initially
    :type ea0: float

    :param em: Isotopic ratio of trace element in standing melt
    :type em: float

    :return: Trace element isotope balance
    :rtype: float
    """

    temp_m = Tm * ldr.params.Tlm
    temp_a = Ta * ldr.params.Tlm

    fmprime = melt_productivity_T_deriv(temp_m, ldr.MAGMA)
    faprime = melt_productivity_T_deriv(temp_a, ldr.ASSIMILANT)
    Cabar = conc_trace_elem_in_anatectic_melt(Ta, Ca0) / Ca0
    Cmbar = Cm / Cm0
    s = Ca0 / Cm0

    term1 = 1 / Mm
    term2 = (s * (Cabar / Cmbar) * (ea0 - em) * Ma0 * faprime * dTa_dTm)
    return term1 * term2


@functools.lru_cache
def oxygen_isotope_comp_ddm_dTm(Tm, Ta, Ma0, Mm, dTa_dTm, dm):
    """
    :param Tm: Normalized temperature of magma
    :type Tm: float

    :param Ta: Normalized temperature of assimilant
    :type Ta: float

    :param Ma0: Mass of country rock
    :type Ma0: float

    :param Mm: Mass of magma
    :type Mm: float

    :param dTa_dTm: Conservation of total mass for restite-magma equilibrium
    :type dTa_dTm: float

    :param dm: Oxygen isotope composition in standing melt
    :type dm: float

    :return: Oxygen isotope balance
    :rtype: float
    """

    temp_m = Tm * ldr.params.Tlm
    temp_a = Ta * ldr.params.Tlm

    faprime = melt_productivity_T_deriv(temp_a, ldr.ASSIMILANT)

    term1 = 1 / Mm
    term2 = (ldr.params.da0 - dm) * ldr.params.koxy * Ma0 * faprime * dTa_dTm
    return term1 * term2


def distribution_coefficient(Tx, component):
    """
    :param Tx: Normalized temperature of magma
    :type Tm: float

    :param component: The component of interest
    :type component: str

    :return: Distribution coefficient
    :rtype: float
    """

    if component == ldr.ASSIMILANT:
        Dx0 = ldr.params.Da0
        dHx = ldr.params.dHa
    elif component == ldr.MAGMA:
        Dx0 = ldr.params.Dm0
        dHx = ldr.params.dHm
    return Dx0 * np.exp(-dHx / (ldr.params.R * ldr.params.Tlm * Tx))


class EquilibrationParams:
    """
    Perform a preliminary equilibration.
    """
    def __init__(self, Tnorm):
        """
        Calculate the various quantities at the given normalized equilibration
        temperature

        :param Tnorm: The normalized temperature
        :type Tnorm: float
        """

        self.Tnorm = Tnorm
        self.Teq = ldr.params.Tlm * Tnorm - ldr.params.K
        self.fm = melt_productivity(self.Teq + ldr.params.K, ldr.MAGMA)
        self.fa = melt_productivity(self.Teq + ldr.params.K, ldr.ASSIMILANT)
        self.Ma0 = mass_country_rock(self.Teq + ldr.params.K)
        self.Mastar = mass_anatectic_melt(self.Ma0, self.fa)
        self.Mct = mass_cumulates(self.fm)
        self.Mm = mass_magma(self.Mastar, self.Mct)
        self.Ma0_Mct = ratio_country_rock_to_cumulates(self.Ma0, self.Mct)
        self.Mastar_Mct = ratio_anatectic_to_cumulates(self.Mastar, self.Mct)

    def printParams(self, header=False):
        """
        Print the equilibrated parameters

        :param header: Whether to print the header
        :type header: bool
        """

        data_len = {
            'Tnorm': 5,
            'Teq': 8,
            'fa': 5,
            'fm': 5,
            'Ma0': 5,
            'Mastar': 6,
            'Mct': 5,
            'Mm': 5,
        }
        if header:
            print_line = '  '.join(
                [f'{data:>{length}s}' for data, length in data_len.items()])
            print(print_line)
        print_line = '  '.join([
            f'{getattr(self, data):>{length}.3f}'
            for data, length in data_len.items()
        ])
        print(print_line)


"""
Below is an example usage of the API.
"""
if __name__ == '__main__':
    import os

    from petrosim.models.ecafc import loader as ldr

    ldr.init(os.path.join(os.path.dirname(ldr.__file__), 'example.in'))

    Teq_norm = 0.8270
    params_eq = EquilibrationParams(Teq_norm)
    params_eq.printParams(header=True)
