"""
These are the helper functions, mostly differential equations used in the EC-AFC
model by Spera & Bohrson J. Petrol., 42, 999–1018, 2001.
"""

import functools

import numpy as np

import loader as ldr


K = 273.16          # Conversion factor from ºC to K
R = 8.31432         # Gas constant


@functools.lru_cache
def normalize_temp(temp):
    """
    Normalize the temperature

    :param Tnorm: The temperature in ºC
    :type Tnorm: float
    """

    return (temp + ldr.params.K) / ldr.params.Tlm


@functools.lru_cache
def unnormalize_temp(Tnorm):
    """
    Unnormalize the temperature

    :param Tnorm: The normalized temperature
    :type Tnorm: float
    """

    return ldr.params.Tlm * Tnorm - ldr.params.K


@functools.lru_cache
def tempK(Tnorm):
    return ldr.params.Tlm * Tnorm


@functools.lru_cache
def Tm_slope(dT=0):
    return dT * 1e3


@functools.lru_cache
def phi(T, component):
    if component == ldr.ASSIMILANT:
        Tlx = ldr.params.Tla
    elif component == ldr.MAGMA:
        Tlx = ldr.params.Tlm
    return (T - ldr.params.Ts) / (Tlx - ldr.params.Ts)


@functools.lru_cache
def phi_prime(component):
    if component == ldr.ASSIMILANT:
        Tlx = ldr.params.Tla
    elif component == ldr.MAGMA:
        Tlx = ldr.params.Tlm
    return ldr.params.Tlm / (Tlx - ldr.params.Ts)


@functools.lru_cache
def melt_productivity(T, component, linear=False):
    loga_x, logb_x = 0, 0
    if component == ldr.ASSIMILANT:
        loga_x = ldr.params.loga_a
        logb_x = ldr.params.logb_a
    elif component == ldr.MAGMA:
        loga_x = ldr.params.loga_m
        logb_x = ldr.params.logb_m

    if linear:
        fx = phi(T, component)
        if fx <= 0:
            fx = 0
        elif fx >= 1:
            fx = 1
    else:
        fx = 1 / (1 + (loga_x * np.exp(logb_x * phi(T, component))))
    return fx


@functools.lru_cache
def melt_productivity_T_deriv(T, component, linear=False):
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
    denom = (loga_x * np.exp(logb_x * phi(T, component)) + 1) ** 2
    fxprime = phi_prime(component) * numer / denom
    return fxprime


@functools.lru_cache
def mass_country_rock(T):
    temp = T / ldr.params.Tlm
    spec_nrg_m = ldr.params.cpm * (ldr.params.Tm0 - ldr.params.Tlm * temp)
    enth_cry_m = ldr.params.dhm * (1 - melt_productivity(T, ldr.MAGMA))
    spec_nrg_a = ldr.params.cpa * (ldr.params.Tlm * temp - ldr.params.Ta0)
    enth_fus_a = ldr.params.dha * melt_productivity(T, ldr.ASSIMILANT)
    return (spec_nrg_m + enth_cry_m) / (spec_nrg_a + enth_fus_a)


def mass_anatectic_melt(Ma0, fa):
    return Ma0 * fa


def mass_cumulates(fm):
    return 1 - fm


def mass_magma_initial(Ma0):
    return -1 / Ma0


def mass_magma(Mastar, Mct):
    return 1 + Mastar - Mct


def ratio_country_rock_to_cumulates(Ma0, Mct):
    return Ma0 / Mct


def ratio_anatectic_to_cumulates(Mastar, Mct):
    return Mastar / Mct


@functools.lru_cache
def conserv_energy_dTa_dTm(Tm, Ta, Ma0):
    temp_m = tempK(Tm)
    temp_a = tempK(Ta)

    fm = melt_productivity(temp_m, ldr.MAGMA)
    fa = melt_productivity(temp_a, ldr.ASSIMILANT)
    fmprime = melt_productivity_T_deriv(temp_m, ldr.MAGMA)
    faprime = melt_productivity_T_deriv(temp_a, ldr.ASSIMILANT)

    term1 = -1 / Ma0
    term2_num = (ldr.params.Tlm * ldr.params.cpm) + (ldr.params.dhm * fmprime) + (Ma0 * ldr.params.cpa * ldr.params.Tlm * fa)
    term2_den = (ldr.params.Tlm * ldr.params.cpa * (1 - fa)) + (ldr.params.dha + ldr.params.cpa * ldr.params.Tlm * (Tm - Ta)) * faprime
    return term1 * term2_num / term2_den


@functools.lru_cache
def conserv_mass_dMm_dTm(Tm, Ta, Ma0):
    temp_m = tempK(Tm)
    temp_a = tempK(Ta)

    fmprime = melt_productivity_T_deriv(temp_m, ldr.MAGMA)
    faprime = melt_productivity_T_deriv(temp_a, ldr.ASSIMILANT)
    dTa_dTm = conserv_energy_dTa_dTm(Tm, Ta, Ma0)
    return Ma0 * faprime * dTa_dTm + fmprime


@functools.lru_cache
def conc_trace_elem_dCm_dTm(Tm, Ta, Ma0, Mm, Cm, Ca0, Cm0, dTa_dTm):
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
    temp_a = tempK(Ta)

    fa = melt_productivity(temp_a, ldr.ASSIMILANT)
    Da = distribution_coefficient(Ta, ldr.ASSIMILANT)
    return (Ca0 / Da) * ((1 - fa) ** ((1 - Da) / Da))


@functools.lru_cache
def isotop_ratio_dem_dTm(Tm, Ta, Ma0, Mm, Cm, Ca0, Cm0, dTa_dTm, ea0, em):
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
    temp_m = Tm * ldr.params.Tlm
    temp_a = Ta * ldr.params.Tlm

    faprime = melt_productivity_T_deriv(temp_a, ldr.ASSIMILANT)

    term1 = 1 / Mm
    term2 = (ldr.params.da0 - dm) * ldr.params.koxy * Ma0 * faprime * dTa_dTm
    return term1 * term2


def distribution_coefficient(Tx, component):
    if component == ldr.ASSIMILANT:
        Dx0 = ldr.params.Da0
        dHx = ldr.params.dHa
    elif component == ldr.MAGMA:
        Dx0 = ldr.params.Dm0
        dHx = ldr.params.dHm
    return Dx0 * np.exp(-dHx / (ldr.params.R * ldr.params.Tlm * Tx))


class EquilibrationParams:

    def __init__(self, Tnorm):
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
            print_line = '  '.join([f'{data:>{length}s}' for data, length in data_len.items()])
            print(print_line)
        print_line = '  '.join([f'{getattr(self, data):>{length}.3f}' for data, length in data_len.items()])
        print(print_line)


if __name__ == '__main__':
    Teq_norm = 0.9040
    params_eq = EquilibrationParams(Teq_norm)
    params_eq.printParams(header=True)
