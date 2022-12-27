"""
These are the helper functions, mostly differential equations used in the EC-AFC
model by Spera & Bohrson J. Petrol., 42, 999–1018, 2001.
"""

import functools

import numpy as np

import settings
#import constants_old as const


K = 273.16          # Conversion factor from ºC to K
R = 8.31432         # Gas constant


@functools.lru_cache
def unnormalize_temp(Tnorm):
    """
    Unnormalize the temperature

    :param Tnorm: The normalized temperature
    :type Tnorm: float
    """

    return settings.const.Tlm * Tnorm - settings.const.K


@functools.lru_cache
def tempK(Tnorm):
    return settings.const.Tlm * Tnorm


def Tm_slope(dT=0):
    return dT * 1e3


@functools.lru_cache
def phi(T, component):
    if component == settings.ASSIMILANT:
        Tlx = settings.const.Tla
    elif component == settings.MAGMA:
        Tlx = settings.const.Tlm
    return (T - settings.const.Ts) / (Tlx - settings.const.Ts)


@functools.lru_cache
def phi_prime(component):
    if component == settings.ASSIMILANT:
        Tlx = settings.const.Tla
    elif component == settings.MAGMA:
        Tlx = settings.const.Tlm
    return settings.const.Tlm / (Tlx - settings.const.Ts)


@functools.lru_cache
def melt_productivity(T, component, linear=False):
    alpha_x, beta_x = 0, 0
    if component == settings.ASSIMILANT:
        alpha_x = settings.const.alpha_a
        beta_x = settings.const.beta_a
    elif component == settings.MAGMA:
        alpha_x = settings.const.alpha_m
        beta_x = settings.const.beta_m

    if linear:
        fx = phi(T, component)
        if fx <= 0:
            fx = 0
        elif fx >= 1:
            fx = 1
    else:
        fx = 1 / (1 + (alpha_x * np.exp(beta_x * phi(T, component))))
    return fx


@functools.lru_cache
def melt_productivity_T_deriv(T, component, linear=False):
    alpha_x, beta_x = 0, 0
    if component == settings.ASSIMILANT:
        alpha_x = settings.const.alpha_a
        beta_x = settings.const.beta_a
        Tlx = settings.const.Tla

    elif component == settings.MAGMA:
        alpha_x = settings.const.alpha_m
        beta_x = settings.const.beta_m
        Tlx = settings.const.Tlm

    numer = -1 * alpha_x * beta_x * np.exp(beta_x * phi(T, component))
    denom = (alpha_x * np.exp(beta_x * phi(T, component)) + 1) ** 2
    fxprime = phi_prime(component) * numer / denom
    return fxprime


@functools.lru_cache
def mass_country_rock(T):
    temp = T / settings.const.Tlm
    spec_nrg_m = settings.const.cpm * (settings.const.Tm0 - settings.const.Tlm * temp)
    enth_cry_m = settings.const.dhm * (1 - melt_productivity(T, settings.MAGMA))
    spec_nrg_a = settings.const.cpa * (settings.const.Tlm * temp - settings.const.Ta0)
    enth_fus_a = settings.const.dha * melt_productivity(T, settings.ASSIMILANT)
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

    fm = melt_productivity(temp_m, settings.MAGMA)
    fa = melt_productivity(temp_a, settings.ASSIMILANT)
    fmprime = melt_productivity_T_deriv(temp_m, settings.MAGMA)
    faprime = melt_productivity_T_deriv(temp_a, settings.ASSIMILANT)

    term1 = -1 / Ma0
    term2_num = (settings.const.Tlm * settings.const.cpm) + (settings.const.dhm * fmprime) + (Ma0 * settings.const.cpa * settings.const.Tlm * fa)
    term2_den = (settings.const.Tlm * settings.const.cpa * (1 - fa)) + (settings.const.dha + settings.const.cpa * settings.const.Tlm * (Tm - Ta)) * faprime
    return term1 * term2_num / term2_den


@functools.lru_cache
def conserv_mass_dMm_dTm(Tm, Ta, Ma0):
    temp_m = tempK(Tm)
    temp_a = tempK(Ta)

    fmprime = melt_productivity_T_deriv(temp_m, settings.MAGMA)
    faprime = melt_productivity_T_deriv(temp_a, settings.ASSIMILANT)
    dTa_dTm = conserv_energy_dTa_dTm(Tm, Ta, Ma0)
    return Ma0 * faprime * dTa_dTm + fmprime


@functools.lru_cache
def conc_trace_elem_dCm_dTm(Tm, Ta, Ma0, Mm, Cm, dTa_dTm):
    temp_m = tempK(Tm)
    temp_a = tempK(Ta)

    fmprime = melt_productivity_T_deriv(temp_m, settings.MAGMA)
    faprime = melt_productivity_T_deriv(temp_a, settings.ASSIMILANT)
    Cabar = conc_trace_elem_in_anatectic_melt(Ta) / settings.const.Ca0
    Cmbar = Cm / settings.const.Cm0
    s = settings.const.Ca0 / settings.const.Cm0
    Dm = distribution_coefficient(Tm, settings.MAGMA)

    term1 = 1 / Mm
    term2 = Ma0 * (s * Cabar - Cmbar) * faprime * dTa_dTm
    term3 = Cmbar * (Dm - 1) * fmprime

    return term1 * (term2 + term3) * settings.const.Cm0


@functools.lru_cache
def conc_trace_elem_in_anatectic_melt(Ta):
    temp_a = tempK(Ta)

    fa = melt_productivity(temp_a, settings.ASSIMILANT)
    Da = distribution_coefficient(Ta, settings.ASSIMILANT)
    return (settings.const.Ca0 / Da) * ((1 - fa) ** ((1 - Da) / Da))


@functools.lru_cache
def isotop_ratio_dem_dTm(Tm, Ta, Ma0, Mm, Cm, dTa_dTm, em):
    temp_m = Tm * settings.const.Tlm
    temp_a = Ta * settings.const.Tlm

    fmprime = melt_productivity_T_deriv(temp_m, settings.MAGMA)
    faprime = melt_productivity_T_deriv(temp_a, settings.ASSIMILANT)
    Cabar = conc_trace_elem_in_anatectic_melt(Ta) / settings.const.Ca0
    Cmbar = Cm / settings.const.Cm0
    s = settings.const.Ca0 / settings.const.Cm0

    term1 = 1 / Mm
    term2 = (s * (Cabar / Cmbar) * (settings.const.ea0 - em) * Ma0 * faprime * dTa_dTm)
    return term1 * term2


@functools.lru_cache
def oxygen_isotope_comp_ddm_dTm(Tm, Ta, Ma0, Mm, dTa_dTm, dm):
    temp_m = Tm * settings.const.Tlm
    temp_a = Ta * settings.const.Tlm

    faprime = melt_productivity_T_deriv(temp_a, settings.ASSIMILANT)

    term1 = 1 / Mm
    term2 = (settings.const.da0 - dm) * settings.const.koxy * Ma0 * faprime * dTa_dTm
    return term1 * term2


def distribution_coefficient(Tx, component):
    if component == settings.ASSIMILANT:
        Dx0 = settings.const.Da0
        dHx = settings.const.dHa
    elif component == settings.MAGMA:
        Dx0 = settings.const.Dm0
        dHx = settings.const.dHm
    return Dx0 * np.exp(-dHx / (settings.const.R * settings.const.Tlm * Tx))


class EquilibrationParams:

    def __init__(self, Tnorm):
        self.Tnorm = Tnorm
        self.Teq = settings.const.Tlm * Tnorm - settings.const.K
        self.fm = melt_productivity(self.Teq + settings.const.K, settings.MAGMA)
        self.fa = melt_productivity(self.Teq + settings.const.K, settings.ASSIMILANT)
        self.Ma0 = mass_country_rock(self.Teq + settings.const.K)
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
