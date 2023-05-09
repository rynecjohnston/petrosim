"""
These are the helper functions, mostly differential equations used in the Energy
-Constrained Recharge, Assimilation, and Fractional Crystallization (EC-RAFC)
model based on work by Spera & Bohrson G3, 3, 12, 2002.
"""

import functools

import numpy as np

from petrosim.models.ecrafc import loader as ldr


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
    elif component == ldr.RECHARGE:
        Tlx = ldr.params.Tlr
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
    elif component == ldr.RECHARGE:
        Tlx = ldr.params.Tlr
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
    elif component == ldr.RECHARGE:
        loga_x = ldr.params.loga_r
        logb_x = ldr.params.logb_r

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
    elif component == ldr.RECHARGE:
        loga_x = ldr.params.loga_r
        logb_x = ldr.params.logb_r
        Tlx = ldr.params.Tlr

    numer = -1 * loga_x * logb_x * np.exp(logb_x * phi(T, component))
    denom = (loga_x * np.exp(logb_x * phi(T, component)) + 1)**2
    fxprime = phi_prime(component) * numer / denom
    return fxprime


@functools.lru_cache
def mass_country_rock_Ma0(T):
    """
    :param T: The temperature
    :type T: float

    :return: Mass of anatectic melt at given a temperature
    :rtype: float
    """

    temp = T / ldr.params.Tlm
    fm = melt_productivity(T, ldr.MAGMA)
    fa = melt_productivity(T, ldr.ASSIMILANT)
    fr = melt_productivity(T, ldr.RECHARGE)

    spec_nrg_m = ldr.params.cpm * (ldr.params.Tm0 - ldr.params.Tlm * temp)
    enth_cry_m = ldr.params.dhm * (1 - fm)
    spec_nrg_a = ldr.params.cpa * (ldr.params.Tlm * temp - ldr.params.Ta0)
    enth_fus_a = ldr.params.dha * fa
    spec_nrg_r = ldr.params.Mr0 * ldr.params.cpr * (ldr.params.Tr0 -
                                                    ldr.params.Tlm * temp)
    enth_cry_r = ldr.params.Mr0 * ldr.params.dhr * (1 - fr)
    return (spec_nrg_m + enth_cry_m + spec_nrg_r + enth_cry_r) / (spec_nrg_a +
                                                                  enth_fus_a)


def total_mass_recharge_RAFC_event(T):
    return melt_productivity(T, ldr.RECHARGE)


def cumulative_recharge_mass_Mr(T):
    if ldr.params.rechg_linear:
        numer = ldr.params.Mr0 * (T - ldr.params.Tm0)
        denom = (ldr.params.Teq_norm * ldr.params.Tlm) - ldr.params.Tm0
        Mr = numer / denom
    else:
        Mr = 0
        for i in range(ldr.params.Nr):
            numer = ldr.params.rechg_dMr[i]
            denom = (
                1 +
                np.exp(ldr.params.rechg_m[i] *
                       (T - ldr.params.rechg_Tr[i])))**ldr.params.rechg_d[i]
            Mr += numer / denom
    return Mr


def mass_anatectic_melt_Mastar(Ma0, fa):
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
    :param fm: Melt productivity of magma
    :type fm: float

    :return: Mass of cumulates
    :rtype: float
    """

    return 1 - fm


def mass_enclaves(fr, Mr):
    """
    :param fa: Melt productivity of recharge
    :type fa: float

    :return: Mass of enclaves
    :rtype: float
    """

    return Mr * (1 - fr)


def mass_magma(Mastar, fr, fm):
    """
    :param Mastar: Mass of anatectic melt
    :type Mastar: float

    :param fm: Melt productivity of magma
    :type fm: float

    :param fr: Melt productivity of recharge
    :type fr: float

    :return: Mass of magma
    :rtype: float
    """

    return Mastar + fr * ldr.params.Mr0 + fm


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
def variation_recharge_mass_dMr_dTm(T):
    if ldr.params.rechg_linear:
        numer = ldr.params.Mr0 * ldr.params.Tlm
        denom = (ldr.params.Teq_norm * ldr.params.Tlm) - ldr.params.Tm0
        dMr_dTm = numer / denom
    else:
        dMr_dTm = 0
        for i in range(ldr.params.Nr):
            term1 = ldr.params.Tlm * ldr.params.rechg_m[
                i] * ldr.params.rechg_d[i] * ldr.params.rechg_dMr[i] * np.exp(
                    ldr.params.rechg_m[i] *
                    (ldr.params.Tlm * T - ldr.params.rechg_Tr[i]))
            term2 = (1 +
                     np.exp(ldr.params.rechg_m[i] *
                            (ldr.params.Tlm * T - ldr.params.rechg_Tr[i])))**(
                                -ldr.params.rechg_d[i] - 1)
            dMr_dTm -= term1 * term2
    return dMr_dTm


@functools.lru_cache
def conserv_energy_dTa_dTm(Tm, Ta, Ma0, Mr, dMr_dTm):
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
    fr = melt_productivity(temp_m, ldr.RECHARGE)
    fmprime = melt_productivity_T_deriv(temp_m, ldr.MAGMA)
    faprime = melt_productivity_T_deriv(temp_a, ldr.ASSIMILANT)
    frprime = melt_productivity_T_deriv(temp_m, ldr.RECHARGE)
    Tr0bar = ldr.params.Tr0 / ldr.params.Tlm

    term1 = -1 / Ma0
    term2_num = (ldr.params.Tlm * ldr.params.cpm) + (
        ldr.params.dhm * fmprime) + (Ma0 * ldr.params.cpa * ldr.params.Tlm *
                                     fa)
    """
    In personal correspondence, F. Spera confirmed that the following term from
        eqn 10 from the Spera 2002 paper is incorrect:
    term3_num = ldr.params.Tlm * ldr.params.cpr + ldr.params.dhr * frprime * Mr
    The whole numerator sum must be multiplied by `Mr`
    """
    term3_num = ((ldr.params.Tlm * ldr.params.cpr) +
                 (ldr.params.dhr * frprime)) * Mr
    term4_num = (ldr.params.Tlm * ldr.params.cpr * (Tr0bar - Tm) +
                 ((1 - fr) * ldr.params.dhr)) * -dMr_dTm
    common_den = (
        ldr.params.Tlm * ldr.params.cpa *
        (1 - fa)) + (ldr.params.dha + ldr.params.cpa * ldr.params.Tlm *
                     (Tm - Ta)) * faprime
    return term1 * (term2_num + term3_num + term4_num) / common_den


@functools.lru_cache
def conserv_mass_dMm_dTm(Tm, Ta, Ma0, Mr, dMr_dTm):
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

    fr = melt_productivity(temp_m, ldr.RECHARGE)
    fmprime = melt_productivity_T_deriv(temp_m, ldr.MAGMA)
    faprime = melt_productivity_T_deriv(temp_a, ldr.ASSIMILANT)
    frprime = melt_productivity_T_deriv(temp_m, ldr.RECHARGE)

    dTa_dTm = conserv_energy_dTa_dTm(Tm, Ta, Ma0, Mr, dMr_dTm)
    return Ma0 * faprime * dTa_dTm + fmprime + frprime * Mr + fr * dMr_dTm


@functools.lru_cache
def conc_trace_elem_dCm_dTm(Tm, Ta, Ma0, Mm, Cm, Ca0, Cm0, Cr0, dTa_dTm, Mr,
                            dMr_dTm):
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

    fr = melt_productivity(temp_m, ldr.RECHARGE)
    fmprime = melt_productivity_T_deriv(temp_m, ldr.MAGMA)
    faprime = melt_productivity_T_deriv(temp_a, ldr.ASSIMILANT)
    frprime = melt_productivity_T_deriv(temp_m, ldr.RECHARGE)
    Dm = distribution_coefficient(Tm, ldr.MAGMA)
    Dr = distribution_coefficient(Tm, ldr.RECHARGE)
    Cabar = conc_trace_elem_in_anatectic_melt(Ta, Ca0) / Ca0
    Cmbar = Cm / Cm0
    s = Ca0 / Cm0
    t = Cr0 / Cm0

    term1 = 1 / Mm
    term2 = Ma0 * (s * Cabar - Cmbar) * faprime * dTa_dTm
    term3 = Cmbar * (Dm - 1) * (fmprime + Mr * frprime)
    """
    Confused yet again:
    This is the term4 from eqn 13 from the Spera 2002 paper, but it produces incorrect results:
    term4 = t * (fr ** (Dr - 1) - Cmbar) * fr * dMr_dTm
    It seems that `t` should be included inside the parentheses, to parallel the
        `(s * Cabar - Cmbar)` in term2 but for recharge
    """
    term4 = (t * fr**(Dr - 1) - Cmbar) * fr * dMr_dTm
    return term1 * (term2 + term3 + term4) * Cm0


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
def isotop_ratio_dem_dTm(Tm, Ta, Ma0, Mm, Cm, Ca0, Cm0, Cr0, dTa_dTm, dMr_dTm,
                         ea0, er0, em):
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

    fr = melt_productivity(temp_m, ldr.RECHARGE)
    fmprime = melt_productivity_T_deriv(temp_m, ldr.MAGMA)
    faprime = melt_productivity_T_deriv(temp_a, ldr.ASSIMILANT)
    Dr = distribution_coefficient(Tm, ldr.RECHARGE)
    Cabar = conc_trace_elem_in_anatectic_melt(Ta, Ca0) / Ca0
    Cmbar = Cm / Cm0
    s = Ca0 / Cm0
    t = Cr0 / Cm0

    term1 = 1 / Mm
    term2 = (s * (Cabar / Cmbar) * (ea0 - em) * Ma0 * faprime * dTa_dTm)
    term3 = (t / Cmbar) * ((fr**(Dr - 1)) * (er0 - em) * fr * dMr_dTm)
    return term1 * (term2 + term3)


@functools.lru_cache
def oxygen_isotope_comp_ddm_dTm(Tm, Ta, Ma0, Mm, dTa_dTm, dMr_dTm, dm):
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

    fr = melt_productivity(temp_m, ldr.RECHARGE)
    faprime = melt_productivity_T_deriv(temp_a, ldr.ASSIMILANT)

    term1 = 1 / Mm
    term2 = (ldr.params.da0 - dm) * ldr.params.koxy * Ma0 * faprime * dTa_dTm
    term3 = (ldr.params.dr0 - dm) * ldr.params.Lamr_Lamm * fr * dMr_dTm
    return term1 * (term2 + term3)


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
    elif component == ldr.RECHARGE:
        Dx0 = ldr.params.Dr0
        dHx = ldr.params.dHr
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
        self.fr = melt_productivity(self.Teq + ldr.params.K, ldr.RECHARGE)
        self.Ma0 = mass_country_rock_Ma0(self.Teq + ldr.params.K)
        self.Mr0 = total_mass_recharge_RAFC_event(self.Teq + ldr.params.K)
        self.Mr = cumulative_recharge_mass_Mr(self.Teq + ldr.params.K)
        self.Mastar = mass_anatectic_melt_Mastar(self.Ma0, self.fa)
        self.Mct = mass_cumulates(self.fm)
        self.Men = mass_enclaves(self.fr, self.Mr)
        self.Ms = self.Men + self.Mct
        self.Mm = mass_magma(self.Mastar, self.fr, self.fm)
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
            'Mm': 5,
            'Ma0': 5,
            'Mastar': 6,
            'Mct': 5,
            'Men': 5,
            'Ms': 5,
            'fa': 5,
            'fm': 5,
            'fr': 5,
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

    from petrosim.models.ecrafc import loader as ldr

    ldr.init(os.path.join(os.path.dirname(ldr.__file__), 'example.in'))

    Teq_norm = 0.904
    params_eq = EquilibrationParams(Teq_norm)
    params_eq.printParams(header=True)
