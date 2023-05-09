import os
import unittest

import pytest

from petrosim.models.ecrafc import equilibration as equil
from petrosim.models.ecrafc import loader as ldr

EXAMPLE_INFILE = os.path.join(os.path.dirname(__file__), 'test.in')


class TestTopLevelFuncs(unittest.TestCase):
    def setUp(self):
        self.Tnorm = 0.827
        self.T = 950
        self.component = ldr.MAGMA
        self.Ma0 = 0.910
        self.Mastar = 0.653
        self.Mct = 0.969
        self.Mr = 0.000
        self.Ta = 1.000
        self.Tm = 0.548
        self.Mm = 1.000
        self.Cm = 700.00
        self.Ca0 = 350.00
        self.Cm0 = 700.00
        self.Cr0 = 1000.00
        self.ea0 = 0.722
        self.er0 = 0.703
        self.em = 0.7035
        self.dm = 1.0
        self.fm = 0.1

        ldr.init(EXAMPLE_INFILE)

    def test_unnormalize_temp(self):
        expect = 1027.84
        assert equil.unnormalize_temp(self.Tnorm) == pytest.approx(
            expect, 0.01)

    def test_tempK(self):
        expect = 1301.0
        assert equil.tempK(self.Tnorm) == pytest.approx(expect, 0.01)

    def test_Tm_slope(self):
        dT = 0.001
        expect = 1
        assert equil.Tm_slope(dT) == expect

    def test_phi(self):
        expect = -0.558
        assert equil.phi(self.T, self.component) == pytest.approx(expect, 0.01)

    def test_phi_prime(self):
        expect = 3.93
        assert equil.phi_prime(self.component) == pytest.approx(expect, 0.01)

    def test_melt_productivity(self):
        expect = 4.8e-6
        assert equil.melt_productivity(self.T,
                                       self.component) == pytest.approx(
                                           expect, 0.01)

    def test_melt_productivity_T_deriv(self):
        expect = 2.08e-4
        assert equil.melt_productivity_T_deriv(
            self.T, self.component) == pytest.approx(expect, 0.01)

    def test_mass_country_rock_Ma0(self):
        expect = -2.98
        assert equil.mass_country_rock_Ma0(self.Tnorm) == pytest.approx(
            expect, 0.01)

    def test_mass_anatectic_melt_Mastar(self):
        fa = 0.717
        expect = 0.652
        assert equil.mass_anatectic_melt_Mastar(self.Ma0, fa) == pytest.approx(
            expect, 0.01)

    def test_mass_cumulates(self):
        fm = 0.031
        expect = 0.969
        assert equil.mass_cumulates(fm) == pytest.approx(expect, 0.01)

    def test_mass_magma(self):
        expect = 1.04
        assert equil.mass_magma(self.Mastar, self.Mct,
                                self.fm) == pytest.approx(expect, 0.01)

    def test_ratio_country_rock_to_cumulates(self):
        expect = 0.939
        assert equil.ratio_country_rock_to_cumulates(
            self.Ma0, self.Mct) == pytest.approx(expect, 0.01)

    def test_ratio_anatectic_to_cumulates(self):
        expect = 0.674
        assert equil.ratio_anatectic_to_cumulates(self.Mastar,
                                                  self.Mct) == pytest.approx(
                                                      expect, 0.01)

    def test_conserv_energy_dTa_dTm(self):
        dMr_dTm = equil.variation_recharge_mass_dMr_dTm(self.Tm)
        expect = -1.255
        assert equil.conserv_energy_dTa_dTm(self.Ta, self.Tm, self.Ma0,
                                            self.Mr, dMr_dTm) == pytest.approx(
                                                expect, 0.01)

    def test_conserv_mass_dMm_dTm(self):
        dMr_dTm = equil.variation_recharge_mass_dMr_dTm(self.Tm)
        expect = 0.320
        assert equil.conserv_mass_dMm_dTm(self.Ta, self.Tm, self.Ma0, self.Mr,
                                          dMr_dTm) == pytest.approx(
                                              expect, 0.01)

    def test_conc_trace_elem_dCm_dTm(self):
        dMr_dTm = equil.variation_recharge_mass_dMr_dTm(self.Tm)
        dTa_dTm = equil.conserv_mass_dMm_dTm(self.Ta, self.Tm, self.Ma0,
                                             self.Mr, dMr_dTm)
        expect = 112.11
        assert equil.conc_trace_elem_dCm_dTm(
            self.Ta, self.Tm, self.Ma0, self.Mm, self.Cm, self.Ca0, self.Cm0,
            self.Cr0, dTa_dTm, self.Mr,
            dMr_dTm) == pytest.approx(expect, 0.01)

    def test_conc_trace_elem_in_anatectic_melt(self):
        expect = 1.14e-92
        assert equil.conc_trace_elem_in_anatectic_melt(
            self.Ta, self.Ca0) == pytest.approx(expect, 0.01)

    def test_isotop_ratio_dem_dTm(self):
        dMr_dTm = equil.variation_recharge_mass_dMr_dTm(self.Tm)
        dTa_dTm = equil.conserv_mass_dMm_dTm(self.Ta, self.Tm, self.Ma0,
                                             self.Mr, dMr_dTm)
        expect = 1.06e-8
        assert equil.isotop_ratio_dem_dTm(self.Ta, self.Tm, self.Ma0, self.Mm,
                                          self.Cm, self.Ca0, self.Cm0,
                                          self.Cr0, dTa_dTm, dMr_dTm, self.ea0,
                                          self.er0, self.em) == pytest.approx(
                                              expect, 0.01)

    def test_oxygen_isotope_comp_ddm_dTm(self):
        dMr_dTm = equil.variation_recharge_mass_dMr_dTm(self.Tm)
        dTa_dTm = equil.conserv_mass_dMm_dTm(self.Ta, self.Tm, self.Ma0,
                                             self.Mr, dMr_dTm)
        expect = 0.002
        assert equil.oxygen_isotope_comp_ddm_dTm(self.Tm, self.Ta, self.Ma0,
                                                 self.Mm, dTa_dTm, dMr_dTm,
                                                 self.dm) == pytest.approx(
                                                     expect, 0.01)

    def test_distribution_coefficient(self):
        expect = 1.5
        assert equil.distribution_coefficient(self.Ta,
                                              self.component) == pytest.approx(
                                                  expect, 0.01)


class TestEquilibrationParams(unittest.TestCase):
    def test_EqulibrationParams(self):
        Teq_norm = 0.8270
        params_eq = equil.EquilibrationParams(Teq_norm)
        expect = 1027.84
        assert params_eq.Teq == pytest.approx(expect, 0.01)
