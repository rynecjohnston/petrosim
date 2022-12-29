import unittest

import pytest

from petrosim.models.ecafc import equilibration as equil
from petrosim.models.ecafc import loader as ldr

EXAMPLE_INFILE = 'test.in'


class TestTopLevelFuncs(unittest.TestCase):
    def setUp(self):
        self.Tnorm = 0.827
        self.T = 950
        self.component = ldr.MAGMA
        self.Ma0 = 0.910
        self.Mastar = 0.653
        self.Mct = 0.969
        self.Ta = 1.000
        self.Tm = 0.548
        self.Mm = 1.000
        self.Cm = 700.00
        self.Ca0 = 350.00
        self.Cm0 = 700.00
        self.ea0 = 0.722
        self.em = 0.7035
        self.dm = 1.0

        ldr.init(EXAMPLE_INFILE)

    def test_unnormalize_temp(self):
        expect = 1044.38
        assert equil.unnormalize_temp(self.Tnorm) == pytest.approx(
            expect, 0.01)

    def test_tempK(self):
        expect = 1317.53
        assert equil.tempK(self.Tnorm) == pytest.approx(expect, 0.01)

    def test_Tm_slope(self):
        dT = 0.001
        expect = 1
        assert equil.Tm_slope(dT) == expect

    def test_phi(self):
        expect = -0.738
        assert equil.phi(self.T, self.component) == pytest.approx(expect, 0.01)

    def test_phi_prime(self):
        expect = 4.31
        assert equil.phi_prime(self.component) == pytest.approx(expect, 0.01)

    def test_melt_productivity(self):
        expect = 9.5e-7
        assert equil.melt_productivity(self.T,
                                       self.component) == pytest.approx(
                                           expect, 0.01)

    def test_melt_productivity_T_deriv(self):
        expect = 4.3e-5
        assert equil.melt_productivity_T_deriv(
            self.T, self.component) == pytest.approx(expect, 0.01)

    def test_mass_country_rock(self):
        expect = -2.28
        assert equil.mass_country_rock(self.Tnorm) == pytest.approx(
            expect, 0.01)

    def test_mass_anatectic_melt(self):
        fa = 0.717
        expect = 0.652
        assert equil.mass_anatectic_melt(self.Ma0,
                                         fa) == pytest.approx(expect, 0.01)

    def test_mass_cumulates(self):
        fm = 0.031
        expect = 0.969
        assert equil.mass_cumulates(fm) == pytest.approx(expect, 0.01)

    def test_mass_magma(self):
        expect = 0.684
        assert equil.mass_magma(self.Mastar,
                                self.Mct) == pytest.approx(expect, 0.01)

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
        expect = -1.282
        assert equil.conserv_energy_dTa_dTm(self.Ta, self.Tm,
                                            self.Ma0) == pytest.approx(
                                                expect, 0.01)

    def test_conserv_mass_dMm_dTm(self):
        expect = 0.547
        assert equil.conserv_mass_dMm_dTm(self.Ta, self.Tm,
                                          self.Ma0) == pytest.approx(
                                              expect, 0.01)

    def test_conc_trace_elem_dCm_dTm(self):
        dTa_dTm = equil.conserv_mass_dMm_dTm(self.Ta, self.Tm, self.Ma0)
        expect = 191.31
        assert equil.conc_trace_elem_dCm_dTm(
            self.Ta, self.Tm, self.Ma0, self.Mm, self.Cm, self.Ca0, self.Cm0,
            dTa_dTm) == pytest.approx(expect, 0.01)

    def test_conc_trace_elem_in_anatectic_melt(self):
        expect = 2.68e5
        assert equil.conc_trace_elem_in_anatectic_melt(
            self.Ta, self.Ca0) == pytest.approx(expect, 0.01)

    def test_isotop_ratio_dem_dTm(self):
        dTa_dTm = equil.conserv_mass_dMm_dTm(self.Ta, self.Tm, self.Ma0)
        expect = 6.34e-15
        assert equil.isotop_ratio_dem_dTm(self.Ta, self.Tm, self.Ma0, self.Mm,
                                          self.Cm, self.Ca0, self.Cm0, dTa_dTm,
                                          self.ea0, self.em) == pytest.approx(
                                              expect, 0.01)

    def test_oxygen_isotope_comp_ddm_dTm(self):
        dTa_dTm = equil.conserv_mass_dMm_dTm(self.Ta, self.Tm, self.Ma0)
        expect = 0.0
        assert equil.oxygen_isotope_comp_ddm_dTm(self.Ta, self.Tm, self.Ma0,
                                                 self.Mm, dTa_dTm,
                                                 self.dm) == expect

    def test_distribution_coefficient(self):
        expect = 1.5
        assert equil.distribution_coefficient(self.Ta,
                                              self.component) == pytest.approx(
                                                  expect, 0.01)


class TestEquilibrationParams(unittest.TestCase):
    def test_EqulibrationParams(self):
        Teq_norm = 0.8270
        params_eq = equil.EquilibrationParams(Teq_norm)
        expect = 1044.383
        assert params_eq.Teq == pytest.approx(expect, 0.01)