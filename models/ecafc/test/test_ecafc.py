import dataclasses
import unittest

import pytest

from petrosim.models.ecafc import ecafc
from petrosim.models.ecafc import equilibration as equil
from petrosim.models.ecafc import loader as ldr

EXAMPLE_INFILE = 'test.in'


class TestTopLevelFuncs(unittest.TestCase):
    def test_get_T_iter(self):
        Tm = 1.000
        dT = 0.001
        slope_iter = 2
        expect = 0.995
        assert ecafc.get_T_iter(Tm, dT,
                                slope_iter) == pytest.approx(expect, 0.01)


class TestParameter(unittest.TestCase):
    def setUp(self):
        ldr.init(EXAMPLE_INFILE)
        self.params_init = ldr.Parameters(**dataclasses.asdict(ldr.params))
        self.params_eq = equil.EquilibrationParams(self.params_init.Teq_norm)
        self.params_sim = {
            'Tm': ecafc.Parameter_Tm('Tm,norm', 1.000),
            'Ta': ecafc.Parameter_Ta('Ta,norm', 0.548),
            'Mm': ecafc.Parameter_Mm('Mm', 1.000)
        }
        self.params_traces = [{
            'Cm': ecafc.Parameter_Cm('Cm', 700.0, decimals=3)
        }]

        self.Tm = 1.000
        self.dT = 0.001

    def test_Parameter(self):
        name = 'Ta_norm'
        name_alt = 'Ta'
        value_old = 0.548
        decimals = 2
        parameter = ecafc.Parameter(name, value_old, decimals, name_alt)
        assert parameter.str.strip() == '0.55'
        assert parameter.slopes == [0.0]

    def test_thisStepValue(self):
        for h in range(4):
            T = ecafc.get_T_iter(self.Tm, self.dT, h)
            param = self.params_sim['Ta']
            param.calcSlope(h, T, self.params_init, self.params_eq,
                            self.params_sim)
            param.addSlope(self.dT)

        slope_iter = 2
        expect = 0.549
        assert self.params_sim['Ta'].thisStepValue(
            slope_iter) == pytest.approx(expect, 0.01)

    def test_CalcSlope(self):
        param = self.params_sim['Ta']

        slope_iters = 2
        for h in range(slope_iters):
            T = ecafc.get_T_iter(self.Tm, self.dT, h)
            param.calcSlope(h, T, self.params_init, self.params_eq,
                            self.params_sim)
            param.addSlope(self.dT)

        expect = -1.28
        assert param.dy_dx == pytest.approx(expect, 0.01)

    def test_CalcSlopeTrace(self):
        param_Ta = self.params_sim['Ta']
        param_Mm = self.params_sim['Mm']
        param_Cm = self.params_traces[0]['Cm']

        slope_iters = 2
        for h in range(slope_iters):
            T = ecafc.get_T_iter(self.Tm, self.dT, h)
            param_Ta.calcSlope(h, T, self.params_init, self.params_eq,
                               self.params_sim)
            param_Ta.addSlope(self.dT)
            param_Mm.calcSlope(h, T, self.params_init, self.params_eq,
                               self.params_sim)
            param_Mm.addSlope(self.dT)
            param_Cm.calcSlopeTrace(h, T, self.params_init, self.params_eq,
                                    self.params_sim, self.params_traces[0], 0)
            param_Cm.addSlope(self.dT)

        expect = 195.60
        assert param_Cm.dy_dx == pytest.approx(expect, 0.01)

    def test_addSlope(self):
        param_Ta = self.params_sim['Ta']
        param_Ta.dy_dx = 5_000
        differential = 0.001
        assert param_Ta.slopes == [0.0]
        param_Ta.addSlope(differential)
        assert param_Ta.slopes == [0.0, 5.0]

    def test_clearSlopes(self):
        param_Ta = self.params_sim['Ta']
        param_Ta.slopes = [0.0, 1.0, 2.0]
        assert len(param_Ta.slopes) == 3
        param_Ta.clearSlopes()
        assert len(param_Ta.slopes) == 1
        assert param_Ta.slopes == [0.0]

    def test_calcNewValue(self):
        param_Ta = self.params_sim['Ta']
        param_Ta.value_old = 1.0
        param_Ta.slopes = [0.0, 0.1, 0.1, 0.1, 0.1]
        param_Ta.calcNewValue()
        assert param_Ta.value_new == 0.9

    def test_updateValue(self):
        param_Ta = self.params_sim['Ta']
        param_Ta.value_old = 1.0
        param_Ta.value_new = 0.9
        assert param_Ta.value_old == 1.0
        param_Ta.updateValue()
        assert param_Ta.value_old == 0.9


class TestECAFC(unittest.TestCase):
    def setUp(self):
        ldr.init(EXAMPLE_INFILE)
        self.params_init = ldr.Parameters(**dataclasses.asdict(ldr.params))
        self.ecafc = ecafc.ECAFC(self.params_init)

    def test_simulate(self):
        self.ecafc.simulate()

        dT = self.params_init.dT
        Tnorm0 = self.params_init.Teq_norm - dT
        Tnorm1 = self.params_init.Tm0_norm
        max_iter = int((Tnorm1 - Tnorm0) / dT + 1)
        assert len(self.ecafc.results.results) == max_iter
