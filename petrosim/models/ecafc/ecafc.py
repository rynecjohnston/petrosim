"""
This is the main simulation routine for the Energy-Constrained Assimilation
Fractional Crystallization (EC-AFC) model based on work by Spera & Bohrson J.
Petrol., 42, 999–1018, 2001.
"""

import copy

import numpy as np

from petrosim.models.ecafc import equilibration as equil
from petrosim.models.ecafc import results


def get_T_iter(Tm, dT, slope_iter):
    """
    Get the incremental temperature at this iteration.

    :param Tm: Temperature of the standing melt
    :type Tm: float

    :param dT: Temperature increment
    :type dT: float

    :param slope_iter: The slope iteration
    :type slope_iter: int

    :return: the incremental temperature at this iteration
    :rtype: float
    """

    if slope_iter == 0:
        T = Tm
    elif slope_iter == 3:
        T = Tm - dT
    else:
        T = Tm - dT / 2
    return T


class Parameter:
    """
    Base class for each parameter.
    """
    def __init__(self, name, value_old, decimals=3, name_alt=None):
        """
        :param value_old: Initial value for this iteration
        :type value_old: float

        :param decimals: Number of decimals to print in the output string
        :type decimals: int
        """

        self.name = name
        self.name_alt = name_alt
        self.value_old = value_old
        self.decimals = decimals
        self.str = f'{self.value_old:8.{self.decimals}f}'
        self.value_new = None
        self.dy_dx = None
        self.slopes = [0.0]

    def thisStepValue(self, slope_iter):
        """
        Calculate the present value for this slope iteration.

        :param slope_iter: The slope iteration
        :type slope_iter: int

        :return: The present value at this slope iteration
        :type: float
        """

        if slope_iter == 3:
            return self.value_old - self.slopes[slope_iter]
        else:
            return self.value_old - self.slopes[slope_iter] / 2

    def calcSlope(self, slope_iter, x, params_init, params_eq, params_sim):
        """
        Calculate the slope at point x using a particular function.

        :param slope_iter: The slope iteration
        :type slope_iter: int

        :param x: The x value for where the new slope is to be calculated
        :type x: float

        :param params_init: The initial values to parameterize the EC-AFC model.
        :type params_init: `settings.Initialization`

        :param params_eq: The equilibrated values to parameterize the EC-AFC model.
        :type params_eq: `equilibration.EquilibrationParams`

        :param params_sim: The parameters for the current simulation
        :type params_sim: dict[str: `Parameter`]
        """

        self.dy_dx = self.slopeFunc(slope_iter, x, params_init, params_eq,
                                    params_sim)

    def calcSlopeTrace(self, slope_iter, x, params_init, params_eq, params_sim,
                       trace, trace_iter):
        """
        Calculate the slope at point x using a particular function.

        :param slope_iter: The slope iteration
        :type slope_iter: int

        :param x: The x value for where the new slope is to be calculated
        :type x: float

        :param params_init: The initial values to parameterize the EC-AFC model.
        :type params_init: `settings.Initialization`

        :param params_eq: The equilibrated values to parameterize the EC-AFC
        model.
        :type params_eq: `equilibration.EquilibrationParams`

        :param params_sim: The parameters for the current simulation
        :type params_sim: dict[str: `Parameter`]

        :param trace: The parameters for the current trace element
        :type trace: dict[str: `Parameter`]

        :param trace: The iteration of the current trace element within the all
        traces elements.
        :type trace: int
        """

        self.dy_dx = self.slopeFuncTrace(slope_iter, x, params_init, params_eq,
                                         params_sim, trace, trace_iter)

    def addSlope(self, differential):
        """
        Calculate the slope at point x using a particular function and append to
        the slopes list.

        :param differential: The differential value to find the slope from
        :type differential: float
        """

        self.slopes.append(self.dy_dx * differential)

    def clearSlopes(self):
        """
        Clear the slopes from the previous iteration.
        """

        self.slopes = [0.0]

    def calcNewValue(self):
        """
        Calculate the new value of this parameter from the present value and its
        slopes.
        """

        k1, k2, k3, k4 = [self.slopes[k] for k in range(1, 5)]
        self.value_new = self.value_old - (k1 + 2 * k2 + 2 * k3 + k4) / 6

    def updateValue(self):
        """
        Replace the present value of this parameter with the new one.
        """

        self.value_old = self.value_new
        self.value_new = None
        self.str = f'{self.value_old:11.{self.decimals}f}'


class Parameter_Tm(Parameter):
    """
    The specific class to calculate the standing melt temperature.
    """
    def __init__(self, name, value_old, decimals=3, name_alt=None):
        """
        :param value_old: Initial value for this iteration
        :type value_old: float

        :param decimals: Number of decimals to print in the output string
        :type decimals: int
        """

        super().__init__(name, value_old, decimals, name_alt)

    def slopeFunc(self, slope_iter, x, params_init, params_eq, params_sim):
        """
        Calculate the slope at point x using a particular function.

        :param slope_iter: The slope iteration
        :type slope_iter: int

        :param x: The x value for where the new slope is to be calculated
        :type x: float

        :param params_init: The initial values to parameterize the EC-AFC model.
        :type params_init: `settings.Initialization`

        :param params_eq: The equilibrated values to parameterize the EC-AFC
        model.
        :type params_eq: `equilibration.EquilibrationParams`

        :param params_sim: The parameters for the current simulation
        :type params_sim: dict[str: `Parameter`]

        :return: the standing melt temperature
        :rtype: float
        """

        kwargs = {'dT': params_init.dT}
        return equil.Tm_slope(**kwargs)


class Parameter_Ta(Parameter):
    """
    The specific class to calculate the wall rock temperature.
    """
    def __init__(self, name, value_old, decimals=3, name_alt=None):
        """
        :param value_old: Initial value for this iteration
        :type value_old: float

        :param decimals: Number of decimals to print in the output string
        :type decimals: int
        """

        super().__init__(name, value_old, decimals, name_alt)

    def slopeFunc(self, slope_iter, x, params_init, params_eq, params_sim):
        """
        Calculate the slope at point x using a particular function.

        :param slope_iter: The slope iteration
        :type slope_iter: int

        :param x: The x value for where the new slope is to be calculated
        :type x: float

        :param params_init: The initial values to parameterize the EC-AFC model.
        :type params_init: `settings.Initialization`

        :param params_eq: The equilibrated values to parameterize the EC-AFC
        model.
        :type params_eq: `equilibration.EquilibrationParams`

        :param params_sim: The parameters for the current simulation
        :type params_sim: dict[str: `Parameter`]

        :return: the wall rock temperature
        :rtype: float
        """

        kwargs = {
            'Tm': x,
            'Ta': self.thisStepValue(slope_iter),
            'Ma0': params_eq.Ma0
        }
        return equil.conserv_energy_dTa_dTm(**kwargs)


class Parameter_Mm(Parameter):
    """
    The specific class to calculate the magma body mass.
    """
    def __init__(self, name, value_old, decimals=3, name_alt=None):
        """
        :param value_old: Initial value for this iteration
        :type value_old: float

        :param decimals: Number of decimals to print in the output string
        :type decimals: int
        """

        super().__init__(name, value_old, decimals, name_alt)

    def slopeFunc(self, slope_iter, x, params_init, params_eq, params_sim):
        """
        Calculate the slope at point x using a particular function.

        :param slope_iter: The slope iteration
        :type slope_iter: int

        :param x: The x value for where the new slope is to be calculated
        :type x: float

        :param params_init: The initial values to parameterize the EC-AFC model.
        :type params_init: `settings.Initialization`

        :param params_eq: The equilibrated values to parameterize the EC-AFC
        model.
        :type params_eq: `equilibration.EquilibrationParams`

        :param params_sim: The parameters for the current simulation
        :type params_sim: dict[str: `Parameter`]

        :return: the magma body mass
        :rtype: float
        """

        kwargs = {
            'Tm': x,
            'Ta': params_sim['Ta'].thisStepValue(slope_iter),
            'Ma0': params_eq.Ma0
        }
        return equil.conserv_mass_dMm_dTm(**kwargs)


class Parameter_dm(Parameter):
    """
    The specific class to calculate the oxygen isotopic composition of standing
    magma
    """
    def __init__(self, name, value_old, decimals=3, name_alt=None):
        """
        :param value_old: Initial value for this iteration
        :type value_old: float

        :param decimals: Number of decimals to print in the output string
        :type decimals: int
        """

        super().__init__(name, value_old, decimals, name_alt)

    def slopeFunc(self, slope_iter, x, params_init, params_eq, params_sim):
        """
        Calculate the slope at point x using a particular function.

        :param slope_iter: The slope iteration
        :type slope_iter: int

        :param x: The x value for where the new slope is to be calculated
        :type x: float

        :param params_init: The initial values to parameterize the EC-AFC model.
        :type params_init: `settings.Initialization`

        :param params_eq: The equilibrated values to parameterize the EC-AFC
        model.
        :type params_eq: `equilibration.EquilibrationParams`

        :param params_sim: The parameters for the current simulation
        :type params_sim: dict[str: `Parameter`]

        :return: the oxygen isotopic composition of standing magma
        :rtype: float
        """

        kwargs = {
            'Tm': x,
            'Ta': params_sim['Ta'].thisStepValue(slope_iter),
            'Ma0': params_eq.Ma0,
            'Mm': params_sim['Mm'].thisStepValue(slope_iter),
            'dTa_dTm': params_sim['Ta'].dy_dx,
            'dm': self.thisStepValue(slope_iter),
        }
        return equil.oxygen_isotope_comp_ddm_dTm(**kwargs)


class Parameter_Cm(Parameter):
    """
    The specific class to calculate the conc. of trace element in standing melt.
    """
    def __init__(self, name, value_old, decimals=3, name_alt=None):
        """
        :param value_old: Initial value for this iteration
        :type value_old: float

        :param decimals: Number of decimals to print in the output string
        :type decimals: int
        """

        super().__init__(name, value_old, decimals, name_alt)

    def slopeFuncTrace(self, slope_iter, x, params_init, params_eq, params_sim,
                       param_trace, trace_iter):
        """
        Calculate the slope at point x using a particular function.

        :param slope_iter: The slope iteration
        :type slope_iter: int

        :param x: The x value for where the new slope is to be calculated
        :type x: float

        :param params_init: The initial values to parameterize the EC-AFC model.
        :type params_init: `settings.Initialization`

        :param params_eq: The equilibrated values to parameterize the EC-AFC
        model.
        :type params_eq: `equilibration.EquilibrationParams`

        :param params_sim: The parameters for the current simulation
        :type params_sim: dict[str: `Parameter`]

        :param param_trace: The parameters for the current trace element
        :type param_trace: dict[str: `Parameter`]

        :param trace: The iteration of the current trace element within the all
        traces elements.
        :type trace: int
        """

        kwargs = {
            'Tm': x,
            'Ta': params_sim['Ta'].thisStepValue(slope_iter),
            'Ma0': params_eq.Ma0,
            'Mm': params_sim['Mm'].thisStepValue(slope_iter),
            'Cm': self.thisStepValue(slope_iter),
            'Ca0': params_init.traces[trace_iter]['Ca0'],
            'Cm0': params_init.traces[trace_iter]['Cm0'],
            'dTa_dTm': params_sim['Ta'].dy_dx,
        }
        return equil.conc_trace_elem_dCm_dTm(**kwargs)


class Parameter_em(Parameter):
    """
    The specific class to calculate the isotopic ratio of trace element in
    standing melt.
    """
    def __init__(self, name, value_old, decimals=3, name_alt=None):
        """
        :param value_old: Initial value for this iteration
        :type value_old: float

        :param decimals: Number of decimals to print in the output string
        :type decimals: int
        """

        super().__init__(name, value_old, decimals, name_alt)

    def slopeFuncTrace(self, slope_iter, x, params_init, params_eq, params_sim,
                       param_trace, trace_iter):
        """
        Calculate the slope at point x using a particular function.

        :param slope_iter: The slope iteration
        :type slope_iter: int

        :param x: The x value for where the new slope is to be calculated
        :type x: float

        :param params_init: The initial values to parameterize the EC-AFC model.
        :type params_init: `settings.Initialization`

        :param params_eq: The equilibrated values to parameterize the EC-AFC
        model.
        :type params_eq: `equilibration.EquilibrationParams`

        :param params_sim: The parameters for the current simulation
        :type params_sim: dict[str: `Parameter`]

        :param param_trace: The parameters for the current trace element
        :type param_trace: dict[str: `Parameter`]

        :param trace: The iteration of the current trace element within the all
        traces elements.
        :type trace: int
        """

        kwargs = {
            'Tm': x,
            'Ta': params_sim['Ta'].thisStepValue(slope_iter),
            'Ma0': params_eq.Ma0,
            'Mm': params_sim['Mm'].thisStepValue(slope_iter),
            'Cm': param_trace['Cm'].thisStepValue(slope_iter),
            'Ca0': params_init.traces[trace_iter]['Ca0'],
            'Cm0': params_init.traces[trace_iter]['Cm0'],
            'dTa_dTm': params_sim['Ta'].dy_dx,
            'ea0': params_init.traces[trace_iter]['ea0'],
            'em': self.thisStepValue(slope_iter),
        }
        return equil.isotop_ratio_dem_dTm(**kwargs)


class ECAFC:
    """
    This is the main EC-AFC simulation routine. The important calculated
    quantities in `self.params` are updated each iteration of the simulation and
    stored by corresponding keys in `self.results`.
    """
    def __init__(self, params_init):
        """
        The EC-AFC simulation is initialized by passing the initialized
        parameters object. Then equilibrate the quantities for starting values
        at the requested normalized equilibration temperature (Teq_norm).

        :param params_init: The initial values to parameterize the EC-AFC model.
        :type params_init: `Initialization`
        """

        self.params_init = params_init
        self.params_eq = equil.EquilibrationParams(self.params_init.Teq_norm)
        self.params_sim = {
            'Tm':
            Parameter_Tm('Tm,norm',
                         params_init.Tm0_norm,
                         decimals=3,
                         name_alt='Tm'),
            'Ta':
            Parameter_Ta('Ta,norm',
                         params_init.Ta0_norm,
                         decimals=3,
                         name_alt='Ta'),
            'Mm':
            Parameter_Mm('Mm', params_init.Mm0, decimals=3),
            'dm':
            Parameter_dm('dm', params_init.dm0, decimals=3),
        }
        self.params_traces = []
        if params_init.traces:
            for trace in params_init.traces:
                params_trace = {
                    'elem': trace['elem'],
                    'isoratio': trace['isoratio'],
                    'Cm': Parameter_Cm('Cm', trace['Cm0'], decimals=3),
                    'em': Parameter_em('em', trace['em0'], decimals=3),
                }
                self.params_traces.append(params_trace)
        self.results = results.Results()

    def simulate(self):
        """
        Solve the system of nonlinear equations for assimilation and fractional
        crystalization using 4th-order Runge-Kutta (RK4) for each normalized
        temperature step in the requested temperature range as the melting temp
        (Tnorm1) cools to the equilibration temperature (Tnorm0).

        :param print_lines: Number of results lines to print (0 to not print) at
        the beginning and end of the simulation
        :type print_lines: int
        """

        dT = self.params_init.dT
        Tnorm0 = self.params_init.Teq_norm - dT
        Tnorm1 = self.params_init.Tm0_norm

        for i, T in enumerate(np.arange(Tnorm1, Tnorm0, -dT), start=1):
            self.params_sim['Tm'].value_old = T

            self.results.store(copy.deepcopy(self.params_sim),
                               copy.deepcopy(self.params_traces))

            for param in self.params_sim.values():
                param.clearSlopes()
            for trace in self.params_traces:
                for name in ['Cm', 'em']:
                    param_trace = trace[name]
                    param_trace.clearSlopes()

            for h in range(4):
                T = get_T_iter(self.params_sim['Tm'].value_old, dT, h)
                for param in self.params_sim.values():
                    param.calcSlope(h, T, self.params_init, self.params_eq,
                                    self.params_sim)
                    param.addSlope(dT)
                for t, trace in enumerate(self.params_traces):
                    for name in ['Cm', 'em']:
                        param_trace = trace[name]
                        param_trace.calcSlopeTrace(h, T, self.params_init,
                                                   self.params_eq,
                                                   self.params_sim, trace, t)
                        param_trace.addSlope(dT)

            for name, param in self.params_sim.items():
                param.calcNewValue()
                param.updateValue()
            for trace in self.params_traces:
                for name in ['Cm', 'em']:
                    param_trace = trace[name]
                    param_trace.calcNewValue()
                    param_trace.updateValue()


"""
Below is an example usage of the API.
"""
if __name__ == '__main__':
    import dataclasses
    import os

    from petrosim.models.ecafc import loader as ldr


    ldr.init(os.path.join(os.path.dirname(ldr.__file__), 'example.in'))
    params_init = ldr.Parameters(**dataclasses.asdict(ldr.params))

    ecafc = ECAFC(params_init)
    ecafc.simulate()
    ecafc.results.print()
    #ecafc.results.write('out.csv')
