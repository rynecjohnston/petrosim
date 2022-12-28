"""
This is the main driver for the Energy-Constrained Assimilation Fractional
Crystallization (EC-AFC) model based on work by Spera & Bohrson J. Petrol., 42,
999–1018, 2001.
"""

import numpy as np

import equilibration as equil


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
    def __init__(self, value_old, decimals=3):
        """
        :param value_old: Initial value for this iteration
        :type value_old: float

        :param decimals: Number of decimals to print in the output string
        :type decimals: int
        """

        self.value_old = value_old
        self.decimals = decimals
        self.str = f'{self.value_old:11.{self.decimals}f}'
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
    def __init__(self, value_old, decimals):
        """
        :param value_old: Initial value for this iteration
        :type value_old: float

        :param decimals: Number of decimals to print in the output string
        :type decimals: int
        """

        super().__init__(value_old, decimals)

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
    def __init__(self, value_old, decimals):
        """
        :param value_old: Initial value for this iteration
        :type value_old: float

        :param decimals: Number of decimals to print in the output string
        :type decimals: int
        """

        super().__init__(value_old, decimals)

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
    def __init__(self, value_old, decimals):
        """
        :param value_old: Initial value for this iteration
        :type value_old: float

        :param decimals: Number of decimals to print in the output string
        :type decimals: int
        """

        super().__init__(value_old, decimals)

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
    def __init__(self, value_old, decimals):
        """
        :param value_old: Initial value for this iteration
        :type value_old: float

        :param decimals: Number of decimals to print in the output string
        :type decimals: int
        """

        super().__init__(value_old, decimals)

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
    def __init__(self, value_old, decimals):
        """
        :param value_old: Initial value for this iteration
        :type value_old: float

        :param decimals: Number of decimals to print in the output string
        :type decimals: int
        """

        super().__init__(value_old, decimals)

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
    def __init__(self, value_old, decimals):
        """
        :param value_old: Initial value for this iteration
        :type value_old: float

        :param decimals: Number of decimals to print in the output string
        :type decimals: int
        """

        super().__init__(value_old, decimals)

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
            'Tm': Parameter_Tm(params_init.Tm0_norm, decimals=4),
            'Ta': Parameter_Ta(params_init.Ta0_norm, decimals=4),
            'Mm': Parameter_Mm(params_init.Mm0, decimals=3),
            'dm': Parameter_dm(params_init.dm0, decimals=3),
        }
        self.params_traces = []
        for trace in params_init.traces:
            params_trace = {
                'elem': trace['elem'],
                'isoratio': trace['isoratio'],
                'Cm': Parameter_Cm(trace['Cm0'], decimals=3),
                'em': Parameter_em(trace['em0'], decimals=3),
            }
            self.params_traces.append(params_trace)
        self.results = []

    def simulate(self, print_lines=3):
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
        max_iter = (Tnorm1 - Tnorm0) / dT

        if print_lines:
            self._printResults(0, max_iter)
        self._storeResults(0)
        for i, T in enumerate(np.arange(Tnorm1, Tnorm0, -dT), start=1):
            self.params_sim['Tm'].value_old = T

            if print_lines:
                self._printResults(i, max_iter, lines_shown=print_lines)
            self._storeResults(i)

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

    def _printResults(self, iter, max_iter, lines_shown=3):
        """
        Print the results for the beginning and end of the simulation, and
        truncate the middle. 

        :param iter: Iteration number
        :type iter: int

        :param max_iter: Maximum number of iterations
        :type max_iter: int

        :param lines_shown: Number of beginning and final lines to print out
        :type lines_shown: int
        """

        if iter == 0:
            fields = [f'{"":11s}', f'{"":11s}', f'{"":11s}', f'{"":11s}',
            f'{"":11s}', f'{"":11s}']
            for trace in self.params_traces:
                for name in ['elem', 'isoratio']:
                    fields.append(f'{trace[name]:>11s}')
            print(' '.join(fields))

            fields = [
                f'{"Tm,norm":>11s}', f'{"Tm":>11s}', f'{"Ta,norm":>11s}',
                f'{"Ta":>11s}', f'{"Mm":>11s}', f'{"dm":>11s}'
            ]
            for trace in self.params_traces:
                for name in ['Cm', 'em']:
                    fields.append(f'{name:>11s}')
            print(' '.join(fields))
            return

        T_m = equil.unnormalize_temp(self.params_sim['Tm'].value_old)
        T_a = equil.unnormalize_temp(self.params_sim['Ta'].value_old)

        if lines_shown == -1:
            values = [
                self.params_sim["Tm"].str, f'{T_m:11.2f}',
                self.params_sim["Ta"].str, f'{T_a:11.2f}',
                self.params_sim["Mm"].str, self.params_sim["dm"].str
            ]
            for trace in self.params_traces:
                for name in ['Cm', 'em']:
                    param_trace = trace[name]
                    values.append(param_trace.str)
            print(' '.join(values))
        else:
            if lines_shown < iter <= max_iter - lines_shown + 1:
                if iter == lines_shown + 1:
                    print(f'{"":.>3s}')
            else:
                values = [
                    self.params_sim["Tm"].str, f'{T_m:11.2f}',
                    self.params_sim["Ta"].str, f'{T_a:11.2f}',
                    self.params_sim["Mm"].str, self.params_sim["dm"].str
                ]
                for trace in self.params_traces:
                    for name in ['Cm', 'em']:
                        param_trace = trace[name]
                        values.append(param_trace.str)
                print(' '.join(values))

    def _storeResults(self, iter):
        """
        Store the results.

        :param iter: Iteration number
        :type iter: int
        """

        if iter == 0:
            fields = ["", "", "", "", "", ""]
            for trace in self.params_traces:
                for name in ['elem', 'isoratio']:
                    fields.append(trace[name])
            self.results.append(fields)
            fields = ["Tm_norm", "Tm", "Ta_norm", "Ta", "Mm", "dm"]
            for trace in self.params_traces:
                for name in ['Cm', 'em']:
                    fields.append(f'{name}')
            self.results.append(fields)
            return

        T_m = equil.unnormalize_temp(self.params_sim['Tm'].value_old)
        T_a = equil.unnormalize_temp(self.params_sim['Ta'].value_old)

        values = [
            self.params_sim["Tm"].value_old, f'{T_m:.2f}',
            self.params_sim["Ta"].value_old, f'{T_a:.2f}',
            self.params_sim["Mm"].value_old, self.params_sim["dm"].value_old
        ]
        for trace in self.params_traces:
            for name in ['Cm', 'em']:
                param_trace = trace[name]
                values.append(param_trace.value_old)
        self.results.append(values)


"""
Below is an example usage of the API.
"""
if __name__ == '__main__':
    import dataclasses

    import loader as ldr


    ldr.init('example.in')
    params_init = ldr.Parameters(**dataclasses.asdict(ldr.params))

    ecafc = ECAFC(params_init)
    ecafc.simulate(print_lines=4)
