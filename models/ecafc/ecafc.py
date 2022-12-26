"""
This is the main driver for the Energy-Constrained Assimilation Fractional
Crystallization (EC-AFC) model based on work by Spera & Bohrson J. Petrol., 42,
999–1018, 2001.
"""

import dataclasses

import numpy as np

#import constants_old as const
import equilibration as equil
import settings


@dataclasses.dataclass
class Params:
    """
    Keep the calculated quantities for each iteration here.
    """

    Tm: float = None
    Ta: float = None
    Mm: float = None
    Cm: float = None
    em: float = None
    dm: float = None


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
        self.params = {
            'Tm': 1.0,
            'Ta': params_init.Ta0_norm,
            'Mm': params_init.Mm0,
            'Cm': params_init.Cm0,
            'em': params_init.em0,
            'dm': params_init.dm0,
        }
        self.results = {
            'Tm_norm': [],
            'Tm': [],
            'Ta_norm': [],
            'Ta': [],
            'Mm': [],
            'Cm': [],
            'em': [],
            'dm': [],
        }

    def _clearSlopes(self):
        """
        Initialize the slopes for each RK iteration.
        """

        self.slopes = {
            'Ta': [0.0],
            'Mm': [0.0],
            'Cm': [0.0],
            'em': [0.0],
            'dm': [0.0],
        }

    def _nextStep(self, param, slope_iter):
        """
        Find the next parameter value from its previous value and slope.

        :param param: The parameter of interest to pull its prior value and slope.
        :type param: str

        :param slope_iter: The slope iteration.
        :type slope_iter: int
        """

        if param == 'Tm':
            return
        if slope_iter == 3:
            return self.params[param] - self.slopes[param][-1]
        else:
            return self.params[param] - self.slopes[param][-1] / 2

    def _addSlope(self, param, differential):
        """
        Calculate and store the new slope for a new interval.

        :param param: The parameter of interest to add the new slope to.
        :type param: str

        :param differential: The differential at an interval
        :type differential: float
        """

        slope = self.params_init.dT * differential
        self.slopes[param].append(slope)

    def _calcSlopes(self, Tm, slope_iter):
        """
        A 4th-order Runge-Kutta (RK4) iteration at a new interval.

        :param Tm: Melting temperature (normalized)
        :type Tm: float

        :param slope_iter: The RK4 slope iteration
        :type slope_iter: int
        """

        if slope_iter == 0:
            Tm = Tm
        elif slope_iter == 3:
            Tm = Tm - self.params_init.dT
        else:
            Tm = Tm - self.params_init.dT / 2

        params = Params(*[self._nextStep(p, slope_iter) for p in self.params])

        dTa_dTm = equil.conserv_energy_dTa_dTm(Tm, params.Ta,
                                               self.params_eq.Ma0)
        dMm_dTm = equil.conserv_mass_dMm_dTm(Tm, params.Ta, self.params_eq.Ma0)
        dCm_dTm = equil.conc_trace_elem_dCm_dTm(Tm, params.Ta,
                                                self.params_eq.Ma0, params.Mm,
                                                params.Cm, dTa_dTm)
        dem_dTm = equil.isotop_ratio_dem_dTm(Tm, params.Ta, self.params_eq.Ma0,
                                             params.Mm, params.Cm, dTa_dTm,
                                             params.em)
        ddm_dTm = equil.oxygen_isotope_comp_ddm_dTm(Tm, params.Ta,
                                                    self.params_eq.Ma0,
                                                    params.Mm, dTa_dTm,
                                                    params.dm)

        self._addSlope('Ta', dTa_dTm)
        self._addSlope('Mm', dMm_dTm)
        self._addSlope('Cm', dCm_dTm)
        self._addSlope('em', dem_dTm)
        self._addSlope('dm', ddm_dTm)

    def _updateValue(self, param):
        """
        Update the stored parameter.

        :param param: The parameter to update
        :type param: str
        """

        old = self.params[param]
        k1, k2, k3, k4 = [self.slopes[param][k] for k in range(1, 5)]
        new = old - (k1 + 2 * k2 + 2 * k3 + k4) / 6
        self.params[param] = new

    def simulate(self, print_lines=3):
        """
        Solve the system of nonlinear equations for assimilation and fractional
        crystalization using a 4th-order Runge-Kutta (RK4) for each normalized
        temperature step in the requested temperature range as the melting temp
        (Tnorm1) cools to the equilibration temperature (Tnorm0).

        :param Tm: Melting temperature (normalized)
        :type Tm: float

        :param slope_iter: The RK4 slope iteration
        :type slope_iter: int
        """

        dT = self.params_init.dT
        Tnorm0 = self.params_init.Teq_norm - dT
        Tnorm1 = self.params_init.T1
        max_iter = (Tnorm1 - Tnorm0) / dT

        if print_lines:
            self._printResults(0, max_iter)
        for i, T in enumerate(np.arange(Tnorm1, Tnorm0, -dT), start=1):
            self.params['Tm'] = T

            if print_lines:
                self._printResults(i, max_iter, lines_shown=print_lines)
            self._saveResults()

            self._clearSlopes()
            for h in range(4):
                self._calcSlopes(T, h)

            for param in self.params.keys():
                if param == 'Tm':
                    continue
                self._updateValue(param)

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
            print(f'{"Tm,norm":>8s}', f'{"Tm":>8s}', f'{"Ta,norm":>8s}',
                  f'{"Tm":>8s}', f'{"Mm":>8s}', f'{"Cm":>8s}',
                  f'{"em":>8s}')  #, f'{self.params["dm"]:>8s}')
            return

        T_m = equil.unnormalize_temp(self.params['Tm'])
        T_a = equil.unnormalize_temp(self.params['Ta'])

        if lines_shown == -1:
            print(f'{self.params["Tm"]:8.4f}', f'{T_m:8.2f}',
                  f'{self.params["Ta"]:8.4f}', f'{T_a:8.2f}',
                  f'{self.params["Mm"]:8.3f}', f'{self.params["Cm"]:8.3f}',
                  f'{self.params["em"]:8.3f}')  #, f'{self.params["dm"]:7.3f}')
        else:
            if lines_shown < iter <= max_iter - lines_shown:
                if iter == lines_shown + 1:
                    print(f'{"":.>3s}')
            else:
                print(f'{self.params["Tm"]:8.4f}', f'{T_m:8.2f}',
                      f'{self.params["Ta"]:8.4f}', f'{T_a:8.2f}',
                      f'{self.params["Mm"]:8.3f}', f'{self.params["Cm"]:8.3f}',
                      f'{self.params["em"]:8.3f}'
                      )  #, f'{self.params["dm"]:7.3f}')

    def _saveResults(self):
        """
        Save the results to a dict of lists keyed by the parameter names.
        """

        self.results['Tm_norm'].append(self.params['Tm'])
        self.results['Ta_norm'].append(self.params['Ta'])
        for param in self.params:
            if param in ['Ta', 'Tm']:
                self.results[param].append(equil.unnormalize_temp(
                    self.params[param]))
            else:
                self.results[param].append(self.params[param])


class Initialization:
    """
    Helper class to read in, validate, and initialize the required parameters.
    """

    def __init__(self, **kwargs):
        """
        Pass in the required starting values for the simulation's parameters as
        a dict. Empty initialization will display which parameters are required.
        """

        self.__dict__.update(kwargs)
        self._validate_args()

    def _validate_args(self):
        """
        Validate the input parameters and make sure none are missing.
        """

        args_reqd = [
            'T1',
            'Teq_norm',
            'Ta0_norm',
            'dT',
            'Mm0',
            'Cm0',
            'Da0',
            'Dm0',
            'dHa',
            'dHm',
            'Ca0',
            'Cm0',
            'ea0',
            'em0',
            'da0',
            'dm0',
            'koxy',
        ]
        args_missing = []
        for arg in args_reqd:
            if not arg in self.__dict__:
                args_missing.append(arg)
        if len(args_missing) > 0:
            exit(
                f'Please initialize the following parameters: {", ".join(args_missing)}'
            )



"""
Below is an example usage of the API.
"""
if __name__ == '__main__':
    settings.init('example.in')
    params_init = settings.Initialization(**dataclasses.asdict(settings.const))
    #settings.InputFileReader('example.in')

#    params_init = Initialization(
#        T1=1,
#        Teq_norm=0.88,
#        Ta0_norm=settings.const.Ta0 / settings.const.Tlm,
#        dT=0.001,
#        Mm0=1,
#        Da0=1.5,
#        Dm0=1.5,
#        dHa=0,
#        dHm=0,
#        Ca0=350,
#        Cm0=700,
#        ea0=0.722,
#        em0=0.7035,
#        da0=1,
#        dm0=1,
#        koxy=1,
#    )

    ecafc = ECAFC(params_init)
    ecafc.simulate(print_lines=4)
