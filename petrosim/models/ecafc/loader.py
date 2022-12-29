import dataclasses
import re

from ruamel.yaml import YAML
from voluptuous import All, Any, Coerce, Optional, Range, Required, Schema


BULK = 'BULK'
MAGMA = 'MAGMA'
ASSIMILANT = 'ASSIMILANT'
OTHER = 'OTHER'
TRACE = 'TRACE'

K = 273.16          # Conversion factor from ºC to K
R = 8.31432         # Gas constant

SCHEMA_COMPONENT_BULK = Schema({
    Required('M0'): Coerce(float),
    Required('T0'): Coerce(float),
    Required('Tl'): Coerce(float),
    Required('cp'): Coerce(float),
    Required('dh'): Coerce(float),
    Required('loga'): Coerce(float),
    Required('logb'): Coerce(float),
    Required('D0'): Coerce(float),
    Required('dH'): Coerce(float),
    Required('d0'): Coerce(float),
})
SCHEMA_OTHER = Schema({
    Required('Ts'): Coerce(float),
    Required('dT'): Coerce(float),
    Required('Teq'): Any(float, int, str, None),
    Required('Teq_norm'): Any(float, int, str, None),
    Required('koxy'): Coerce(float),
})
SCHEMA_BULK = Schema({
    Required(MAGMA): SCHEMA_COMPONENT_BULK,
    Required(ASSIMILANT): SCHEMA_COMPONENT_BULK,
    Required(OTHER): SCHEMA_OTHER,
})
SCHEMA_COMPONENT_TRACE = Schema({
    Required('C0'): Coerce(float),
    Required('e0'): Coerce(float),
})
SCHEMA_TRACE_ENTRY = Schema({
    Required('elem'): str,
    Required('isoratio'): str,
    Required(MAGMA): SCHEMA_COMPONENT_TRACE,
    Required(ASSIMILANT): SCHEMA_COMPONENT_TRACE,
})
SCHEMA_TRACE = Schema([SCHEMA_TRACE_ENTRY])
SCHEMA_INPUT = Schema({
    Required(BULK): SCHEMA_BULK,
    Optional(TRACE): SCHEMA_TRACE,
})

DEST_DICT_BULK = {
    'M0': {MAGMA: 'Mm0', ASSIMILANT: 'Ma0'},
    'T0': {MAGMA: 'Tm0', ASSIMILANT: 'Ta0'},
    'Tl': {MAGMA: 'Tlm', ASSIMILANT: 'Tla'},
    'cp': {MAGMA: 'cpm', ASSIMILANT: 'cpa'},
    'dh': {MAGMA: 'dhm', ASSIMILANT: 'dha'},
    'loga': {MAGMA: 'loga_m', ASSIMILANT: 'loga_a'},
    'logb': {MAGMA: 'logb_m', ASSIMILANT: 'logb_a'},
    'D0': {MAGMA: 'Dm0', ASSIMILANT: 'Da0'},
    'dH': {MAGMA: 'dHm', ASSIMILANT: 'dHa'},
    'd0': {MAGMA: 'dm0', ASSIMILANT: 'da0'},
}
DEST_DICT_TRACE = {
    'C0': {MAGMA: 'Cm0', ASSIMILANT: 'Ca0'},
    'e0': {MAGMA: 'em0', ASSIMILANT: 'ea0'},
}


def _read_input_file(input_file):
    """
    Read and parse the input YAML file.
    """

    with open(input_file) as filehandle:
        yaml = YAML()
        return _validate_input_data(yaml.load(filehandle.read()))

def _validate_input_data(input_data):
    """
    Validate the input data.
    """

    return SCHEMA_INPUT(input_data)


class ParametersReader:
    """
    Class to parse the input YAML file and assign the parameters to variables.
    """

    def __init__(self, input_file):
        """
        :param input_file: The name of the input file in YAML format.
        :type input_file: str
        """

        self.parameters = Parameters()
        self.initial_values = _read_input_file(input_file)
        self._assignParameters()

    def _assignParameters(self):
        for param_fam, dest in DEST_DICT_BULK.items():
            for component, param in dest.items():
                setattr(self.parameters, param,
                        self.initial_values[BULK][component][param_fam])

        self.parameters.Ts = self.initial_values[BULK][OTHER]['Ts']
        self.parameters.dT = self.initial_values[BULK][OTHER]['dT']
        self.parameters.koxy = self.initial_values[BULK][OTHER]['koxy']

        Teq = self.initial_values[BULK][OTHER]['Teq']
        Teq_norm = self.initial_values[BULK][OTHER]['Teq_norm']
        if all([Teq not in ('None', None), Teq_norm not in ('None', None)]):
            msg = 'Cannot have both Teq and Teq_norm set at the same time. Either remove one value or set one to "None" and then resubmit.'
            exit(msg)
        if Teq not in ('None', None):
            self.parameters.Teq_norm = float(
                self.initial_values[BULK][OTHER]['Teq']) / self.parameters.Tlm
        else:
            self.parameters.Teq_norm = float(
                self.initial_values[BULK][OTHER]['Teq_norm'])

        for param in ['Tlm', 'Tla', 'Tm0', 'Ta0', 'Ts']:
            uncorrected_temp = getattr(self.parameters, param)
            setattr(self.parameters, param,
                    uncorrected_temp + self.parameters.K)

        self.parameters.Tm0_norm = self.parameters.Tm0 / self.parameters.Tlm
        self.parameters.Ta0_norm = self.parameters.Ta0 / self.parameters.Tlm

        values_traces = self.initial_values.get(TRACE)
        if values_traces:
            for entry in values_traces:
                parameters_trace = ParametersTrace()
                parameters_trace.elem = entry['elem']
                parameters_trace.isoratio = entry['isoratio']
                for param_fam, dest in DEST_DICT_TRACE.items():
                    for component, param in dest.items():
                        setattr(parameters_trace, param,
                                entry[component][param_fam])
                self.parameters.traces.append(parameters_trace)


@dataclasses.dataclass
class ParametersTrace:
    elem: str = None
    isoratio: str = None
    Ca0: float = None
    Cm0: float = None
    ea0: float = None
    em0: float = None


@dataclasses.dataclass
class Parameters:
    K: float = K
    R: float = R
    Mm0: float = None
    Ma0: float = None
    Ts: float = None
    T1: float = 1.0
    Ta0: float = None
    Tla: float = None
    Tlm: float = None
    Teq_norm: float = None
    Ta0_norm: float = None
    Tm0_norm: float = None
    dT: float = None
    Mm0: float = None
    Da0: float = None
    Dm0: float = None
    dHa: float = None
    dHm: float = None
    traces: list[ParametersTrace] = dataclasses.field(default_factory=list)
    da0: float = None
    dm0: float = None
    koxy: float = None


def init(INFILE):
    global params
    params = ParametersReader(INFILE).parameters


if __name__ == '__main__':
    import os

    INFILE = os.path.join(os.path.dirname(__file__), 'example.in')
    params = ParametersReader(INFILE)
    print(params.parameters)
