import dataclasses
import re

ASSIMILANT = 'assimilant'
MAGMA = 'magma'
OTHER = 'other'

K = 273.16          # Conversion factor from ºC to K
R = 8.31432         # Gas constant


def _read_param_line(line):
    if re.match(r'^\s+', line):
        if len(line.split()) >= 2:
            param, value = line.strip().split()[:2]
            return param, value
        else:
            return
    else:
        return


class InputFileReader:

    def __init__(self, input_file):
        self.params_init = {MAGMA.upper(): {}, ASSIMILANT.upper(): {}, OTHER.upper(): {}}
        self.initial_values = self._readInputFile(input_file)


    def _readInputFile(self, input_file):
        with open(input_file) as opened:
            self._splitBlocks(opened)
            self._readParameters()
            self._assignParameters()

    def _splitBlocks(self, file_handle):
        self.blocks = {MAGMA.upper(): [], ASSIMILANT.upper(): [], OTHER.upper(): []}
        block = None
        for line in file_handle.readlines():
            if len(line.strip()) == 0:
                block = None
                continue
            elif line.strip() in [MAGMA.upper(), ASSIMILANT.upper(), OTHER.upper()]:
                block = line.strip().upper()
            else:
                self.blocks[block].append(line)

    def _readParameters(self):
        for block in [MAGMA.upper(), ASSIMILANT.upper(), OTHER.upper()]:
            for line in self.blocks[block]:
                data = _read_param_line(line)
                if data:
                    param, value = data
                    if value.lstrip("-")[0].isdigit():
                        self.params_init[block][param] = float(value)
                    else:
                        print(block, line.strip())

    def _assignParameters(self):
        self.constants = Constants()
        dest_dict = {
            'T0': {MAGMA: 'Tm0', ASSIMILANT: 'Ta0'},
            'Tl': {MAGMA: 'Tlm', ASSIMILANT: 'Tla'},
            'cp': {MAGMA: 'cpm', ASSIMILANT: 'cpa'},
            'dh': {MAGMA: 'dhm', ASSIMILANT: 'dha'},
            'alpha': {MAGMA: 'alpha_m', ASSIMILANT: 'alpha_a'},
            'beta': {MAGMA: 'beta_m', ASSIMILANT: 'beta_a'},
            'D0': {MAGMA: 'Dm0', ASSIMILANT: 'Da0'},
            'dH': {MAGMA: 'dHm', ASSIMILANT: 'dHa'},
            'C0': {MAGMA: 'Cm0', ASSIMILANT: 'Ca0'},
            'e0': {MAGMA: 'em0', ASSIMILANT: 'ea0'},
            'd0': {MAGMA: 'dm0', ASSIMILANT: 'da0'},
        }
        for param, dest in dest_dict.items():
            for block, variant in dest.items():
                setattr(self.constants, variant, self.params_init[block.upper()][param])

        self.constants.Ts = self.params_init[OTHER.upper()]['Ts']
        self.constants.dT = self.params_init[OTHER.upper()]['dT']
        self.constants.Mm0 = self.params_init[OTHER.upper()]['Mm0']
        self.constants.koxy = self.params_init[OTHER.upper()]['koxy']
        self.constants.Teq_norm = self.params_init[OTHER.upper()]['Teq']
        self.constants.Ta0_norm = self.constants.Ta0 / self.constants.Tlm

        for param in ['Tlm', 'Tla', 'Tm0', 'Ta0', 'Ts']:
            uncorrected_temp = getattr(self.constants, param)
            setattr(self.constants, param, uncorrected_temp + self.constants.K)
        self.constants.T1 = 1


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


@dataclasses.dataclass
class Constants:
        K: float = K
        R: float = R
        Ts: float = None
        T1: float = None
        Ta0: float = None
        Tla: float = None
        Tlm: float = None
        Teq_norm: float = None
        Ta0_norm: float = None
        dT: float = None
        Mm0: float = None
        Cm0: float = None
        Da0: float = None
        Dm0: float = None
        dHa: float = None
        dHm: float = None
        Ca0: float = None
        Cm0: float = None
        ea0: float = None
        em0: float = None
        da0: float = None
        dm0: float = None
        koxy: float = None


def init(INFILE):
    global const
    const = InputFileReader(INFILE).constants


if __name__ == '__main__':
    INFILE = 'example.in'
    InputFileReader(INFILE)
