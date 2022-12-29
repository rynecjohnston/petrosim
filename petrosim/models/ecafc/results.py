"""
This module collects the results of the EC-AFC simulation and enables printing
to screen and writing to file.
"""


import csv

from petrosim.models.ecafc import equilibration as equil


MAXLEN = 8
PARAMETERS_TRACE = {'Cm': 'elem', 'em': 'isoratio'}


class Results:
    """
    Class to store simulation results and methods to print or write them to file
    """
    def __init__(self):
        self.results = []

    def store(self, results_bulk, results_traces=None):
        """
        Store the results for each iteration in the `self.results` list

        :param results_bulk: Dictionary of the parameters for the bulk
        :type results_bulk: dict[str: `ecafc.Parameter`]

        :param results_traces: List of dictionary of the parameters for the
        trace elements
        :type results_traces: list[dict[str: `ecafc.Parameter`]]
        """

        if results_traces:
            self.results.append([results_bulk, results_traces])
        else:
            self.results.append([results_bulk])

    def print(self, lines_shown=3):
        """
        Print the results for the beginning & end of the simulation and
        truncate the middle.

        :param lines_shown: Number of beginning and final lines to print out
        :type lines_shown: int
        """

        print_lines = self._getHeader()
        for i, row in enumerate(self.results):
            print_line = []
            for name, param in row[0].items():
                maxlen = self._getMaxStrLength(name)
                print_line.append(
                    f'{param.value_old:>{maxlen}.{param.decimals}f}')
                if param.name_alt:
                    value_alt = equil.unnormalize_temp(param.value_old)
                    print_line.append(
                        f'{value_alt:>{maxlen}.{param.decimals}f}')
            if len(row) == 2:
                for j, trace in enumerate(row[1], start=1):
                    for param_name in PARAMETERS_TRACE:
                        param = trace[param_name]
                        maxlen = self._getMaxStrLength(param, j)
                        print_line.append(
                            f'{param.value_old:>{maxlen}.{param.decimals}f}')

            if lines_shown == -1:
                print_lines.append(print_line)
            else:
                if lines_shown <= i <= len(self.results) - lines_shown - 1:
                    if i == lines_shown + 1:
                        print_lines.append(f'{"":.>3s}')
                else:
                    print_lines.append(print_line)
        if lines_shown:
            for l in print_lines:
                print(' '.join(l))

    def _getHeader(self, pad=True):
        """
        Construct the table header for the `self.print` and `self.write` methods

        :param pad: Whether to pad the cell values with leading whitespace (for
        printing)
        :type pad: bool

        :return: The header lines for printing or writing to file
        :rtype: list[str]
        """

        fields = []
        if len(self.results[0]) == 2:
            fields_line = []
            for parameter in self.results[0][0]:
                maxlen = self._getMaxStrLength(parameter)
                if pad:
                    fields_line.append(f'{"":>{maxlen}s}')
                else:
                    fields_line.append("")
                if self.results[0][0][parameter].name_alt:
                    if pad:
                        fields_line.append(f'{"":>{maxlen}s}')
                    else:
                        fields_line.append("")
            for i, trace in enumerate(self.results[0][1], start=1):
                for parameter_name in PARAMETERS_TRACE.keys():
                    parameter = trace[parameter_name]
                    maxlen = self._getMaxStrLength(parameter, i)
                    if parameter.name == 'Cm':
                        if pad:
                            fields_line.append(f'{trace["elem"]:>{maxlen}s}')
                        else:
                            fields_line.append(trace["elem"])
                    if parameter.name == 'em':
                        if pad:
                            fields_line.append(
                                f'{trace["isoratio"]:>{maxlen}s}')
                        else:
                            fields_line.append(trace["isoratio"])
            fields.append(fields_line)

        fields_line = []
        for name, parameter in self.results[0][0].items():
            maxlen = self._getMaxStrLength(name)
            if pad:
                fields_line.append(f'{parameter.name:>{maxlen}s}')
            else:
                fields_line.append(parameter.name)
            if parameter.name_alt:
                if pad:
                    fields_line.append(f'{parameter.name_alt:>{maxlen}s}')
                else:
                    fields_line.append(parameter.name_alt)
        if len(self.results[0]) == 2:
            for i, trace in enumerate(self.results[0][1], start=1):
                for parameter_name in PARAMETERS_TRACE.keys():
                    parameter = trace[parameter_name]
                    maxlen = self._getMaxStrLength(parameter, i)
                    if parameter.name in PARAMETERS_TRACE.keys():
                        if pad:
                            fields_line.append(f'{parameter_name:>{maxlen}s}')
                        else:
                            fields_line.append(parameter_name)
            fields.append(fields_line)
        return fields

    def _getMaxStrLength(self, parameter, trace_iter=0):
        """
        Get the maximum string length for this particular column.

        :param parameter: The parameter of interest
        :type parameter: `ecafc.Parameter`

        :return: The maximum string length
        :rtype: int
        """

        maxlen = MAXLEN
        if trace_iter > 0:
            trace = self.results[0][1][trace_iter - 1]
            if len(trace[PARAMETERS_TRACE[parameter.name]]) > maxlen:
                maxlen = len(trace[PARAMETERS_TRACE[parameter.name]])
            if len(trace[PARAMETERS_TRACE[parameter.name]]) > maxlen:
                maxlen = len(trace[PARAMETERS_TRACE[parameter.name]])
        else:
            param = self.results[0][0][parameter]
            if len(param.name) > maxlen:
                maxlen = len(param.name)
            if param.name_alt:
                if len(param.name_alt) > maxlen:
                    maxlen = len(param.name_alt)
            if len(param.str) > maxlen:
                maxlen = len(param.str)
        return maxlen

    def write(self, outname):
        """
        Write the results to an output .csv file.

        :param outname: Name of the output .csv file
        :type outname: str
        """

        write_lines = self._getHeader(pad=False)
        for line in write_lines:
            for cell in line:
                cell = cell.strip()
        for i, row in enumerate(self.results):
            write_line = []
            for name, param in row[0].items():
                maxlen = self._getMaxStrLength(name)
                write_line.append(param.value_old)
                if param.name_alt:
                    value_alt = equil.unnormalize_temp(param.value_old)
                    write_line.append(value_alt)
            if len(row) == 2:
                for j, trace in enumerate(row[1], start=1):
                    for param_name in PARAMETERS_TRACE:
                        param = trace[param_name]
                        maxlen = self._getMaxStrLength(param, j)
                        write_line.append(param.value_old)
            write_lines.append(write_line)

        with open(outname, 'w', newline='') as file_handle:
            writer = csv.writer(file_handle, quotechar='"')
            for l in write_lines:
                writer.writerow(l)