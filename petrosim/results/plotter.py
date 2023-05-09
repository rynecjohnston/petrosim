import argparse
import csv
import re

import matplotlib as mpl
import matplotlib.pyplot as plt


NUM_HEADER_ROWS = 2
REPLACEMENTS = {
    'Tm': r'$T_m$',
    'Tm,norm': r'$T_m,norm$',
    'Ta': r'$T_a$',
    'Ta,norm': r'$T_a,norm$',
    'Mm': r'$M_m$',
    'dm': r'$δ_m$',
    'Cm': r'$C_m$',
    'em': r'$ε_m$',
    'Mr': r'$M_r$',
}
UNITS = {
    'Tm': 'K',
    'Ta': 'K',
    'Cm': 'ppm',
}


def parse_cmdline():
    p = argparse.ArgumentParser(
        description='Make a 2D plot comparing two quantities.')
    p.add_argument('incsv', help='Filename of the output CSV file to analyze.')
    p.add_argument(
        '-x',
        required=True,
        help=
        'Column name of the x-axis data from the output CSV file. (If column name spans two rows, then combine them with a space, e.g., "87Sr/86Sr em"'
    )
    p.add_argument(
        '-y',
        required=True,
        help=
        'Column name of the y-axis data from the output CSV file. (If column name spans two rows, then combine them with a space, e.g., "87Sr/86Sr em"'
    )
    p.add_argument('-title', help='Title of the plot.')
    p.add_argument('-outname', help='Filename for the output PNG file.')
    p.add_argument('-nosave',
                   dest='save',
                   action='store_false',
                   help='Do not save the plot to file.')
    p.add_argument('-noshow',
                   dest='show',
                   action='store_false',
                   help='Do not show the plot on screen.')
    p.add_argument('-rev_x', action='store_true', help='Reverse x-axis.')
    p.add_argument('-rev_y', action='store_true', help='Reverse y-axis.')
    return p.parse_args()


def combine_header_rows(reader, num_lines_header=2):
    """
    If trace elements are requested in the simulation, then the CSV file header
    will span two lines. In order to register the requested column with the
    correct data, zip the two lines to make the new header fields.

    :param reader: The CSV reader
    :type reader: `_csv.reader`

    :param num_lines_header: Number of lines to combine in the header
    :type num_lines_header: int

    :return: The combined header fields
    :rtype: list[str]
    """

    header = []
    new_header = []
    for i, row in enumerate(reader):
        if i < num_lines_header:
            header.append(row)
    for pair in zip(*header):
        if pair[0]:
            new_header.append(' '.join(pair))
        else:
            new_header.append(pair[1])
    return new_header


def csv_reader(filename):
    """
    :param filename: The filename of the input CSV file
    :type filename: str

    :yield: The lines of the CSV parsed into a dict by field keys
    :rtype: generator[dict[str: float]]
    """

    with open(filename, 'r', newline='') as file_handle:
        reader = csv.reader(file_handle)
        new_header = combine_header_rows(reader, NUM_HEADER_ROWS)

        file_handle.seek(0)

        for i, row in enumerate(reader):
            if i < NUM_HEADER_ROWS:
                continue
            row_dict = dict(zip(new_header, [float(r) for r in row]))
            yield row_dict


def get_data(reader, col):
    """
    Collect the data of the specified column name from the reader

    :param reader: The CSV lines
    :type reader: list[dict[str: float]]

    :param col: Name of the column of interest
    :type col: str

    :return: The data in the column col
    :rtype: list[float]
    """

    data = []
    for row in reader:
        value = row.get(col)
        if value is not None:
            data.append(value)
        else:
            exit(f'Could not find "{col}" column in output file.')
    return data


def plot_cols(reader,
              x_col,
              y_col,
              title=None,
              filename=None,
              save=True,
              show=True,
              rev=(False, False)):
    """
    Generate the 2D plot of the specified x and y values.

    :param reader: The CSV lines
    :type reader: list[dict[str: float]]

    :param x_col: Name of the column for the x-axis
    :type x_col: str

    :param y_col: Name of the column for the y-axis
    :type y_col: str

    :param title: Title of the plot
    :type title: str

    :param filename: Name to give the output PNG file.
    :type filename: str

    :param save: Whether to save the plot to file.
    :type save: bool

    :param show: Whether to show the plot on screen.
    :type show: bool

    :param rev: Which axis to reverse.
    :type rev: Tuple[bool]
    """

    x = get_data(reader, x_col)
    y = get_data(reader, y_col)

    fig, ax = plt.subplots(figsize=(6.667, 6))
    ax.plot(x, y)
    if title:
        ax.set_title(title, fontsize=16)
    if rev[0]:
        ax.invert_xaxis()
    if rev[1]:
        ax.invert_yaxis()
    ax.set_xlabel(pretty_label(x_col), fontsize=14)
    ax.set_ylabel(pretty_label(y_col), fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=14)
    if not filename:
        filename = construct_filename(x_col, y_col)
    plt.tight_layout()
    if save:
        plt.savefig(filename, bbox_inches='tight')
    if show:
        plt.show()


def pretty_label(string):
    """
    Format the string to look nice, e.g., use subscripts, correct unicode chars,
    and detect isotopes.

    :param string: The string to format
    :type string: str

    :return: The formatted string
    :rtype: str
    """

    new_strings = []
    if len(string.split()) == 2:
        specific, param = string.split()
        if specific.find('/') != -1:
            isotopes = []
            for elem in specific.split('/'):
                isotopes.append(isotopize(elem))
            specific = '/'.join(isotopes)
            new_strings.append(specific)
        else:
            isotope = isotopize(specific)
            new_strings.append(isotope)
    else:
        param = string
    param_fmt = REPLACEMENTS[param] if param in REPLACEMENTS else param
    new_strings.append(param_fmt)
    if param in UNITS:
        new_strings.append(f'({UNITS[param]})')
    return ' '.join(new_strings)


def isotopize(string):
    """
    Correct string to display isotopes properly.

    :param string: The string to format
    :type string: str

    :return: The formatted string
    :rtype: str
    """

    isotope = ''
    atom_wt = ''
    elem = ''
    for s in string:
        if s.isdigit():
            isotope += rf'$^{s}$'
        else:
            isotope += s
    return isotope


def construct_filename(x_col, y_col):
    """
    Construct the filename for the 2D plot from the x- and y-column names.

    :param x_col: Name of the column for the x-axis
    :type x_col: str

    :param y_col: Name of the column for the y-axis
    :type y_col: str

    :return: Name of the output file
    :rtype: str
    """

    cols_fmt = []
    for col in (x_col, y_col):
        col = re.sub('/', '-', col)
        col = re.sub(' ', '_', col)
        cols_fmt.append(col)
    return '-'.join(cols_fmt)


def main():
    args = parse_cmdline()

    reader = list(csv_reader(args.incsv))
    plot_cols(reader,
              args.x,
              args.y,
              title=args.title,
              filename=args.outname,
              save=args.save,
              show=args.show,
              rev=(args.rev_x, args.rev_y))


if __name__ == '__main__':
    main()