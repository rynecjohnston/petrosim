"""
This is the driver for the Energy-Constrained Recharge, Assimilation, and
Fractional Crystallization (EC-RAFC) model based on work by Spera & Bohrson G3,
3, 12, 2002.
"""

import argparse
import csv
import dataclasses
import os

from petrosim.models.ecrafc import ecrafc
from petrosim.models.ecrafc import loader as ldr


def parse_cmdline():
    p = argparse.ArgumentParser(
        description=
        'EC-RAFC: Energy-Constrained Recharge, Assimilation, and Fractional Crystallization petrological model\nbased on Spera & Bohrson G3, 3, 12, 2002.'
    )
    p.add_argument('infile', help='Input YAML file.')
    p.add_argument('outcsv', help='Output CSV file.')
    p.add_argument(
        '-print',
        type=int,
        default=0,
        help=
        'Number of lines of results to print. (n >= 1 to print n lines, n = 0 for no print, and n = -1 to print all lines)'
    )
    p.add_argument('-Teq',
                   type=float,
                   default=None,
                   help='Equilibrate to this alternative temperature (ºC).')
    p.add_argument(
        '-Teq_norm',
        type=float,
        default=None,
        help='Equilibrate to this alternative normalized temperature.')
    return p.parse_args()


def validate_args(args):
    """
    Do some basic argument validation.

    :param args: User arguments
    :type args: `argparse.ArgumentParser`
    """

    if not os.path.exists(args.infile):
        exit(f'Could not find input file {args.infile}')
    if args.Teq_norm == 0.0 or (args.Teq_norm and args.Teq_norm <= 0.0):
        exit(f'Teq_norm must be a positive value.')
    if args.Teq and args.Teq <= -273.16:
        exit(f'Teq must be above -273.16 ºC.')


def read_input_file(infile):
    """
    Read the input file.

    :param infile: The input filename.
    :type infile: str
    """

    ldr.init(infile)
    return ldr.Parameters(**dataclasses.asdict(ldr.params))


def main():
    args = parse_cmdline()
    validate_args(args)

    params = read_input_file(args.infile)

    if args.Teq_norm is not None:
        params.Teq_norm = args.Teq_norm
    if args.Teq is not None:
        params.Teq = args.Teq
        params.Teq_norm = params.Teq / params.Tlm

    rafc = ecrafc.ECRAFC(params)
    rafc.simulate()
    rafc.results.print(args.print)
    rafc.results.write(args.outcsv)


if __name__ == '__main__':
    main()