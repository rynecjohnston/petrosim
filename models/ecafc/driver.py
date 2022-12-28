"""
This is the driver for the Energy-Constrained Assimilation Fractional
Crystallization (EC-AFC) model based on work by Spera & Bohrson J. Petrol., 42,
999–1018, 2001.
"""

import argparse
import csv
import dataclasses
import os

import ecafc
import equilibration as equil
import loader as ldr


def parse_cmdline():
    p = argparse.ArgumentParser(description='EC-AFC: Energy-Constrained Assimilation Fractional Crystallization petrological model\nbased on Spera & Bohrson J. Petrol., 42, 999–1018, 2001.')
    p.add_argument('infile', help='Input YAML file.')
    p.add_argument('outcsv', help='Output CSV file.')
    p.add_argument('-print', type=int, default=0, help='Number of lines of results to print. (n >= 1 to print n lines, n = 0 for no print, and n = -1 to print all lines)')
    p.add_argument('-Teq', type=float, default=None, help='Equilibrate to this alternative temperature (ºC).')
    p.add_argument('-Teq_norm', type=float, default=None, help='Equilibrate to this alternative normalized temperature.')
    return p.parse_args()


def validate_args(args):
    if not os.path.exists(args.infile):
        exit(f'Could not find input file {args.infile}')
    if args.Teq_norm == 0.0 or (args.Teq_norm and args.Teq_norm <= 0.0):
        exit(f'Teq_norm must be a positive value.')
    if args.Teq and args.Teq <= -273.16:
        exit(f'Teq must be above -273.16 ºC.')


def read_input_file(infile):
    ldr.init(infile)
    return ldr.Parameters(**dataclasses.asdict(ldr.params))


def write_output_csv(outcsv, results):
    with open(outcsv, 'w', newline='') as filehandle:
        writer = csv.writer(filehandle)
        writer.writerows(results)


def main():
    args = parse_cmdline()
    validate_args(args)

    params = read_input_file(args.infile)

    if args.Teq_norm is not None:
        params.Teq_norm = args.Teq_norm
    if args.Teq is not None:
        params.Teq = args.Teq
        params.Teq_norm = equil.normalize_temp(params.Teq)

    afc = ecafc.ECAFC(params)
    afc.simulate(print_lines=args.print)
    write_output_csv(args.outcsv, afc.results)


if __name__ == '__main__':
    main()