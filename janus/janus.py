from .driver import run_janus
import argparse
import sys

def main():

    parser = argparse.ArgumentParser(prog='janus')
    parser.add_argument('input_file', metavar='i', type=str, help='input file name')
    parser.add_argument('-o', type=str,help='output file name')
    args = parser.parse_args()
    file_in = args.input_file
    file_out = args.o
    if file_out is None:
        file_out = 'output.dat'

    sys.stdout = open(file_out, 'w')

    print('running janus')
    print('input file is {}'.format(file_in))
    print('output file is {}'.format(file_out))
    run_janus(filename=file_in)

if __name__ == '__main__':
    main()
