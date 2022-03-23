import pyximport; pyximport.install()
from powtnrka import make_parser, process_args, run_powtnrka

LICENSE = """
PowTNRka   Copyright (C) 2024   Alvin Hsu
This program comes with ABSOLUTELY NO WARRANTY!
"""

def main():
    print(LICENSE)
    parser = make_parser()
    args = parser.parse_args()
    run_powtnrka(**process_args(args))


if __name__ == '__main__':
    main()
