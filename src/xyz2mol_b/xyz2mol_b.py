#!/usr/bin/env python3
"""
name   : xyz2mol_b.py
author : nbehrnd@yahoo.com
license: GPL2
date   : [2023-08-08 Tue]
edit   : [2025-04-15 Tue]
purpose: a modernized xyz2mol interface to RDKit 2023.03.2, or greater
"""

import argparse

import rdkit
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds


def get_args():
    """get arguments from the CLI"""

    parser = argparse.ArgumentParser(
        description="""
Convert a .xyz file into a .sdf file.  This is a moderator to Jan Jensen's
xyz2mol, a functionality which eventually was included into RDKit.  Processing
requires RDKit version 2023.03.2, or greater.  By default, small molecules with
less than 999 atoms are reported in the syntax of sdf (V2000); the optional
flag `--x3` allows to use the syntax of sdf (V3000) instead.""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "file",
        help=".xyz input file to process",
        metavar="FILE",
        type=argparse.FileType("rt"),
        default=None,
    )

    parser.add_argument(
        "-c",
        "--charge",
        help="user assigned total charge of the input structure",
        metavar="int",
        type=int,
        default=0,
    )

    parser.add_argument(
        "--x3",
        help="report structure in syntax of .sdf (V3000)",
        action="store_true",
    )

    return parser.parse_args()


def main():
    """join the functionalities"""

    args = get_args()

    raw_mol = Chem.MolFromXYZFile(args.file.name)
    mol = Chem.Mol(raw_mol)

    try:
        rdkit.Chem.rdDetermineBonds.DetermineBonds(mol, charge=args.charge)
        if args.x3:
            print(Chem.MolToMolBlock(mol, forceV3000=True))
        else:
            print(Chem.MolToMolBlock(mol))
    except ValueError:
        print("The total charge -c is incompatible with bond orders assigned.")


if __name__ == "__main__":
    main()
