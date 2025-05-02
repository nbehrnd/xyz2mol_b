#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-2.0-only

# name:    test_xyz2mol_b.py
# author:  nbehrnd@yahoo.com
# license: GPL v2, 2025
# date:    [2025-04-14 Mon]
# edit:    [2025-05-02 Fri]

"""pytest script of xyz2mol_b.py

Black box-tests to probe the wrapper `xyz2mol_b.py`, i.e. to check from the
outside if the xyz2sdf conversion by Jan Jensen's conversion function still
is available and functional.
"""

import os
import subprocess

import pytest

PRG = "src/xyz2mol_b/xyz2mol_b.py"


def test_script_exists():
    """check for the script's presence"""
    assert os.path.isfile(PRG), f"script {PRG} was not found"


# section of black box-tests
#
# Jan Jensen's repository <https://github.com/jensengroup/xyz2mol> includes
# many tests about the inner working of `xyz2mol`.  This is why, at present,
# I constrain the scope of testing the moderator script to black box-tests
# without proper import of of xyz2mol's functions.
#
# The checks probe the recovery of the neutral ethane molecule, acetate anion,
# and tetramethylammonium cation (chiefly `NMe4_cation`) the .xyz files store
# without bond information as either neutral molecule, anion, or cation.


@pytest.mark.ethane
def test_neutral_ethane():
    """probe recovery of ethane as an implicitly neutral molecule"""
    input_file = os.path.join("tests", "ethane.xyz")
    expected_output = """

     RDKit          3D

  8  7  0  0  0  0  0  0  0  0999 V2000
   -4.5873    0.9270    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1105    0.9270    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9379    1.7888    0.5806 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9379   -0.0068    0.4561 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9379    0.9989   -1.0367 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7600    0.8550    1.0367 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7600    1.8607   -0.4561 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7600    0.0651   -0.5806 H   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0
  3  1  1  0
  4  1  1  0
  5  1  1  0
  6  2  1  0
  7  2  1  0
  8  2  1  0
M  END

    """
    expected_lines = expected_output.strip().splitlines()

    result = subprocess.run(
        f"python {PRG} {input_file}",
        shell=True,
        check=True,
        capture_output=True,
        text=True,
    )
    reported_lines = result.stdout.strip().splitlines()

    assert (
        "The total charge -c is incompatible with bond orders assigned."
        not in result.stdout
    ), "error note issued by the script"
    assert reported_lines == expected_lines, "mismatch sdf output"


@pytest.mark.ethane
def test_neutral_ethane_b():
    """probe recovery of ethane as an explicitly neutral molecule"""
    input_file = os.path.join("tests", "ethane.xyz")
    result = subprocess.run(
        f"python {PRG} {input_file} -c 0",
        shell=True,
        check=True,
        capture_output=True,
        text=True,
    )
    assert (
        "The total charge -c is incompatible with bond orders assigned."
        not in result.stdout
    )


@pytest.mark.ethane
def test_neutral_ethane_c():
    """probe recovery of ethane as neutral molecule"""
    input_file = os.path.join("tests", "ethane.xyz")
    result = subprocess.run(
        f"python {PRG} {input_file} --charge 0",
        shell=True,
        check=True,
        capture_output=True,
        text=True,
    )
    assert (
        "The total charge -c is incompatible with bond orders assigned."
        not in result.stdout
    )


@pytest.mark.ethane
def test_negative_ethane():
    """probe recovery of ethane as an anion"""
    input_file = os.path.join("tests", "ethane.xyz")
    result = subprocess.run(
        f"python {PRG} {input_file} -c -1",
        shell=True,
        check=True,
        capture_output=True,
        text=True,
    )
    assert (
        "The total charge -c is incompatible with bond orders assigned."
        in result.stdout
    )


@pytest.mark.acetate
def test_neutral_acetate():
    """probe recovery of acetate as a neutral molecule"""
    input_file = os.path.join("tests", "acetate.xyz")
    result = subprocess.run(
        f"python {PRG} {input_file}",
        shell=True,
        check=True,
        capture_output=True,
        text=True,
    )
    assert (
        "The total charge -c is incompatible with bond orders assigned."
        in result.stdout
    )


@pytest.mark.acetate
def test_negative_acetate():
    """probe recovery of acetate as an anion"""
    input_file = os.path.join("tests", "acetate.xyz")
    expected_output = """

     RDKit          3D

  7  6  0  0  0  0  0  0  0  0999 V2000
   -4.7169    0.8992    0.0571 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2490    0.9840   -0.2283 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0417    1.7438    0.6786 H   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0171   -0.0221    0.5634 H   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2108    0.9687   -0.9121 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6591    2.0570   -0.3402 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6341   -0.1870   -0.4868 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0
  3  1  1  0
  4  1  1  0
  5  1  1  0
  6  2  2  0
  7  2  1  0
M  CHG  1   7  -1
M  END

    """
    expected_lines = expected_output.strip().splitlines()

    result = subprocess.run(
        f"python {PRG} {input_file} -c -1",
        shell=True,
        check=True,
        capture_output=True,
        text=True,
    )
    reported_lines = result.stdout.strip().splitlines()

    assert (
        "The total charge -c is incompatible with bond orders assigned."
        not in result.stdout
    ), "error note issued by the script"
    assert reported_lines == expected_lines, "mismatch sdf output"


@pytest.mark.acetate
def test_negative_acetate_sdfV3000():
    """probe recovery of acetate as an anion"""
    input_file = os.path.join("tests", "acetate.xyz")
    expected_output = """

     RDKit          3D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -4.716860 0.899190 0.057140 0
M  V30 2 C -3.248980 0.984000 -0.228300 0
M  V30 3 H -5.041670 1.743840 0.678620 0
M  V30 4 H -5.017100 -0.022050 0.563440 0
M  V30 5 H -5.210760 0.968740 -0.912080 0
M  V30 6 O -2.659090 2.057020 -0.340250 0
M  V30 7 O -2.634130 -0.187020 -0.486790 0 CHG=-1
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 3 1
M  V30 3 1 4 1
M  V30 4 1 5 1
M  V30 5 2 6 2
M  V30 6 1 7 2
M  V30 END BOND
M  V30 END CTAB
M  END

    """
    expected_lines = expected_output.strip().splitlines()

    result = subprocess.run(
        f"python {PRG} {input_file} -c -1 --x3",
        shell=True,
        check=True,
        capture_output=True,
        text=True,
    )
    reported_lines = result.stdout.strip().splitlines()

    assert (
        "The total charge -c is incompatible with bond orders assigned."
        not in result.stdout
    ), "error note issued by the script"
    assert reported_lines == expected_lines, "mismatch sdf output"


@pytest.mark.acetate
def test_positive_acetate():
    """probe recovery of acetate as a cation"""
    input_file = os.path.join("tests", "acetate.xyz")
    result = subprocess.run(
        f"python {PRG} {input_file} -c +1",
        shell=True,
        check=True,
        capture_output=True,
        text=True,
    )
    assert (
        "The total charge -c is incompatible with bond orders assigned."
        in result.stdout
    )


@pytest.mark.NMe4
def test_neutral_NMe4_cation():
    """probe recovery of NMe4_cation as a neutral molecule"""
    input_file = os.path.join("tests", "NMe4_cation.xyz")
    result = subprocess.run(
        f"python {PRG} {input_file}",
        shell=True,
        check=True,
        capture_output=True,
        text=True,
    )
    assert (
        "The total charge -c is incompatible with bond orders assigned."
        in result.stdout
    )


@pytest.mark.NMe4
def test_negative_NMe4_cation():
    """probe recovery of NMe4_cation as an anion"""
    input_file = os.path.join("tests", "NMe4_cation.xyz")
    result = subprocess.run(
        f"python {PRG} {input_file} -c -1",
        shell=True,
        check=True,
        capture_output=True,
        text=True,
    )
    assert (
        "The total charge -c is incompatible with bond orders assigned."
        in result.stdout
    )


@pytest.mark.NMe4
def test_positive_NMe4_cation():
    """probe recovery of NMe4_cation as a cation"""
    input_file = os.path.join("tests", "NMe4_cation.xyz")
    expected_output = """

     RDKit          3D

 17 16  0  0  0  0  0  0  0  0999 V2000
    1.0681    0.0726    0.0097 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.5628   -1.3045    0.3921 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5628    1.0923    1.0111 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5628    0.4299   -1.3741 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5840    0.0726    0.0097 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9374   -1.5441    1.3922 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5315   -1.2844    0.3866 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.9374   -2.0255   -0.3412 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.9374    2.0782    0.7186 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5315    1.0775    0.9965 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.9374    0.8177    2.0021 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.9374    1.4255   -1.6319 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.9374   -0.3164   -2.0817 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5315    0.4247   -1.3539 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.9291   -0.1871    1.0154 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.9291   -0.6685   -0.7180 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.9291    1.0734   -0.2682 H   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0
  3  1  1  0
  4  1  1  0
  5  1  1  0
  6  2  1  0
  7  2  1  0
  8  2  1  0
  9  3  1  0
 10  3  1  0
 11  3  1  0
 12  4  1  0
 13  4  1  0
 14  4  1  0
 15  5  1  0
 16  5  1  0
 17  5  1  0
M  CHG  1   1   1
M  END

    """
    expected_lines = expected_output.strip().splitlines()

    result = subprocess.run(
        f"python {PRG} {input_file} -c 1",
        shell=True,
        check=True,
        capture_output=True,
        text=True,
    )
    reported_lines = result.stdout.strip().splitlines()

    assert (
        "The total charge -c is incompatible with bond orders assigned."
        not in result.stdout
    )
    assert reported_lines == expected_lines, "mismatch sdf output"
