#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-2.0-only

# name:    test_xyz2mol_b.py
# author:  nbehrnd@yahoo.com
# license: GPL v2, 2025
# date:    [2025-04-14 Mon]
# edit:    [2025-04-15 Tue]

"""pytest script of xyz2mol_b.py

Black box-tests to probe the wrapper `xyz2mol_b.py`, i.e. to check from the
outside if the xyz2sdf conversion by Jan Jensen's conversion function still
is available and functional.
"""

import os
import subprocess

import pytest

PRG = "xyz2mol_b/__init__.py"


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
    """probe recovery of ethane as neutral molecule"""
    input_file = os.path.join("tests", "ethane.xyz")
    result = subprocess.run(
        f"python {PRG} {input_file}",
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
def test_neutral_ethane_b():
    """probe recovery of ethane as neutral molecule"""
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
    result = subprocess.run(
        f"python {PRG} {input_file} -c -1",
        shell=True,
        check=True,
        capture_output=True,
        text=True,
    )
    assert (
        "The total charge -c is incompatible with bond orders assigned."
        not in result.stdout
    )


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
    result = subprocess.run(
        f"python {PRG} {input_file} -c 1",
        shell=True,
        check=True,
        capture_output=True,
        text=True,
    )
    assert (
        "The total charge -c is incompatible with bond orders assigned."
        not in result.stdout
    )
