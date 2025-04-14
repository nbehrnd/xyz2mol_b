#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-2.0-only

# name:    test_xyz2mol_b.py
# author:  nbehrnd@yahoo.com
# license: GPL v2, 2025
# date:    [2025-04-14 Mon]
# edit:

"""pytest script of xyz2mol_b.py

Black box-tests to probe the wrapper `xyz2mol_b.py`, i.e. to check from the
outside if the xyz2sdf conversion by Jan Jensen's conversion function still
is available and functional.
"""

import os
import subprocess

PRG = "xyz2mol_b.py"


def test_script_exists():
    """check for the script's presence"""
    assert os.path.isfile(PRG), f"script {PRG} was not found"
