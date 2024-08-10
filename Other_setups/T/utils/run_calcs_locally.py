#!/usr/bin/env python3
"""Runs calculations locally."""

import os
from tqdm import tqdm

MOLPRO = "/home/trushin/Molpro-dev/molpro-swap/bin/molpro"
NUM_OMP_THREADS = 16
MEM = 2000

calc_dirs = []

for root, _, files in os.walk("."):
    for name in files:
        if name == "input":
            calc_dirs.append(root)

pbar = tqdm(calc_dirs)

os.system(f"export OMP_NUM_THREADS={NUM_OMP_THREADS}")

script_dir = os.getcwd()

for calc in pbar:
    os.chdir(os.path.join(script_dir, calc))
    os.system(
        f"{MOLPRO} -t {NUM_OMP_THREADS} -m {MEM} --no-xml-output < input > output"
    )

os.chdir(script_dir)
