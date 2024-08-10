#!/usr/bin/env python3
"""Clean calculation directories."""

import os

calc_dirs = []

for root, _, files in os.walk("."):
    for name in files:
        if name == "input":
            calc_dirs.append(root)

for calc_dir in calc_dirs:
    os.system(f"rm {calc_dir}/run.sh")
    os.system(f"rm {calc_dir}/slurm*")
    os.system(f"rm {calc_dir}/error")
