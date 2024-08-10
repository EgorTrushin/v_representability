#!/usr/bin/env python3
"""Determines progress for a bunch of Molplro calculations running in
given directory."""

import os

calc_dirs = []

nstarted = 0
nfinished = 0
for root, _, files in os.walk("."):
    for name in files:
        if name == "input":
            calc_dirs.append(root)
        if name == "output":
            nstarted += 1
            lfinished = False
            lconverged = True
            with open(os.path.join(root, "output"), encoding="utf-8") as file_obj:
                outtext = file_obj.read()
                if "Molpro calculation terminated" in outtext:
                    nfinished += 1
                    lfinished = True
                if "SCF NOT converged" in outtext:
                    lconverged = False
            if lfinished is False:
                print("UNFINISHED:", root)
            if lconverged is False:
                print("NOT CONVERGED:", root)


print("Number of calculations:", len(calc_dirs))
print("Number of started calculations:", nstarted)
print("Number of finished calculations:", nfinished)
