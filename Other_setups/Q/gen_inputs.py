#!/usr/bin/env python3
"""Generation of input files for calculations."""

# pylint: disable=E1101

import argparse
import json
import os
from pprint import pprint
import yaml
from utils.slurm_runsh import create_runsh


def get_args():
    """Gets the command-line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "config_path",
        nargs="?",
        type=str,
        default="configs/config.yaml",
        help="Path to config file",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="increase output verbosity")
    return parser.parse_args()


def get_spin(atom):
    """Determines spin for given atom"""
    if atom in ["C", "O", "Si", "S"]:
        spin = 2
    elif atom in ["B", "F", "Al", "Cl"]:
        spin = 1
    elif atom in ["N", "P"]:
        spin = 3
    return spin


def input_exx(system, method_dict):
    """Generate input text for USCEXX calculations."""

    # Determine OEP basis
    if system in ["B", "C", "N", "O", "F"]:
        oep_basis = "default,aug-cc-pVDZ/mp2fit"
    elif system in ["Al", "Si", "P", "S", "Cl"]:
        oep_basis = "default,aug-cc-pVTZ/mp2fit"

    # Determine spin
    spin = get_spin(system)

    # Represent method dict as string
    method = ""
    for key in method_dict.keys():
        method += f",{key}={method_dict[key]}"

    return f"""basis={{
default,aug-cc-pwCVQZ
set,oep
{oep_basis}
}}

symmetry,nosym

angstrom
geometry={{
1

{system} 0.0 0.0 0.0
}}

spin={spin}

uhf,maxit=0

acfd;uscexx,dfit=0,maxit=100,energy=1d-10,vhoep=1,thr_fai_oep=1.7d-2{method}
"""


def gen_inputs_exx(cfg):
    """Generates batch of inputs for EXX calculations."""
    for system in cfg.systems:
        for symmetry in cfg.symmetry:
            if len(cfg.symmetry) > 1:
                subdir = os.path.join(cfg.output_dir, system, symmetry)
            else:
                subdir = os.path.join(cfg.output_dir, system)
            os.makedirs(subdir)

            # Prepare dictionary with parameters for calculation
            method_dict = {}
            if symmetry in ["SpinSym", "FullSym"]:
                method_dict["spin_sym"] = 1
            if symmetry in ["SpaceSym", "FullSym"]:
                method_dict["thr_sym"] = "1d-10"
            if symmetry == "SpinSym":
                if system in ["O", "F"]:
                    method_dict["swap"] = 1
                    method_dict["iswap1"] = 3
                    method_dict["iswap2"] = 5
                elif system in ["S", "Cl"]:
                    method_dict["swap"] = 1
                    method_dict["iswap1"] = 7
                    method_dict["iswap2"] = 9
            if cfg.plot:
                method_dict["plot_x"] = 1
                method_dict["plot_y"] = 1
                method_dict["plot_z"] = 1

            # Create input file
            molpro_input = os.path.join(subdir, "input")
            with open(molpro_input, "w", encoding="utf8") as infile_obj:
                print(input_exx(system, method_dict), file=infile_obj)

            # Create run.sh script, if required
            if cfg.create_runsh:
                create_runsh(cfg.molpro_path, cfg.partition, subdir)


def input_cahf(system):
    """Generate input text for Molpro."""

    if system == "B":
        method = "{cahf;closed,1;shell,3,2.1,3.1,4.1,5.1}"
    elif system == "C":
        method = "{cahf;closed,2;shell,2,3.1,4.1,5.1}"
    elif system == "N":
        method = "{cahf;closed,2;shell,3,3.1,4.1,5.1}"
    elif system == "O":
        method = "{cahf;closed,2;shell,4,3.1,4.1,5.1}"
    elif system == "F":
        method = "{cahf;closed,2;shell,5,3.1,4.1,5.1}"
    elif system == "Al":
        method = "{cahf;closed,5;shell,3,6.1,7.1,8.1,9.1}"
    elif system == "Si":
        method = "{cahf;closed,6;shell,2,7.1,8.1,9.1}"
    elif system == "P":
        method = "{cahf;closed,6;shell,3,7.1,8.1,9.1}"
    elif system == "S":
        method = "{cahf;closed,6;shell,4,7.1,8.1,9.1}"
    elif system == "Cl":
        method = "{cahf;closed,6;shell,5,7.1,8.1,9.1}"

    return f"""basis={{
default,aug-cc-pwCVQZ
}}

symmetry,nosym

angstrom
geometry={{
1

{system} 0.0 0.0 0.0
}}

{method}
"""


def gen_inputs_cahf(cfg):
    """Generates bunch of CAHF inputs."""
    os.mkdir(cfg.output_dir)

    for system in cfg.systems:
        subdir = os.path.join(cfg.output_dir, system)
        os.makedirs(subdir)

        molpro_input = os.path.join(subdir, "input")
        with open(molpro_input, "w", encoding="utf8") as file_obj:
            print(input_cahf(system), file=file_obj)


def input_rohf(system):
    """Generate input text for Molpro."""

    # Determine spin
    spin = get_spin(system)

    return f"""basis={{
default,aug-cc-pwCVQZ
}}

symmetry,nosym

angstrom
geometry={{
1

{system} 0.0 0.0 0.0
}}

spin={spin}

hf
"""


def gen_inputs_rohf(cfg):
    """Generates bunch of ROHF inputs."""
    os.mkdir(cfg.output_dir)

    for system in cfg.systems:
        subdir = os.path.join(cfg.output_dir, system)
        os.makedirs(subdir)

        molpro_input = os.path.join(subdir, "input")
        with open(molpro_input, "w", encoding="utf8") as file_obj:
            print(input_rohf(system), file=file_obj)


def input_ksinv_rohf(system, mb_method, method_dict):
    """Generate input text for USCEXX calculations."""

    # Determine OEP basis
    if system in ["B", "C", "N", "O", "F"]:
        oep_basis = "default,aug-cc-pVDZ/mp2fit"
    elif system in ["Al", "Si", "P", "S", "Cl"]:
        oep_basis = "default,aug-cc-pVTZ/mp2fit"

    # Determine spin
    spin = get_spin(system)

    # Represent method dict as string
    method = ""
    for key in method_dict.keys():
        method += f",{key}={method_dict[key]}"

    if method_dict["swap"] == 1:
        if system in ["B", "C", "O", "F"]:
            method += ",iswap1=3,iswap2=5"
        else:
            method += ",iswap1=7,iswap2=9"

    with open("json/atoms.json", "r", encoding="utf-8") as file_obj:
        atoms = json.load(file_obj)
    noa = atoms[system]["na"]
    nob = atoms[system]["nb"]

    if mb_method == "RS2":
        mb_method_ = "{rs2;dm,2140.2;core}"
    elif mb_method == "CISD":
        mb_method_ = "{ci;save,density=2140.2,spinden;core}"
    elif mb_method == "AQCC":
        mb_method_ = "{aqcc;save,density=2140.2,spinden;core}"

    return f"""basis={{
default,aug-cc-pwCVQZ
set,oep
{oep_basis}
}}

symmetry,nosym

angstrom
geometry={{
1

{system} 0.0 0.0 0.0
}}

spin={spin}
hf,maxit=30

{mb_method_}
e_ref=energy

acfd;ksinv,maxit=100,refden=2140.2,e_ref=e_ref,vhoep=1,thr_fai_oep=1.7d-2,density_test=1,noa={noa},nob={nob}{method}
"""


def input_ksinv_cahf(system, mb_method, method_dict):
    """Generate input text for USCEXX calculations."""

    # Determine OEP basis
    if system in ["B", "C", "N", "O", "F"]:
        oep_basis = "default,aug-cc-pVDZ/mp2fit"
    elif system in ["Al", "Si", "P", "S", "Cl"]:
        oep_basis = "default,aug-cc-pVTZ/mp2fit"

    # Determine spin
    spin = get_spin(system)

    # Represent method dict as string
    method = ""
    for key in method_dict.keys():
        method += f",{key}={method_dict[key]}"

    if method_dict["swap"] == 1:
        if system in ["B", "C", "O", "F"]:
            method += ",iswap1=3,iswap2=5"
        else:
            method += ",iswap1=7,iswap2=9"

    with open("json/atoms.json", "r", encoding="utf-8") as file_obj:
        atoms = json.load(file_obj)
    noa = atoms[system]["na"]
    nob = atoms[system]["nb"]

    if mb_method == "RS2":
        mb_method_ = "{rs2;dm,2140.2;core}"
    elif mb_method == "CISD":
        mb_method_ = "{ci;save,density=2140.2,spinden;core}"
    elif mb_method == "AQCC":
        mb_method_ = "{aqcc;save,density=2140.2,spinden;core}"

    if system == "B":
        cahf_method = "{cahf;closed,1;shell,3,2.1,3.1,4.1,5.1}"
    elif system == "C":
        cahf_method = "{cahf;closed,2;shell,2,3.1,4.1,5.1}"
    elif system == "N":
        cahf_method = "{cahf;closed,2;shell,3,3.1,4.1,5.1}"
    elif system == "O":
        cahf_method = "{cahf;closed,2;shell,4,3.1,4.1,5.1}"
    elif system == "F":
        cahf_method = "{cahf;closed,2;shell,5,3.1,4.1,5.1}"
    elif system == "Al":
        cahf_method = "{cahf;closed,5;shell,3,6.1,7.1,8.1,9.1}"
    elif system == "Si":
        cahf_method = "{cahf;closed,6;shell,2,7.1,8.1,9.1}"
    elif system == "P":
        cahf_method = "{cahf;closed,6;shell,3,7.1,8.1,9.1}"
    elif system == "S":
        cahf_method = "{cahf;closed,6;shell,4,7.1,8.1,9.1}"
    elif system == "Cl":
        cahf_method = "{cahf;closed,6;shell,5,7.1,8.1,9.1}"

    return f"""basis={{
default,aug-cc-pwCVQZ
set,oep
{oep_basis}
}}

symmetry,nosym

angstrom
geometry={{
1

{system} 0.0 0.0 0.0
}}

{cahf_method}

spin={spin}

{mb_method_}
e_ref=energy

acfd;ksinv,maxit=100,refden=2140.2,e_ref=e_ref,vhoep=1,thr_fai_oep=1.7d-2,density_test=1,noa={noa},nob={nob}{method}
"""


def gen_inputs_ksinv(cfg):
    """Generates batch of KSINV inputs."""
    for system in cfg.systems:
        for mb_method in cfg.mb_methods:
            subdir = os.path.join(cfg.output_dir, mb_method, system)
            os.makedirs(subdir)

            # Prepare dictionary with parameters for calculation
            method_dict = {}
            if mb_method == "RS2":
                if cfg.reference == "ROHF" and system in ["O", "F", "S", "Cl"]:
                    method_dict["swap"] = 1
                elif cfg.reference == "CAHF" and system in [
                    "B",
                    "C",
                    "O",
                    "F",
                    "Al",
                    "Si",
                    "S",
                    "Cl",
                ]:
                    method_dict["swap"] = 1
                else:
                    method_dict["swap"] = 0
            elif mb_method == "AQCC" or mb_method == "CISD" and system in ["B", "C", "O", "F", "Al", "Si", "S", "Cl"]:
                method_dict["swap"] = 1
            else:
                method_dict["swap"] = 0

            if cfg.plot:
                method_dict["plot_x"] = 1
                method_dict["plot_y"] = 1
                method_dict["plot_z"] = 1
                method_dict["plot_vxc"] = 1
                method_dict["plot_vref"] = 1
                method_dict["plot_rho_ks"] = 1
                method_dict["plot_rho_ref"] = 1
                method_dict["gridsize"] = 12000

            # Create input file
            molpro_input = os.path.join(subdir, "input")
            with open(molpro_input, "w", encoding="utf8") as infile_obj:
                if cfg.reference == "ROHF":
                    print(
                        input_ksinv_rohf(system, mb_method, method_dict),
                        file=infile_obj,
                    )
                elif cfg.reference == "CAHF":
                    print(
                        input_ksinv_cahf(system, mb_method, method_dict),
                        file=infile_obj,
                    )

                # Create run.sh script, if required
                if cfg.create_runsh:
                    create_runsh(cfg.molpro_path, cfg.partition, subdir)


# Parse command-line arguments
args = get_args()

# Read config file
with open(args.config_path, "r", encoding="utf-8") as cfg_file:
    config = yaml.safe_load(cfg_file)
config = argparse.Namespace(**config)

if args.verbose:
    pprint(vars(config))

if config.calc_type == "EXX":
    gen_inputs_exx(config)
elif config.calc_type == "ROHF":
    gen_inputs_rohf(config)
elif config.calc_type == "CAHF":
    gen_inputs_cahf(config)
elif config.calc_type == "KSINV":
    gen_inputs_ksinv(config)
