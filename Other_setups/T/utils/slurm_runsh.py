#!/usr/bin/env python3
"""Generation of run.sh script for calculations of cluster."""
import os


def slurm_text(path="/home/trushin/Molpro/molpro/bin/", partition="mem256"):
    """Generate input text for Slurm script."""
    if partition == "mem256":
        cpu_per_task = 16
        mem = 20000
    elif partition == "mem128_s1":
        cpu_per_task = 16
        mem = 10000
    return f"""#!/bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task={cpu_per_task}
#SBATCH --job-name=testjob
#SBATCH --partition={partition}
#SBATCH --mail-user=egor.trushin@fau.de

export I_MPI_DEBUG=5
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so.0

MOLPRO={path}

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_STACKSIZE=3000M
ulimit -s unlimited
export MV2_CPU_MAPPING=0-$(( ${{num}} - 1))

cd $SLURM_SUBMIT_DIR
echo $SLURM_JOB_NODELIST > $SLURM_SUBMIT_DIR/machines_jobid

TMPDIR=/scratch/trushin/${{SLURM_JOBID}}

$MOLPRO/molpro -t {cpu_per_task} -m {mem} --no-xml-output --no-helper-server -d $TMPDIR < input 1> $SLURM_SUBMIT_DIR/output 2> $SLURM_SUBMIT_DIR/error

rm *.wfu *vxdiff machines_jobid
scp $TMPDIR/* tcsv020:$SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR
"""


def create_runsh(molpro_path, partition, subdir):
    """Creates run.sh scripit."""
    slurm_input = os.path.join(subdir, "run.sh")
    with open(slurm_input, "w", encoding="utf8") as runsh_file:
        print(slurm_text(molpro_path, partition), file=runsh_file)
