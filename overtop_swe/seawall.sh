#!/bin/bash
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -c 1 # specify 6 threads per process
#SBATCH -t 12:00:00
#SBATCH -p workq
#SBATCH -A loni_proteus01s
#SBATCH -o o.out # optional, name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e e.err # optional, name of the stderr, using job and first node values
#SBATCH -J seawall_swe
date

module purge
module load proteus/1.8.1

mkdir -p $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cd $WORK/$SLURM_JOB_NAME.$SLURM_JOBID 
cp $SLURM_SUBMIT_DIR/*.py .
cp $SLURM_SUBMIT_DIR/*.sh .


srun parun -l5 --SWEs seawall.py
#srun parun --SWEs seawall.py -F -l 5
#parun --SWEs --path " + self.path + " "
#                  "-l1 -v seawall.py

date
exit 0

