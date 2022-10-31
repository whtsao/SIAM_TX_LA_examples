#!/bin/bash
#SBATCH -N 4
#SBATCH -n 192
#SBATCH -c 1 # specify 6 threads per process
#SBATCH -t 24:00:00
#SBATCH -p workq
#SBATCH -A loni_proteus01s
#SBATCH -o o.out # optional, name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e e.err # optional, name of the stderr, using job and first node values

date

module purge
module load proteus/1.8.1
#which python
#unset PYTHONPATH
#which python
#export PROTEUS_ARCH=linux
#export PYTHONPATH=/project/cekees/cekees/proteus/1.7.5/lib/python3.7/site-packages:$PYTHONPATH
#export PATH=${PROJECT}/bin:${PROJECT}/proteus/1.7.5/bin:${PATH}
#export LD_LIBRARY_PATH=${PROJECT}/proteus/1.7.5/lib:${PROJECT}/proteus/1.7.5/lib64:${LD_LIBRARY_PATH}

#. ~/.bashrc
mkdir -p $WORK/beji_absorp_flat.$SLURM_JOBID
cd $WORK/beji_absorp_flat.$SLURM_JOBID 
cp $SLURM_SUBMIT_DIR/*.py .
cp $SLURM_SUBMIT_DIR/petsc.options.superlu_dist .
cp $SLURM_SUBMIT_DIR/petsc.options.asm .
cp $SLURM_SUBMIT_DIR/beji_battjes_qbc.sh .

#conda activate proteus-dev
#unset PYTHONPATH

#srun parun beji_battjes_so.py -F -l 5 -C "he=0.001 T=300.0" -O petsc.options.superlu_dist
srun parun beji_battjes_so.py -F -l 5 -C "he=0.008 T=200. wave_period=2.01975" -O petsc.options.asm

date

exit 0

