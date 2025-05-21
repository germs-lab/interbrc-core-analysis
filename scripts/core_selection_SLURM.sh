#!/bin/bash

#SBATCH --nodes=1   # Number of nodes to use
#SBATCH --ntasks-per-node=32   # Use 32 processor cores per node 
#SBATCH --time=2-0:0:0   # Walltime limit (DD-HH:MM:SS)
#SBATCH --mem=256G   # Maximum memory per node
#SBATCH --job-name="interbrc_core_sel"   # Job name to display in squeue
#SBATCH --mail-user=bolivar@iastate.edu   # Email address
#SBATCH --mail-type=BEGIN   # Send an email when the job starts
#SBATCH --mail-type=END   # Send an email when the job ends
#SBATCH --mail-type=FAIL   # Send an email if the job fails
#SBATCH --output="/work/adina/bolivar/interbrc-core-analysis/hpc/slurm-%j-004_core_selection.out"   # Job standard output file (%j will be replaced by the slurm job id)
#SBATCH --error="/work/adina/bolivar/interbrc-core-analysis/hpc/slurm-%j-004_core_selection.error"   # Job standard error file (%j will be replaced by the slurm job id)


#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK # Set OMP_NUM_THREADS to the number of CPUs per task we asked for.

##Modules/Singularity
module purge
module load micromamba/1.4.2-lcemqbe # Latest in Nova

# Initialize micromamba
eval "$(micromamba shell hook --shell=bash)"
micromamba activate interbrc_env

# Basic session info
echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

#Run R script
R CMD BATCH R/analysis/004_core_selection_HPC.R # Add your script here

module purge

echo End Job
