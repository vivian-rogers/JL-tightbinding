#!/bin/bash

#SBATCH -J JL.p.weyl        # Job Name
#SBATCH -o ../outputs/jl-tb.o%j  # output and error file name (%j expands to jobID) Name of the output file (eg. myMPI.oJobID)
#SBATCH -n 128         # number of tasks per node
#SBATCH -N 1         # total number of mpi tasks requested
#SBATCH -p normal              # Queue name
#SBATCH -t 24:00:00       # Run time (hh:mm:ss) - 0.25 hours (max is 48 hours)
#SBATCH --mail-user=8326909459@tmomail.net    # Email notification address MAKE THIS YOURS, PLEASE
#SBATCH --mail-type=ALL                  # Email at Begin/End of job  (UNCOMMENT)
#SBATCH -A OTH21017


n=128 # number of tasks per node
N=1 # number of tasks
ntot=$((n * N))

export LD_LIBRARY_PATH="" 
export LD_PRELOAD=""
export JULIA_NUM_THREADS=ntot

vncserver
export DISPLAY=:1
julia runs.jl

#modul load gcc
#module load cuda/9.0


### staging files for copy to scratch
#FOLDERPATH=$SCRATCH/mumax3/$FOLDERNAME 
#mkdir $FOLDERPATH
#cp -r * $FOLDERPATH
#cd $FOLDERPATH


### running desired calculations, postprocessing
#ibrun mumax3 $FILENAME


### copy files back to output folder
#cd ../
#cp -r $FOLDERPATH $mumaxOutputPath
