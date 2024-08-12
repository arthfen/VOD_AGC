#PBS -N 2020
#PBS -q mediump
#PBS -l nodes=1:ppn=1
#PBS -e /dev/null
#PBS -o /dev/null
#PBS -M /dev/null

# Loads all modules
module purge
module load gcc nco
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/users/afendric/.local/lib

# Runs all files
cd $PBS_O_WORKDIR

export OMPI_MCA_rmaps_base_oversubscribe=true
export mpi_warn_on_fork=0
Rscript rds2tif_year.R $PBS_JOBNAME

