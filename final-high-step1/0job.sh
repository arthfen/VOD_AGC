#PBS -N highstep1
#PBS -q longp
#PBS -l nodes=1:ppn=13
#PBS -e /dev/null
#PBS -o /dev/null

# Loads all modules
module purge
module load gcc nco
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/users/afendric/.local/lib

# Runs all files
cd $PBS_O_WORKDIR

it=1
cd ../final-high-step1/ &&
/home/users/afendric/.local/bin/R --no-save < 0call.R > output_ML1 2>&1 &&
mv output_ML1 output_ML1-step"$it" &&
cd ../final-high-step4/ &&
/home/users/afendric/.local/bin/R --no-save < 0call.R > output_ML1 2>&1 &&
mv output_ML1 output_ML1-step"$it"

