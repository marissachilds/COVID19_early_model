#!/bin/bash
export IFS=","

cat fit_names.csv | while read a b; do 

job_file="filter_traj_${a}.job"



echo "#!/bin/bash

#SBATCH  -p normal
#SBATCH --nodes=1              
#SBATCH --ntasks-per-node=20      
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --output filter_traj_${a}.log


ml R/4.0.2

Rscript ./COVID_filtering_traj.R "$a" " > $job_file

    sbatch $job_file

done



