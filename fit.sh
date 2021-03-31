#!/bin/bash
export IFS=","

cat fit_dates.csv | while read a b c d; do 

job_file="fit_${a}_${c}.job"



echo "#!/bin/bash

#SBATCH  -p normal
#SBATCH --nodes=1              
#SBATCH --ntasks-per-node=${c}      
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --output fit_${a}_${c}.log


ml R/4.0.2

Rscript ./COVID_fit.R "$a" " > $job_file

    sbatch $job_file

done



