#!/bin/bash
export IFS=","

cat sim-args.csv | while read a b; do 

job_file="covidsim_${a}.job"



echo "#!/bin/bash

#SBATCH  -p bigmem
#SBATCH --nodes=1             
#SBATCH --ntasks-per-node=1      
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --output covidsim_${a}.log
#SBATCH --mem=192GB

ml R/4.0.2

Rscript /home/users/mharris9/COVID19_early_model-main/sims-model_assess.R "$a" " > $job_file

    sbatch $job_file

done



