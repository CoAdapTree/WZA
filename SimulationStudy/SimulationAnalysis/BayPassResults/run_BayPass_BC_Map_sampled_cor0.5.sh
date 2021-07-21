#!/bin/bash
#SBATCH --array=1-20
#SBATCH --job-name=TomSlim
#SBATCH --time=20:00:00
#SBATCH --mem=6000
#SBATCH --output=runBayPass_BC_map.cor0.5.%A%a.out
#SBATCH --error=runBayPass_BC_map.cor0.5.%A%a.err


/home/booker/bin/baypass_2.2/sources/i_baypass -npop 40 -gfile BC_Map_sampled_cor0.5/${SLURM_ARRAY_TASK_ID}_0.5_192_d40n50_cor0.5.bayPass.txt -efile BC_Map_sampled_cor0.5/${SLURM_ARRAY_TASK_ID}_0.5_192_d40n50_cor0.5.bayPass.pc1 -outprefix BC_Map_sampled_cor0.5/${SLURM_ARRAY_TASK_ID}_0.5_192_d40n50_cor0.5.bayPass
