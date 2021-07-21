# Running BayPass on the simulation data

I ran BayPass on the simulated data as it is a widely used GEA package that controls for the confounding effects of population structure.

I downloaded BayPass from: [http://www1.montpellier.inra.fr/CBGP/software/baypass/](http://www1.montpellier.inra.fr/CBGP/software/baypass/)

I ran the BayPass analysis on Compute Canada servers, here's an example of how I ran the program, following worked example 5.1.2 in the manual that comes with BayPass:

```bash
#!/bin/bash
#SBATCH --array=1-20
#SBATCH --job-name=TomSlim
#SBATCH --time=20:00:00
#SBATCH --mem=6000
#SBATCH --output=runBayPass_BC_map.%A%a.out
#SBATCH --error=runBayPass_BC_map.%A%a.err


/home/booker/bin/baypass_2.2/sources/i_baypass -npop 40 \
      -gfile BC_Map_sampled/${SLURM_ARRAY_TASK_ID}_0.5_192_d40n50_i10000.bayPass.txt \
      -efile BC_Map_sampled/${SLURM_ARRAY_TASK_ID}_0.5_192_d40n50_i10000.bayPass.pc1 \
      -outprefix BC_Map_sampled/${SLURM_ARRAY_TASK_ID}_0.5_192_d40n50_i10000.bayPass
```

When I was working on the server, I loaded the data in using a different architechture than I used for this repository, so, again, the paths would need to be corrected to run on another machine. 
