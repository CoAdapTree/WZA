# Simulation Study

In this directory you'll find the code and some processed data files for performing and analyzing simulations of local adaptation that I implemented to test the performance of the WZA. The actual data files (tree sequences recorded from SLiM) are far too big to store on GitHub. In due course, I'll deposit the tree sequences on Dryad or some other place.

The first thing to do is to actually perform the simulations. See the pre-print for a full description of the simulations. The idea of the simulations was to model local adaptation to heterogeneous environments using fairly large genomes. The simulations were performed in SLiM (version 3.4). We used three different maps of environmental heterogeneity - these are specified at the command line.

The SLiMulation configuration files are contained in [slim_configs/](slim_configs/).

I ran the simulations on Compute Canada servers, so submitted the simulations to the cluster with a SLURM array job script. Here's what that looked like for simulations assuming the BC map:

```
#!/bin/bash
#SBATCH --array=1-20
#SBATCH --job-name=TomSlim
#SBATCH --time=25:00:00
#SBATCH --mem=20000
#SBATCH --output=runDirectionalSelection_BC_map.%A%a.out
#SBATCH --error=runDirectionalSelection_BC_map.%A%a.err

time ~/bin/build/slim -d REP=$SLURM_ARRAY_TASK_ID -d map=2  ~/projects/def-whitlock/booker/GEA/slim_configs/localAdaptation_stabilisingSelection_chainedGenes.2DsteppingStone.slim

```
I included the UNIX ```time``` command to print a record of how long the simulation took to STDOUT. 
