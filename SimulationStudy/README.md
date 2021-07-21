# Simulation Study

In this directory you'll find the code and some processed data files for performing and analyzing simulations of local adaptation that I implemented to test the performance of the WZA. The actual data files (tree sequences recorded from SLiM) are far too big to store on GitHub. In due course, I'll deposit the tree sequences on Dryad or some other place.

## Running the simulations
The first thing to do is to actually perform the simulations. See the pre-print for a full description of the simulations. The idea of the simulations was to model local adaptation to heterogeneous environments using fairly large genomes. The simulations were performed in SLiM (version 3.4). We used three different maps of environmental heterogeneity - these are specified at the command line.

I ran 20 replicates for each map of environmental heterogeneity using either stabilising or directional selection and 20 replicates evolving neutrally. Additionally, I ran 20 replicates of a population evolving according to the island model. In total there were 160 simulations.

The SLiMulation configuration files are contained in [slim_configs/](slim_configs/).

I ran the simulations on Compute Canada servers, so submitted the simulations to the cluster with a SLURM array job script. Here's what that looked like for simulations assuming the BC map and stabilizing selection:

```
#!/bin/bash
#SBATCH --array=1-20
#SBATCH --job-name=TomSlim
#SBATCH --time=25:00:00
#SBATCH --mem=20000
#SBATCH --output=runDirectionalSelection_BC_map.%A%a.out
#SBATCH --error=runDirectionalSelection_BC_map.%A%a.err

time ~/bin/build/slim -d REP=$SLURM_ARRAY_TASK_ID -d map=1 \
 slim_configs/localAdaptation_stabilisingSelection_chainedGenes.2DsteppingStone.slim

# I included the time command to print a record of how long the simulation took to STDOUT.

```

In this example, you can see that I ran 20 simulation replicates using the script modelling stabilising selection.  I passed the SLURM_ARRAY_TASK_ID to each replicate to keep track of the different reps. I also specified that "map=1", this was my way of using the same script for multiple maps. The map numbers were, 1 for the BC Map, 2 for the Truncated Map and 3 for the Gradient map. The simulations were fairly greedy and took about 25 hours and needed access to about 20Gb of memory.

Each of the maps has a file in the ```slim_configs/``` directory that contains 196 rows, with the phenotypic optimum of each deme in the stepping-stone model. The order of these files is very important if you want to incorporate a particular pattern of spatial heterogeneity. The same files can be used for the island model simulations, but in that case the order is irrelevant.

The SLURM scripts for each set of simulations are given here, but you may need to use something tailored to your set-up.

## Analyzing the simulations
