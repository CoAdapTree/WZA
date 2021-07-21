# Analysis of the simulation data

The first thing that I did is to use the tree-sequence files to perform a simple genotype-environment association analysis. The python script [bin/parseSimulations.py](bin/parseSimulations.py) takes a tree-seq file (the output from SLiM) and uses the PySLiM, msprime, tskit workflow to add neutral mutations to the simulated populations. Once that is done, it takes the file of phenotypic optima used to simulate local adaptation and performs Kendall's tau (a simple rank correlation) between allele frequency and the environment for each SNP in the genome.

Here's an example for a single file:

```

python bin/parseSimulations.py --trees 1_0.003_1.12Loci.directionalSelection.trees \
                --optima slim_configs/BC_Map_environments.14x14.txt \
               --output 1_0.003_directional  \
               --nPops 40 \
               --nInds 50 \
               --directional \
               --bayPass

```
This will use the tree-sequence flie ```1_0.003_1.12Loci.directionalSelection.trees``` and the phenotypic optima for the BC Map (```slim_configs/BC_Map_environments.14x14.txt```). GEA will be performed using a sample of 40 demes, with 50 individuals taken in each location to estimate allele frequencies. The ```--directional``` flag is used to rell the script that the simulation modelled directional selection and to determine the genes involved in local adaptation accordingly. The output prefix is given as ```1_0.003_directional```. Since the ```--bayPass``` flag was used, the script will perform a GEA on each SNP *and* use the data to make a config file for BayPass.

I wrote long and difficult to read shell scripts to run these analyses for each of the simulated datasets that I generated.

For example, the script [directionalSelection/run.sh](directionalSelection/run.sh) runs the above script in a variety of ways for the data simulated under the directional selection model.
