## Analysing CoAdapTree data with the WZA

Here's a little bit of test data to analyse CoAdapTree GEA results.

For each annotation in a specified file (either GFF or BED), the script ```bin/WZA_CoAdapTree.py``` performs both the WZA and the Top-Candidate test (Yeaman et al 2016) on the results from GEA.

For each set of GEA analyses performed on individual SNPs, we need to calculated the empirical pValue. The empirical p-value is just the rank of each SNP divided by the total number of SNPs. E.g. the smallest p-value would be 1/k, where k is the total number of SNPs tested.

The following little script adds empirical p-values to the analysis files:

```
python bin/calculateEmpiricalP.py --correlations file_of_per_SNP_correlation_stats.txt --output file_of_per_SNP_correlation_stats_empiricalP.txt

```

The ```calculateEmpiricalP.py``` will make a dataframe for each environment, with empirical p-Values added as an additional column.

Once you've added the empirical-p values to the individual SNP files, you can run the WZA on them. Here's that step for analysing a particular environment:

```
# The first step is to add gene names to the SNP files for each environment

python bin/annotateSNPfiles.py \
    --correlations file_of_per_SNP_correlation_stats_empiricalP.txt \
    --annotations gff.sorted.filteredforWZA.formatted \
    --output file_of_per_SNP_correlation_stats_empiricalP.txt_geneNames \
    --bed

# The next step is to perform the WZA on the annotated files
python bin/WZA_CoAdapTree_annotatedSNPtables.py \
    --correlations file_of_per_SNP_correlation_stats_empiricalP.txt_geneNames \
    --output WZA.csv \
    --sample_snps 20 \
    --resamples 100
```

____________________
## Deprecated



```
python bin/WZA_CoAdapTree.py --correlations MAT_file_of_per_SNP_correlation_stats_empiricalP.txt --annotations DFtransv20202_DFrefedit.filtered.genes.collapsed.bed --output DF_interior_AHM.csv --env MAT --bed
```

However, there is good reason to make sure that each window has a fairly similar number of SNPs. High heterogeniety in the number of SNPs may generate erratic results from the WZA, so I (TRB) added support for SNP sampling to the WZA scripts.

It is simple to perform. If you have a particular number of SNPs that you want to downsample to, simply add the flag ```--sample X``` to the ```python bin/WZA_CoAdapTree.py``` command above, where ```X``` is maximum number of SNPs you want to analyse for each gene. What to pick for this number? We want to make as much use of the data we have as possible, so one strategy would be to use the median number of SNPs. That way, you are using the exact same number of SNPs for half the genes. I added specific support for using the median to the script. To use the median give ```--sample -1```.
