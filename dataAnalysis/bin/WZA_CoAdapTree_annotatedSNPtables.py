import pandas as pd
import scipy.stats
import numpy as np
import sys, glob, argparse
from GEA_functions import *

def correlationThreshold( corData, targetEnv, percentile_threshold = 99.9):
## Make an empty container to dump the pVals into
    pValues = []
    with open( corData ) as cor:
        for c in cor:
            if c.startswith("snp") or c.startswith("contig"):continue
            currentLine = corLine(c)
#            print( currentLine.pVal )
            if currentLine.env != targetEnv: continue
            else:
                pValues.append(currentLine.pVal)
#    print( pValues )
    if len(pValues) == 0:
        return None
    else:
        return( 1 - np.percentile( 1 - np.array(pValues), percentile_threshold ) )
#        return( np.percentile( np.array(pValues), 100-percentile_threshold ) )

def main():

## Define command line args
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("--correlations", "-c",
            required = True,
            dest = "correlations",
            type = str,
            help = "The file containing the correlations")

    parser.add_argument("--output",
            required = True,
            dest = "output",
            type = str,
            help = "The name of the output file (the environment will be prepended to the file name so be sure to write to this dir!)")

    parser.add_argument("--sample_snps",
            required = False,
            dest = "sample_snps",
            type = int,
            help = "[OPTIONAL] Give the number of SNPs you want to downsample to. Give -1 if you want to use the median number of SNPs. Note that calculating the median within the script is slow, so you may want to run a dummy analysis, get the median number of SNPs then use that explicitly.",
            default = 0)

    parser.add_argument("--resamples",
            required = False,
            dest = "resamples",
            type = int,
            help = "[OPTIONAL] Give the number of times you want to resample WZA scores when the number of SNPs exceeds the sample_snps threshold. [100]",
            default = 100)

    args = parser.parse_args()

    csv = pd.read_csv(args.correlations, sep = "\t")

    print(csv)

    if args.sample_snps == -1:
## Get a list of the number of SNPs per contig
        num_SNP_list = np.array([contigSnpTable(SNPs, annotations[annotations["seqname"] == contig]).shape[0] for contig,SNPs in contigGenerator(args.correlations, env)])
## Calculate the median number of SNPs per gene
        max_SNP_count = int(np.median(num_SNP_list[num_SNP_list!=0]))
        print("Using the median number of SNPs as the maximum in each gene:", max_SNP_count)

    elif args.sample_snps == 0:
        max_SNP_count = int(1e6) # This is just a large number that is never going to be the number of SNPs within a gene
        print("The maximum number of SNPs in each gene:", max_SNP_count)

    else:
        max_SNP_count = int(args.sample_snps)
        print("The maximum number of SNPs in each gene:", max_SNP_count)

    csv_genes = csv[csv["attribute"]!="none"]
    csv_gb_gene = csv_genes.groupby("attribute")

    all_genes = []

    count = 0

    for g in csv_gb_gene:
        count += 1
        gene = g[0]
        gene_df = g[1].copy()

## Perform the WZA on the annotations in the contig using parametric ps
        wza_pVal = WZA(gene_df, "pVal", varName = "Z_p")

        if gene_df.shape[0] <=max_SNP_count:
## Perform the WZA on the annotations in the contig using empirical ps
            wza_emp_pVal = WZA(gene_df, "emp_pVal", varName = "Z_empP")
        else:
            wza_emp_pVal = np.array( [ WZA(gene_df.sample(max_SNP_count), "emp_pVal", varName = "Z_empP") for i in range(args.resamples ) ] ).mean()

#        print(wza_emp_pVal)

        print( count )

        output = {"contig": list(gene_df.contig)[0],
                    "env":list(gene_df.env)[0],
                    "SNPs":gene_df.shape[0],
                    "gene":gene,
                    "Z_empP":wza_emp_pVal,
                    "Z_para":wza_pVal,
                    "pos":gene_df.pos.mean()}

        all_genes.append( output )

#            input("\nPress Enter to continue...")
    out_DF =  pd.DataFrame( all_genes )
    out_DF.to_csv(args.output, index = False)


main()
