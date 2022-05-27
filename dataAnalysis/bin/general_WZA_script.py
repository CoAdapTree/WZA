import pandas as pd
import scipy.stats
import numpy as np
import sys, glob, argparse


def WZA( gea , statistic , MAF_filter = 0.05):
## gea - the name of the pandas dataFrame with the gea results
## statistic - the name of the column with your p-values
## MAF_filter - the lowest MAF you wil tolerate
## NOTE, this function assumes that the DataFrame has a column named pbar_qbar
#	print("before sample\n", gea)

## Very small p-values throw Infinities when converted to z_scores, so I convert them to small numbers (i.e. 1e-15)

	gea[statistic] = gea[statistic].clip(lower = 1e-15 )

	gea[statistic] = gea[statistic].replace(1 , 1-1e-3)

# convert the p-values into 1-sided Z scores (hence the 1 - p-values)
	gea["z_score"] = scipy.stats.norm.ppf(1 - np.array( gea[statistic], dtype = float))

	gea["pbar_qbar"] = gea["MAF"]*(1-gea["MAF"])

## Apply the MAF filter
	gea_filt = gea[ gea["MAF"] > MAF_filter ].copy()

## Calculate the numerator and the denominator for the WZA

	gea_filt["weiZ_num"] = gea_filt["pbar_qbar"] * gea_filt["z_score"]

	gea_filt["weiZ_den"] =  gea_filt["pbar_qbar"]**2

	numerator = gea_filt["weiZ_num"].sum()

	denominator = np.sqrt(gea_filt["weiZ_den"].sum())

## We've calculated the num. and the den., let's take the ratio
	weiZ  = numerator/denominator

## Return the final dataframe
	return weiZ


def top_candidate( gea, thresh):

## gea - the name of the pandas dataFrame with the gea results
## thresh - the p_value threshold for determining hits
## MAF_filter - the lowest MAF you wil tolerate
## prop_hits - the average probility of getting a hit - this should be the quantile threshold
##top_candidate_threshold - the probability point for calculating the expected number of genes
## Identifty the hits
	hits = ( gea["pVal"] < thresh ).sum()

	snps = gea.shape[0]

	top_candidate_p = scipy.stats.binom_test(hits, snps, thresh, alternative = "greater" )


	return top_candidate_p, hits


def main():

## Define command line args
    parser = argparse.ArgumentParser(description="A script that implements the WZA, a method for combining evidence across closely linked SNPs in GEA studies.")

    parser.add_argument("--correlations", "-c",
            required = True,
            dest = "correlations",
            type = str,
            help = "The file containing the correlations")

    parser.add_argument("--summary_stat", "-s",
            required = True,
            dest = "summary_stat",
            type = str,
            help = "The name of the column you are analysing")

    parser.add_argument("--window", "-w",
            required = True,
            dest = "window",
            type = str,
            help = "The name of column containing the windows you want to analyse")

    parser.add_argument("--output",
            required = True,
            dest = "output",
            type = str,
            help = "The name of the output file")

    parser.add_argument("--sample_snps",
            required = False,
            dest = "sample_snps",
            type = int,
            help = "[OPTIONAL] Give the number of SNPs you want to downsample to. Give -1 if you want to use the 75th percentile of the number of SNPs. Note that calculating the median within the script is slow, so you may want to run a dummy analysis, get the median number of SNPs then use that explicitly.",
            default = 0)

    parser.add_argument("--resamples",
            required = False,
            dest = "resamples",
            type = int,
            help = "[OPTIONAL] Give the number of times you want to resample WZA scores when the number of SNPs exceeds the sample_snps threshold. [100]",
            default = 100)

    parser.add_argument("--min_snps",
            required = False,
            dest = "min_snps",
            type = int,
            help = "[OPTIONAL] Give the minimum number of SNPs you are willing to analyse per window. Default is 3",
            default = 3)

    parser.add_argument("--large_i_small_p",
            required = False,
            action = "store_true",
            help = "[OPTIONAL] Give this flag if extreme values of the summary stat you're using are large values.")

    parser.add_argument("--top_candidate_threshold",
            required = False,
            dest = "top_candidate_threshold",
            type = float,
            help = "[OPTIONAL] Give the percentile threshold you want to use for the top-candidate test (from Yeaman et al 2016 - Science)",
            default = 99)

    parser.add_argument("--verbose", "-v",
            required = False,
            action = "store_true",
            help = "[OPTIONAL] Give this flag if you want to run the script in verbose mode")

    parser.add_argument("--MAF",
            required = False,
            dest = "MAF",
            type = str,
            help = "[OPTIONAL] Give the name of the MAF colummn - if it's not explicitly given in the input dataframe")


    args = parser.parse_args()

    csv = pd.read_csv(args.correlations, sep = "\t")

    if args.window not in list(csv):
        print("The window variable you provided is not in the dataframe you gave")
        return
    if args.summary_stat not in list(csv):
        print("The summary statistic variable you provided is not in the dataframe you gave")
        return

    if "MAF" in list(csv):
        pass
    else:
        csv["MAF"] = csv[args.MAF].copy()

    if args.large_i_small_p:
        csv["pVal"] =  1 - csv[args.summary_stat].rank()/csv.shape[0]
    else:
        csv["pVal"] = csv[args.summary_stat].rank()/csv.shape[0]

    csv_genes = csv[csv[args.window]!="none"]

    if args.verbose:
        print("here's a peek at the input data")
        print(csv_genes.head)

    print(list(csv_genes))

    csv_gb_gene = csv_genes.groupby(args.window)

    if args.sample_snps == -1:
## Get a list of the number of SNPs per contig

        csv_gb_gene_SNP_count = csv_genes.groupby(args.window)

        num_SNP_list = np.array( [s[1].shape[0] for s in csv_gb_gene_SNP_count if s[1].shape[0] > args.min_snps]  )

## Calculate the 75th percentile of SNPs per gene
        max_SNP_count = int(np.percentile(num_SNP_list[num_SNP_list!=0], 75))

        if args.verbose:
            print("Using the 75th percentile number of SNPs as the maximum in each gene:", max_SNP_count)

    elif args.sample_snps == 0:
        max_SNP_count = int(1e6) # This is just a large number that is never going to be the number of SNPs within a gene
        if args.verbose:
            print("The maximum number of SNPs in each gene:", max_SNP_count)

    else:
        max_SNP_count = int(args.sample_snps)
        if args.verbose:
            print("The maximum number of SNPs in each gene:", max_SNP_count)


    all_genes = []

    count = 0

    for g in csv_gb_gene:
        count += 1
        gene = g[0]
        gene_df = g[1].copy()

## Perform the WZA on the annotations in the contig using parametric ps

        if gene_df.shape[0] <=max_SNP_count:
## Perform the WZA on the annotations in the contig using empirical ps
            wza = WZA(gene_df, "pVal")
        else:
            wza = np.array( [ WZA(gene_df.sample(max_SNP_count), "pVal") for i in range(args.resamples ) ] ).mean()

#        print(wza)

        top_candidate_p, hits = top_candidate( gene_df, 1-(args.top_candidate_threshold/100))

#        print(top_candidate_p, hits, gene_df.shape[0], 1-(args.top_candidate_threshold/100))
#        print(wza_emp_pVal)

        if args.verbose:
            print("\ngene #:", count, "\tgene:", gene,"\tWZA:",wza,"\tTC:", top_candidate_p)

        output = {"gene": gene,
                    "SNPs":gene_df.shape[0],
                    "hits":hits,
                    "gene":gene,
                    "Z":wza,
                    "top_candidate_p":top_candidate_p
                    }

        all_genes.append( output )
#            input("\nPress Enter to continue...")
    out_DF =  pd.DataFrame( all_genes )
    out_DF.to_csv(args.output, index = False)


main()
