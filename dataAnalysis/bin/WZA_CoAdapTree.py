import pandas as pd
import scipy.stats
import numpy as np
import sys, glob, argparse
from GEA_functions import *

## Here's a script to calculate the Weighted Z score (and the top-candidate test from the CoAdapTree data structures)

class corLine:

	def __init__(self, c):

		dat = c.strip().split()

		self.contig = dat[0]

		self.pos = int(dat[1])

## The data in the Correlations file are stored as contig-pos, so I'll parse that out
## Note that this is VERY sensitive to the file structure staying the same
		self.env = dat[2]

		self.rho = float(dat[3])

		self.pVal = float(dat[4])

		MAF = float(dat[5])

		self.pbar_qbar = MAF * (1 - MAF)

		self.emp_pVal = float(dat[6])


## A generator function to spit out the consecutive lines for the contig
def contigGenerator(correlationFile, envFilter):

	current_contig = ""

	contig_data = []

	with open( correlationFile ) as cor:
		for c in cor:
## Ignore the headera
			if c.startswith("snp") or c.startswith("contig"):continue

			currentLine = corLine(c)

			if currentLine.env != envFilter: continue

			if currentLine.contig != current_contig:

				if current_contig == "":

					current_contig = currentLine.contig

					contig_data.append( currentLine )

				else:

					yield current_contig, contig_data

					contig_data = [currentLine]

					current_contig = currentLine.contig

			elif currentLine.contig == current_contig:

				contig_data.append( currentLine )

	yield current_contig, contig_data

def contigSnpTable(snps, contig_dat):

	data_for_table = []

	for s in snps:

		gene_bool = (s.pos>=contig_dat.start)&(s.pos<=contig_dat.end)

		if sum(gene_bool) == 0:continue

		gene_name = ":".join(list(contig_dat.attribute[gene_bool]))

		gene_start = list(contig_dat.start[gene_bool] )[0]

		gene_end = list( contig_dat.end[gene_bool] )[0]

		if sum(gene_bool) > 1:
			print("more than one gene overlapping with SNP:", s.contig, s.pos)

		data_for_table.append({"contig":s.contig,
									"pos":s.pos,
									"rho":s.rho,
									"pVal":s.pVal,
									"empirical_pVal":s.emp_pVal,
									"pbar_qbar":s.pbar_qbar,
									"gene":gene_name,
									"gene_start":gene_start,
									"gene_end":gene_end})

	return pd.DataFrame(data_for_table)


def correlationThreshold( corData, targetEnv, percentile_threshold = 99.9):
## Make an empty container to dump the pVals into
	pValues = []
	with open( corData ) as cor:
		for c in cor:
			if c.startswith("snp") or c.startswith("contig"):continue
			currentLine = corLine(c)
#			print( currentLine.pVal )
			if currentLine.env != targetEnv: continue
			else:
				pValues.append(currentLine.pVal)
#	print( pValues )
	if len(pValues) == 0:
		return None
	else:
		return( 1 - np.percentile( 1 - np.array(pValues), percentile_threshold ) )
#		return( np.percentile( np.array(pValues), 100-percentile_threshold ) )


def main():

## Define command line args

	parser = argparse.ArgumentParser(description="")

	parser.add_argument("--correlations", "-c",

			required = True,

			dest = "correlations",

			type = str,

			help = "The file containing the correlations")

	parser.add_argument("--annotations","-a",

			required = True,

			dest = "annotations",

			help = "The file of annotations. The script asssumes GFF as default, set the '--bed' flag if using that format")

	parser.add_argument("--output",

			required = True,

			dest = "output",

			type = str,

			help = "The name of the output file (the environment will be prepended to the file name so be sure to write to this dir!)")

	parser.add_argument("--env",

			required = False,

			dest = "env",

			type = str,

			help = "[OPTIONAL] If you want to analyse just a single environment, give it here")

	parser.add_argument("--bed",

			required = False,

			dest = "bed",

			action = "store_true",

			help = "[OPTIONAL] Give this flag if the analysis files are in BED format. Otherwise the script assumes GFF format")

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

	if args.bed:
		annotations = pd.read_csv(args.annotations ,

					sep = "\t",

						names = ["seqname",

								"start",

								"end",

								"attribute"])
## Add 1 to the positions to make correct for 0-based BedTools
		annotations["start"] +=1

	else:
## GFF header from ENSEMBL webpage
		annotations = pd.read_csv(args.annotations ,

					sep = "\t",

						names = ["seqname",

								"source",

								"feature",

								"start",

								"end",

								"score",

								"strand",

								"frame",

								"attribute"])

## For all environmental variables:
## MAT MWMT MCMT TD MAP MSP AHM SHM DD_0 DDS NFFD bFFP FFP PAS EMT EXT Eref CMD

	envs =["MAT","MWMT","MCMT","TD",
			"MAP","MSP","AHM","SHM",
			"DD_0","DD5","NFFD","bFFP",
			"FFP","PAS","EMT","EXT","Eref","CMD",
			"Latitude", "Elevation"]
## If you specify just a single environment (which is recommended for parallelisation)
	if args.env:
## Make sure that the environment specified is actually valid
		if args.env not in envs:
			print("you did not specify a valid environment. Take a look at your command, bozo")
			return
		envs = [args.env]

## We're going to analyse each env. separately
	for env in envs:
		all_contigs = []

## Now let's calculate the outlier threshold (for the TC test) from the data - Make sure it's a percentile!
		threshold_99th = correlationThreshold( args.correlations, env, percentile_threshold = 99)
		threshold_95th = correlationThreshold( args.correlations, env, percentile_threshold = 95)
		threshold_90th = correlationThreshold( args.correlations, env, percentile_threshold = 90)
		threshold_50th = correlationThreshold( args.correlations, env, percentile_threshold = 50)
		if threshold_99th == None:
			print("Something went wrong when identifying the outlier threshold")
			return
		print("99th percentile:",threshold_99th)
		print("95th percentile:",threshold_95th)
		print("90th percentile:",threshold_90th)
		print("50th percentile:",threshold_50th)

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


		print("Analysing:",env)



## Iterate over contigs spat out by the
		for contig,SNPs in contigGenerator(args.correlations, env):
## For testing:
#			if contig != "jcf7190000000051":continue
## Grab all the genes present on this contig
			contigDF = contigSnpTable(SNPs, annotations[annotations["seqname"] == contig])

		#	print(contigDF)
## If there are no annotations on the current contig, move to the next
			if contigDF.shape == (0,0):
				continue

## Get the average position of each annotation - not used for the analysis, just for downstream plotting

			position = contigDF.groupby(["gene"])["pos"].mean().to_frame()

			anno_start = contigDF.groupby(["gene"])["gene_start"].mean().to_frame()

			anno_end = contigDF.groupby(["gene"])["gene_end"].mean().to_frame()


## Perform the WZA on the annotations in the contig using parametric ps
			wza_pVal = WZA(contigDF, "pVal", varName = "Z_p", SNP_count = max_SNP_count)
#			print(wza_pVal)
## Perform the WZA on the annotations in the contig using empirical ps
#			for k in range(resamples):
#				print(k)
			wza_emp_pVal_samples = [WZA(contigDF, "empirical_pVal", varName = "Z_empP", SNP_count = max_SNP_count) for  k in range(args.resamples)]

			wza_concat = pd.concat( wza_emp_pVal_samples )

			wza_by_row_index = wza_concat.groupby(wza_concat.index)

			wza_df_means = wza_by_row_index.mean()

			print( wza_df_means )

## Perform the top-candidate for the annotations in the contig
			TopCan = top_candidate( contigDF, threshold_99th, 0.01, "pVal", 0.9999, MAF_filter = 0.05 )

## Combine the results
			result = pd.concat([ position, wza_pVal, wza_df_means, TopCan, anno_start, anno_end] , axis = 1, sort = True ).reset_index()

			result["contig"] = contig

			result["env"] = env

			result["max_WZA_snps"] = max_SNP_count
			all_contigs.append( result )

#			input("\nPress Enter to continue...")

		if len( all_contigs ) == 0: continue

## Combine all contig-specific dataframes into a single big one
		outputDF = pd.concat(all_contigs)

## Get the expected proportion of hits per hit-bearing gene
		expected = outputDF[outputDF["hits"]!=0]["hits"].sum() / outputDF[outputDF["hits"]!=0]["SNPs"].sum()

		print("expected proportion", expected)

		top_candidate_p = [ scipy.stats.binom_test(h, s, expected, alternative = "greater" ) for h, s in zip(outputDF.hits, outputDF.SNPs)]

		outputDF["top_candidate_p"] = top_candidate_p

#		expected_hits = [ scipy.stats.binom.ppf( 0.9999 , s, expected )  for s in outputDF.SNPs]

		outputDF["expected_hits"] =  scipy.stats.binom.ppf( 0.9999 , outputDF.SNPs, expected )


## Write the dataframe to an output file
		outputDF.to_csv(env + "_" + args.output,index = False)


main()
