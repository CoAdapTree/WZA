import pandas as pd
import scipy.stats
import numpy as np
import sys, glob, argparse

## Here's a script to calculate the Weighted Z score (and the top-candidate test from the CoAdapTree data structures)


def WZA( gea , statistic , MAF_filter = 0.05, varName = "Z"):
## gea - the name of the pandas dataFrame with the gea results
## statistic - the name of the column with your p-values
## MAF_filter - the lowest MAF you wil tolerate
## varName - the column name for the weigted Z results

## Very small p-values throw Infinities when converted to z_scores, so I convert them to small numbers (i.e. 1e-15)

	gea[statistic] = gea[statistic].clip(lower = 1e-15 )

	gea[statistic] = gea[statistic].replace(1 , 1-1e-3)

## convert the p-values into 1-sided Z scores (hence the 1 - p-values)
	gea["z"] = scipy.stats.norm.ppf(1 - gea[statistic])

## Apply the MAF filter
#	gea_filt = gea[ gea["maf"] > MAF_filter ].copy()

	gea_filt = gea.copy()

## Calculate the numerator and the denominator for the WZA

	gea_filt[varName+"weiZ_num"] = gea_filt["pbar_qbar"] * gea_filt["z"]

	gea_filt[varName+"weiZ_den"] =  gea_filt["pbar_qbar"]**2 	

	numerator = gea_filt.groupby(["gene"])[varName+"weiZ_num"].sum().to_frame()

	denominator = np.sqrt(gea_filt.groupby(["gene"])[varName+"weiZ_den"].sum()).to_frame()

## We've calculated the num. and the den., let's make a dataframe that has both 
	weiZ  = pd.concat([numerator,denominator], axis = 1, sort = False)

## Actually calculate the Z scores for each gene
	weiZ[varName] = weiZ[varName+"weiZ_num"] / weiZ[varName+"weiZ_den"]

## One might be interested in calculating a p_value from the Z-scores (though this only works if the data are normal, which they won't be if there's populations structure).
	weiZ[varName+"_hits"] = (weiZ[varName] > scipy.stats.norm.ppf(1 - 0.05/50)).astype(int)

## Return the final dataframe
	return weiZ


## A function for performing the top-candidate test

def top_candidate( gea, thresh, statistic, prop_hits, MAF_filter = 0.05):

## gea - the name of the pandas dataFrame with the gea results
## thresh - the p_value threshold for determining hits
## MAF_filter - the lowest MAF you wil tolerate
## prop_hits - the average probility of getting a hit

## Identifty the hits
	gea["hits"] = ( -np.log10(gea[statistic]) > thresh).astype(int)

## Apply the MAF filter
#	gea_filt = gea[ gea["maf"] > MAF_filter ]
	gea_filt = gea.copy()

## Count the hits per gene
	num_hits = gea_filt.groupby(["gene"])["hits"].sum().to_frame()

## Count the SNPs per gene
	num_SNPs = gea_filt.groupby(["gene"])["hits"].count().to_frame()

## Make a single DF with the hits and the SNPs
	TC  = pd.concat([num_hits, num_SNPs], axis = 1, sort = False) 

## Name the cols
	TC.columns = ["hits", "SNPs"]

## Init an empty vector for p_values (the top-candidate index)
	p_vals = []
	for index, row in TC.iterrows():
		p_vals.append( scipy.stats.binom_test(row.hits, row.SNPs, prop_hits, alternative = "greater" ) )

## Add the TC p_vals to the genes
	TC["top_candidate_p"] = p_vals

## Return the resulting dataFrame
	return(TC)


class corLine:

	def __init__(self, c):

		dat = c.strip().split()

		self.contig = dat[0].split("-")[0]

		self.pos = int(dat[0].split("-")[1])

## The data in the Correlations file are stored as contig-pos, so I'll parse that out

		self.env = dat[1]

		self.rho = float(dat[2])

		self.pVal = float(dat[3])


## A generator function to spit out the consecutive lines for the contig
def contigGenerator(correlationFile, envFilter = "MAT"):

	current_contig = ""

	contig_data = []

	with open( correlationFile ) as cor:
		for c in cor:
## Ignore the header
			if c.startswith("snp"):continue
			
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

		if sum(gene_bool) > 1:
			print("more than one gene overlapping with SNP:", s.contig, s.pos)

		data_for_table.append({"contig":s.contig,
									"pos":s.pos,
									"rho":s.rho,
									"pVal":s.pVal,
									"gene":gene_name})

	return pd.DataFrame(data_for_table)	

def main():

## Define command line args

	parser = argparse.ArgumentParser(description="")

	parser.add_argument("--correlations", "-c", 

			required = True,

			dest = "correlations",

			type = str, 

			help = "The file containing the correlations")

	parser.add_argument("--GFF","-g", 

			required = True,

			dest = "GFF", 

			help = "The GFF file (gene annotations)?")

	parser.add_argument("--output", 

			required = True,

			dest = "output",

			type = str, 

			help = "The name of the output file")

	args = parser.parse_args()
	
## GFF header from ENSEMBL webpage
	gff = pd.read_csv(args.GFF , 

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

	all_contigs = []

	for contig,SNPs in contigGenerator(args.correlations):

		contigDF = contigSnpTable(SNPs, gff[gff["seqname"] == contig])

		if contigDF.shape == (0,0):
			continue

		contigDF["pbar_qbar"] = 0.25

		position = contigDF.groupby(["gene"])["pos"].mean().to_frame()

		wza = WZA(contigDF, "pVal")

		TopCan = top_candidate( contigDF, 0.05, "pVal", 0.05, MAF_filter = 0.05 )

		result = pd.concat([ position, wza, TopCan] , axis = 1, sort = True ).reset_index()

		result["contig"] = contig

		all_contigs.append( result )

	pd.concat(all_contigs).to_csv(args.output,index = False)

main()
