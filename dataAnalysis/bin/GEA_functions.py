import pandas as pd
import scipy.stats
import numpy as np
from GEA_functions import *

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

## We've calculated the num. and the den., let's make a dataframe that has both
	weiZ  = numerator/denominator

## Deprecated
## One might be interested in calculating a p_value from the Z-scores (though this only works if the data are normal, which they won't be if there's populations structure).
##	weiZ[varName+"_hits"] = (weiZ[varName] > scipy.stats.norm.ppf(1 - 0.05/50)).astype(int)

## Return the final dataframe
	return weiZ


def WZA_group( gea , statistic , MAF_filter = 0.05, varName = "Z", SNP_count = 1e6):
## gea - the name of the pandas dataFrame with the gea results
## statistic - the name of the column with your p-values
## MAF_filter - the lowest MAF you wil tolerate
## varName - the column name for the weigted Z results
## NOTE, this function assumes that the DataFrame has a column named pbar_qbar
#	print("before sample\n", gea)
	if SNP_count != 1e6:
		gea_gby = gea.groupby(["gene"])
		gea = gea_gby.apply(lambda x: x.sample(n= SNP_count) if len(x)> SNP_count else x).reset_index(drop=True)

## Very small p-values throw Infinities when converted to z_scores, so I convert them to small numbers (i.e. 1e-15)

	gea[statistic] = gea[statistic].clip(lower = 1e-15 )

	gea[statistic] = gea[statistic].replace(1 , 1-1e-3)

# convert the p-values into 1-sided Z scores (hence the 1 - p-values)
	gea["z"] = scipy.stats.norm.ppf(1 - np.array( gea[statistic], dtype = float))

## Apply the MAF filter
#	gea_filt = gea[ gea["maf"] > MAF_filter ].copy()

	gea_filt = gea.copy()

## Calculate the numerator and the denominator for the WZA

	gea_filt[varName+"_weiZ_num"] = gea_filt["pbar_qbar"] * gea_filt["z"]

	gea_filt[varName+"_weiZ_den"] =  gea_filt["pbar_qbar"]**2

	numerator = gea_filt.groupby(["gene"])[varName+"_weiZ_num"].sum().to_frame()

	denominator = np.sqrt(gea_filt.groupby(["gene"])[varName+"_weiZ_den"].sum()).to_frame()

## We've calculated the num. and the den., let's make a dataframe that has both
	weiZ  = pd.concat([numerator,denominator], axis = 1, sort = False)

## Actually calculate the Z scores for each gene
	weiZ[varName] = weiZ[varName+"_weiZ_num"] / weiZ[varName+"_weiZ_den"]

## Deprecated
## One might be interested in calculating a p_value from the Z-scores (though this only works if the data are normal, which they won't be if there's populations structure).
##	weiZ[varName+"_hits"] = (weiZ[varName] > scipy.stats.norm.ppf(1 - 0.05/50)).astype(int)

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


## A function for performing the top-candidate test on grouped dataframes

def top_candidate_group( gea, thresh, threshQuant, statistic, top_candidate_threshold, MAF_filter = 0.05):

## gea - the name of the pandas dataFrame with the gea results
## thresh - the p_value threshold for determining hits
## MAF_filter - the lowest MAF you wil tolerate
## prop_hits - the average probility of getting a hit - this should be the quantile threshold
##top_candidate_threshold - the probability point for calculating the expected number of genes
## Identifty the hits
	gea["hits"] = ( -np.log10( np.array( gea[statistic] , dtype = float)) > -np.log10(thresh)).astype(int)

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

	return(TC)


### Init an empty vector for expected hits at the "top-candidate" threshold
#	expectedHits = []
#	for index, row in TC.iterrows():
#		print(row.hits, row.SNPs, top_candidate_threshold )
#		p_vals.append( scipy.stats.binom_test(row.hits, row.SNPs, threshQuant, alternative = "greater" ) )
#		expectedHits.append( scipy.stats.binom.ppf( top_candidate_threshold , row.SNPs, threshQuant ) )


## Add the TC p_vals to the genes DF
#	TC["top_candidate_p"] = p_vals

## Add the TC expected number of hits to the genes DF
#	TC["expected_outliers"] = expectedHits


## Return the resulting dataFrame
#	return(TC)



# A function for
#

def WZA_SNP_bins(input_gea, wza_column, snp_column, num_bins):

#   WZA_SNP_bins(out_DF, wza_col = "Z_empP", snp_col = "SNPs"):
	print(" hi tommy")

	num_SNP_list = np.array( input_gea[snp_column] )
	print(num_SNP_list.max())

# I use Numpy's histogram function to get the bin edges
	bin_edges = np.histogram( num_SNP_list, bins = num_bins )[1]

	bin_edges[-1] = 1e6


	intervals =  pd.IntervalIndex.from_breaks(bin_edges, closed = "left")

	bin_names =np.array([str(k) for k in bin_edges ][0:-1])

	gene_bins =  [ bin_names[intervals.contains(n)][0] for n in num_SNP_list ]


#	return
	input_gea["bin"] = gene_bins

	approx_pVals = []
	for gs in input_gea.groupby("bin"):
		mean_Z = gs[1][wza_column].mean()
		sd_Z = gs[1][wza_column].std()
		temp =  -1*np.log10( 1-scipy.stats.norm.cdf( gs[1][wza_column], mean_Z, sd_Z) )
		new = gs[1].copy()
		new["approx_pVal"] = temp
		approx_pVals.append( new )
		print( new )

	return( pd.concat(approx_pVals))
