import pandas as pd
import scipy.stats
import numpy as np
import sys, glob, argparse

def calculateGeneMidpoints(SNPs_in_genes, gene_annotations, posName):

    genes_with_SNPs = set(SNPs_in_genes["attribute"])
    genes_with_SNPs_waypoints =  SNPs_in_genes.groupby('attribute')[posName].median().to_dict()

    genes_without_SNPs_slice = gene_annotations[~gene_annotations["attribute"].isin(genes_with_SNPs) ]
    genes_without_SNPs_waypoints =  pd.Series(genes_without_SNPs_slice.midpoint.values,index=genes_without_SNPs_slice.attribute).to_dict()

#    print("genes with SNPs", genes_with_SNPs_waypoints)
#    print("genes without SNPs", genes_without_SNPs_waypoints)
#    print()
    genes_with_SNPs_waypoints.update(genes_without_SNPs_waypoints)
    return( genes_with_SNPs_waypoints )



def calcGeneDistances(snp_series, waypoints ):

# Make an array of the gene waypoints
    gene_pos_array = np.array(list(waypoints.values()))

# Make an array of the gene waypoints
    gene_name_array = np.array(list(waypoints.keys()))

# For each SNP, calculate the distance from each waypoint
    dist_arrays =  snp_series.apply( lambda x: abs(x-gene_pos_array))

# For each SNP, get the index of the nearest array
    min_dist_index_series = dist_arrays.apply( lambda r: np.argmin(r))

# For each SNP, get the index of the nearest array
    min_dist_index_series = dist_arrays.apply( lambda r: np.argmin(r))

# For each SNP, get the minimum distance
    gene_dists = [ e[g] for e, g in zip(dist_arrays, min_dist_index_series)]

# For each SNP, get the name of the gene with the minium distance
    gene_names = [ gene_name_array[g] for g in  min_dist_index_series]

    return gene_dists, gene_names

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

    parser.add_argument("--bed",
    required = False,
    dest = "bed",
    action = "store_true",
    help = "[OPTIONAL] Give this flag if the analysis files are in BED format. Otherwise the script assumes GFF format")


    parser.add_argument("--Fst",
    required = False,
    dest = "Fst",
    action = "store_true",
    help = "[OPTIONAL] Is this an Fst file?")
    args = parser.parse_args()

    if args.bed:
        annotations = pd.read_csv(args.annotations ,
        sep = "\t",
        names = ["seqname",
        "start",
        "end",
        "attribute"])
        ## Add 1 to the positions to make correct for 0-based BedTools
#        annotations["start"] +=1
#   Commented out because Pooja did not make the BED files 0-based
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

    annotations["midpoint"] =  (annotations["start"] +  annotations["end"])/2

#    header = ["contig","pos","env","rho","pVal","MAF","emp_pVal"]
    if args.Fst:
        csv = pd.read_csv(args.correlations)
    else:
        csv = pd.read_csv(args.correlations,
#                        header = None,
                        delim_whitespace=True,
#                        names = header,
                        index_col=False)

    if args.Fst:
        chrom_name = "CHROM"
        pos_name = 'POS'
    else:
        chrom_name = "contig"
        pos_name = 'pos'
## Group the SNP DF by contig
    csv_gb = csv.groupby(chrom_name)

    output_DFs = []
    count = 0
    for i in csv_gb:
        count +=1
    ## Grab the chunk of the annotation dict that matches
        relevant_annotations = annotations[annotations["seqname"] == i[0]].copy()

#        if relevant_annotations.shape[0] !=0:
#            print(relevant_annotations)
#            print(i[0])
#            print(i[1])
#            print(relevant_annotations)
#            print()
        SNP_df = i[1].copy()
        ii = pd.IntervalIndex.from_arrays(relevant_annotations['start'], relevant_annotations['end'],
        closed="neither")

        relevant_annotations['Interval'] = ii
        relevant_annotations['length'] = relevant_annotations.end-relevant_annotations.start

## The following condition will be met if there are no genes on the contig/chromosome
        if relevant_annotations.shape[0] == 0:
            out = SNP_df
            out["attribute"] = "NA"
#            print(out)
            out["geneDist"] = "NA"
            out["geneName"] = "NA"
            output_DFs.append(out)
        else:
## This event will be triggered if there are genes on the relevant contig/chromosome
            out = SNP_df.assign(Interval=pd.cut(SNP_df[pos_name], bins=ii))

            ## SNPs with no genes
            noGeneSNPs = out.loc[~out.index.isin(out.dropna(subset=['Interval']).index)].drop(columns='Interval')
            noGeneSNPs["attribute"] = "none"
            noGeneSNPs["length"] = 0
#            print(noGeneSNPs)

            ## SNPs with genes
            geneSNPs = out.dropna(subset=['Interval'])
            geneSNPs = geneSNPs.merge(relevant_annotations[['attribute', "length",'Interval']], on='Interval', how='left') \
            .drop(columns='Interval')


            geneMidpoints = calculateGeneMidpoints(geneSNPs, relevant_annotations,pos_name)

            if len(geneMidpoints.keys()) != relevant_annotations.shape[0]:
                print("What the fuck?")
                return

            contig_SNPs = pd.concat([noGeneSNPs, geneSNPs ])

            geneDists, geneNames = calcGeneDistances( contig_SNPs[pos_name], geneMidpoints )

            contig_SNPs["geneDist"] = geneDists
            contig_SNPs["geneName"] = geneNames
            output_DFs.append(contig_SNPs)


    pd.concat(output_DFs).to_csv(args.output, index = False)



main()
