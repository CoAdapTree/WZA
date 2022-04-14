import pandas as pd
import scipy.stats
import numpy as np
import sys, glob, argparse


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



    header = ["contig","pos","env","rho","pVal","MAF","emp_pVal"]
    csv = pd.read_csv(args.correlations,
                        header = None,
                        delim_whitespace=True,
                        names = header,
                        index_col=False)
## Group the SNP DF by contig
    csv_gb = csv.groupby('contig')

    output_DFs = []

    for i in csv_gb:
    ## Grab the chunk of the annotation dict that matches
        relevant_annotations = annotations[annotations["seqname"] == i[0]].copy()

        SNP_df = i[1].copy()
        ii = pd.IntervalIndex.from_arrays(relevant_annotations['start'], relevant_annotations['end'],
        closed="neither")

        relevant_annotations['Interval'] = ii

        if relevant_annotations.shape[0] == 0:
            out = SNP_df
            out["attribute"] = "none"
#            print(out)
            output_DFs.append(out)
        else:


            out = SNP_df.assign(Interval=pd.cut(SNP_df['pos'], bins=ii))

            ## SNPs with no genes
            noGeneSNPs = out.loc[~out.index.isin(out.dropna(subset=['Interval']).index)].drop(columns='Interval')
            noGeneSNPs["attribute"] = "none"
#            print(noGeneSNPs)
            output_DFs.append(noGeneSNPs)

            ## SNPs with genes
            geneSNPs = out.dropna(subset=['Interval'])
            geneSNPs = geneSNPs.merge(relevant_annotations[['attribute', 'Interval']], on='Interval', how='left') \
            .drop(columns='Interval')
#            print(geneSNPs)
            output_DFs.append(geneSNPs)


    pd.concat(output_DFs).to_csv(args.output, sep = "\t", index = False)



main()
