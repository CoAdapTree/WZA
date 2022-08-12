import argparse
import pandas as pd

def main():
	parser = argparse.ArgumentParser(description="")

	parser.add_argument("--anno",
		required = True,
		dest = "anno",
		type = str,
		help = "The name of the file containing information on the genes and orthogroups")

	parser.add_argument("--wza",
		required = True,
		dest = "wza",
		type = str,
		help = "The name of the file containing the wza scores per gene")

	parser.add_argument("--species",
		required = True,
		dest = "species",
		type = str,
		help = "The initials of the annotation species you will use")

	parser.add_argument("--output",
		required = True,
		dest = "output",
		type = str,
		help = "Give a name for the output file")

	args = parser.parse_args()


	species_dict = {}

# This dict stores the position index of each species in the annotation file
	species_index = {}
## It is assumed that the annotation file has the following structure
#Df	Pa	Pi
#DF      JP      LP      PA      WL
#	species_list = [ "Df" , "Pa", "Pi"]
	species_list = [ "DF" , "JP", "LP", "PA", "WL"]
## Orthogroup_name     species_1      species_2      species_3      ... species_k

	for focal_species in species_list:
		count = 0

		focal_species_annotation_dict = {}

		for line in open( args.anno ):
			count +=1

			if count == 1:
				for i in range(len(line.strip().split("\t"))):
					## ignore "Orthogroup" in the header
					if i ==0: continue
					## save the index of each species in the
					species_index[ line.strip().split()[i] ] = i
				print( species_index )
				continue

## Need to be specific with split and strip - Orthofinder doesn't give a missing data character
			line_data = line.strip("\n").split("\t")
			if len(line_data) != len(species_index.keys())+1:
				print("missing data in the orthogroup table at line:", count)
				break

			orthogroup_name = line_data[0]

## make a list of all the genes within the current orthogroup for the focal species
			species_anno_list = line_data[ species_index[ focal_species ] ].split(", ")

## For each gene belonging to the current orthogroup, make an entry in the dict
			for s_i in species_anno_list:
				focal_species_annotation_dict[ s_i ] = orthogroup_name

		species_dict[ focal_species ] = focal_species_annotation_dict


	wza_csv = pd.read_csv( args.wza )

	OG_dict = {}
	for index, row in wza_csv.iterrows():
		try:
			OG =  species_dict[ args.species ][row['gene']]  # print(row['gene'])
		except KeyError:
## It would be better if we did not have to do the following, if there was a sinlge name for each gene
			splitter = row['gene'].split(";")[-1]#.split("=")[1])
			print(splitter)
			if len( splitter.split("=") ) == 2:
				try:
					OG =  species_dict[ args.species ][splitter.split("=")[1]]  # print(row['gene'])
				except KeyError:
					OG = "NA"
				print("YP",splitter.split("=")[1])
				print(OG)

			else:
				OG = "NA"
		## do some string manipulation thing here
		if OG != "NA":print(OG)

		OG_dict[row['gene']] = OG


## Flip the dictionary around to get the number of genes in each orthogroup for my species:
	flipped = {}

	for key, value in OG_dict.items():
	    if value not in flipped:
	        flipped[value] = [key]
	    else:
	        flipped[value].append(key)

	number_OG_dict = {}

	for key, value in flipped.items():
		print(key, value)
		number_OG_dict[key] = len(value)
# Set it so that genes not assigned to orthogroups have 0 paralogs
	number_OG_dict["NA"] = 0
## Make a dict that contains the number of members of each orthogroup

	## Now map the OG_dict to the WZA file using pandas
	wza_csv["orthogroup"] = wza_csv["gene"].map(OG_dict)
	wza_csv["num_genes_in_OG"] = wza_csv["orthogroup"].map(number_OG_dict)

	wza_csv.to_csv( args.output , index = False, sep = "\t")
main()
