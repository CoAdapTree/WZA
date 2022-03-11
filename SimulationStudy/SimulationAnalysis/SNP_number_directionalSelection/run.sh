
## Run the WZA on the output of the GEA analyses, but use windows that are 2, 4 or 8 SNPs long

## Takes a little while...
parallel "python ../../../bin/windowSizeExploration.py --csv ../s0.003/BC_Map_sampled/ --output BC_Map_sampled.WZA.n{}.csv --nSNPs {}" ::: 2 4 8

parallel "python ../../../bin/windowSizeExploration.py --csv ../s0.003/cline_sampled/ --output cline_sampled.WZA.n{}.csv --nSNPs {}" ::: 2 4 8

parallel "python ../../../bin/windowSizeExploration.py --csv ../s0.003/trunc_sampled/ --output trunc_sampled.WZA.n{}.csv --nSNPs {}" ::: 2 4 8


## Summarise the results:

parallel "python ../../../bin/summariseWZA.py --wza {1}_sampled.WZA.n{2}.csv --LA 0.005 --output {1}_sampled.WZA.n{2} --SNP_intervals" ::: BC_Map cline trunc ::: 2 4 8

## I then plot the reults using the  Plots/SNP_number.R script
