parallel "python ../../../bin/WZA.py --csv ../s0.003/{}_sampled/ --bayPass  ../BayPassResults/s0.003/{}_sampled_BayPass/ --output {}_sampled.ZA.csv --unweighted" ::: BC_Map cline trunc

parallel "python ../../../bin/summariseWZA.py --wza {}_sampled.ZA.csv --LA 0.005 --output {}_sampled.ZA --BayPass" ::: BC_Map cline trunc
