parallel "python ../../../../../bin/WZA.py --csv ../../../s0.003/{}_sampled/ --bayPass ../{}_sampled_BayPass/ --output {}_sampled.WZA.csv" ::: BC_Map cline trunc

parallel "python ../../../../../bin/summariseWZA.py --wza {}_sampled.WZA.csv --LA 0.005 --output {}_sampled.WZA --BayPass" ::: BC_Map cline trunc
