parallel "python ../../../../../bin/WZA.py --csv ../../../Vs192/{}_sampled/ --bayPass ../{}_sampled/ --output {}_sampled.WZA.csv" ::: BC_Map cline trunc

parallel "python ../../../../../bin/summariseWZA.py --wza {}_sampled.WZA.csv --LA 0.005 --output {}_sampled.WZA --BayPass" ::: BC_Map cline trunc

parallel "python ../../../../../bin/WZA.py --csv ../../../Vs192/BC_Map_sampled_cor{}/ --bayPass ../correlatedEnvironments/BC_Map_sampled_cor{}/ --output BC_Map_sampled_cor{}.WZA.csv" ::: 0.1 0.3 0.5 0.8

parallel "python ../../../../../bin/summariseWZA.py --wza BC_Map_sampled_cor{}.WZA.csv  --LA 0.005 --output BC_Map_sampled_cor{}.WZA --BayPass" ::: 0.1 0.3 0.5 0.8
