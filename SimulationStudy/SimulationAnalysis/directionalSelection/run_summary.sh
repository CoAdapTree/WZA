


parallel "python ../../../bin/WZA.py --csv {1}_sampled/ --output {1}_sampled.WZA.csv" ::: BC_Map cline trunc
parallel "python ../../../bin/summariseWZA.py --wza {1}_sampled.WZA.csv --output {1}_sampled.WZA --LA 0.005" ::: BC_Map cline trunc

parallel "python ../../../bin/WZA.py --csv {1}_sampled_d{2}/ --output {1}_sampled_d{2}.WZA.csv" ::: BC_Map cline trunc ::: 10 20
parallel "python ../../../bin/summariseWZA.py --wza {1}_sampled_d{2}.WZA.csv --output {1}_sampled_d{2}.WZA --LA 0.005" ::: BC_Map cline trunc ::: 10 20


parallel "python ../../../bin/WZA.py --csv {1}/ --output {1}.WZA.csv" ::: BC_Map cline trunc
parallel "python ../../../bin/summariseWZA.py --wza {1}.WZA.csv --output {1}.WZA --LA 0.005" ::: BC_Map cline trunc
