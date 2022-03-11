

for d in 10 20
	do
	mkdir cline_sampled_d${d}
	cd cline_sampled_d${d}
	parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/s0.003/cline/{}_0.003_3.12Loci.directionalSelection.trees --optima ../../../../slim_configs/oneDCline_environments.14x14.txt --output {}_0.003_directional_i10000   --nPops $d --nInds 50 --intervals 10000 --directional" ::: $(seq 1 20)
	cd ../

	mkdir trunc_sampled_d${d}
	cd trunc_sampled_d${d}
	parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/s0.003/trunc/{}_0.003_2.12Loci.directionalSelection.trees --optima ../../../../slim_configs/trunc_environments.14x14.txt --output {}_0.003_directional_i10000   --nPops $d --nInds 50 --intervals 10000 --directional" ::: $(seq 1 20)
	cd ../
	
	mkdir BC_Map_sampled_d${d}
	cd BC_Map_sampled_d${d}
	parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/s0.003/BC_Map/{}_0.003_1.12Loci.directionalSelection.trees --optima ../../../../slim_configs/BC_Map_environments.14x14.txt --output {}_0.003_directional_i10000   --nPops $d --nInds 50 --intervals 10000 --directional" ::: $(seq 1 20)
	cd ../

	done



for n in 5 10 20
	do
	for d in 10 20 40
		do
		mkdir BC_Map_sampled_d${d}_n${n}
		cd BC_Map_sampled_d${d}_n${n}
		parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/GEA/G.2.3/s0.003/BC_Map/{}_0.003_1.12Loci.directionalSelection.trees --optima ../../../../slim_configs/BC_Map_environments.14x14.txt --output {}_0.003_directional_i10000   --nPops $d --nInds $n --intervals 10000 --directional" ::: $(seq 1 20)
		cd ../
		python ../../../bin/WZA.py --csv BC_Map_sampled_d${d}_n${n}/ --output BC_Map_sampled_d${d}_n${n}.WZA.csv 
		python ../../../bin/summariseWZA.py --wza BC_Map_sampled_d${d}_n${n}.WZA.csv --LA 0.005 --output BC_Map_sampled_d${d}_n${n}.WZA
		done
	done

exit 0



parallel "python ../../../bin/WZA.py --csv {1}_sampled_d{2}/ --output {1}_sampled_d{2}.WZA.csv " ::: BC_Map cline trunc ::: 10 20

parallel "python ../../../bin/summariseWZA.py --wza {1}_sampled_d{2}.WZA.csv --LA 0.005 --output {1}_sampled_d{2}.WZA.csv" ::: BC_Map trunc cline ::: 10 20


parallel "python ../../../bin/WZA.py --csv {1}/ --output {1}.WZA.csv " ::: BC_Map cline trunc

parallel "python ../../../bin/summariseWZA.py --wza {1}.WZA.csv --LA 0.005 --output {1}.WZA.csv" ::: BC_Map trunc cline 

parallel "python ../../../bin/WZA.py --csv {1}_sampled/ --output {1}_sampled.WZA.csv " :::  BC_Map cline trunc

parallel "python ../../../bin/summariseWZA.py --wza {1}_sampled.WZA.csv --LA 0.005 --output {1}_sampled.WZA" :::  BC_Map cline trunc


exit 0

for n in 5 10 20
	do
	for d in 10 20 40
		do
		mkdir BC_Map_sampled_d${d}_n${n}
		cd BC_Map_sampled_d${d}_n${n}
		parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/GEA/G.2.3/s0.003/BC_Map/{}_0.003_1.12Loci.directionalSelection.trees --optima ../../../../slim_configs/BC_Map_environments.14x14.txt --output {}_0.003_directional_i10000   --nPops $d --nInds $n --intervals 10000 --directional" ::: $(seq 1 20)
		cd ../
		
		done
	done


exit 0
#################
### TRUNC
#################
i=10000

mkdir trunc_sampled/
cd trunc_sampled/
parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/s0.003/trunc/{}_0.003_2.12Loci.directionalSelection.trees --optima ../../../../slim_configs/trunc_environments.14x14.txt --output {}_0.003_directional_d40n50_i$i   --nPops 40 --nInds 50  --bayPass --directional" ::: $(seq 1 20)
cd ../

mkdir trunc/
cd trunc/
parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/s0.003/trunc/{}_0.003_2.12Loci.directionalSelection.trees --optima ../../../../slim_configs/trunc_environments.14x14.txt --output {}_0.003_directional_i$i  --bayPass --directional" ::: $(seq 1 20)
cd ../

#for i in 1000 100 10
#	do
#		mkdir trunc_sampled_$i
#		cd trunc_sampled_$i
#		parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/s0.003/trunc/{}_0.003_2.12Loci.directionalSelection.trees --optima ../../../../slim_configs/trunc_environments.14x14.txt --output {}_0.003_directional_d40n50_i$i  --intervals $i  --nPops 40 --nInds 50 --directional" ::: $(seq 1 20)
#		cd ../
		
#	done




	
#################
### CLINE
#################

i=10000

mkdir cline_sampled/
cd cline_sampled/
parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/s0.003/cline/{}_0.003_3.12Loci.directionalSelection.trees --optima ../../../../slim_configs/oneDCline_environments.14x14.txt --output {}_0.003_directional_d40n50_i$i   --nPops 40 --nInds 50  --bayPass --directional" ::: $(seq 1 20)
cd ../

mkdir cline/
cd cline/
parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/GEA/G.2.3/s0.003/cline/{}_0.003_3.12Loci.directionalSelection.trees --optima ../../../../slim_configs/oneDCline_environments.14x14.txt --output {}_0.003_directional_i$i  --directional" ::: $(seq 1 20)
cd ../

#for i in 1000 100 10
#	do
#		mkdir cline_sampled_$i
#		cd cline_sampled_$i
#		parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/GEA/G.2.3/s0.003/cline/{}_0.003_3.12Loci.directionalSelection.trees --optima ../../../../slim_configs/oneDCline_environments.14x14.txt --output {}_0.003_directional_d40n50_i$i  --intervals $i  --nPops 40 --nInds 50 --directional" ::: $(seq 1 20)
#		cd ../
		
#	done
	
	
	
#################
### BC_Map
#################

i=10000

mkdir BC_Map_sampled/

cd BC_Map_sampled/
parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/s0.003/BC_Map/{}_0.003_1.12Loci.directionalSelection.trees --optima ../../../../slim_configs/BC_Map_environments.14x14.txt --output {}_0.003_directional_d40n50_i$i   --nPops 40 --nInds 50  --bayPass --directional" ::: $(seq 1 20)
cd ../

mkdir BC_Map/
cd BC_Map/
parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/GEA/G.2.3/s0.003/BC_Map/{}_0.003_1.12Loci.directionalSelection.trees --optima ../../../../slim_configs/BC_Map_environments.14x14.txt --output {}_0.003_directional_i$i  --directional" ::: $(seq 1 20)
cd ../

for i in 1000 100 10
	do
		mkdir BC_Map_sampled_$i
		cd BC_Map_sampled_$i
		parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/s0.003/BC_Map/{}_0.003_1.12Loci.directionalSelection.trees --optima ../../../../slim_configs/BC_Map_environments.14x14.txt --output {}_0.003_directional_d40n50_i$i   --nPops 40 --nInds 50 --intervals $i --directional" ::: $(seq 1 20)
		cd ../
		
	done
	


