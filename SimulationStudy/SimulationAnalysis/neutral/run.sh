# Island

mkdir Island_sampled/
cd Island_sampled/
parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/neutral/Island/{}_0.12Loci.directionalSelection.trees --optima ../../../../slim_configs/BC_Map_environments.14x14.txt --output {}_0.003_1.12Loci.directionalSelection_d40n50_i$i --directional --nPops 40 --nInds 50  --bayPass" ::: $(seq 1 20)
cd ../

exit 0 

#mkdir Island/
#cd Island/
#parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/neutral/Island/{}_0_1.12Loci.directionalSelection.trees --optima ../../../../slim_configs/BC_Map_environments.14x14.txt --output {}_0.003_1.12Loci.directionalSelection_i$i --directional  --bayPass" ::: $(seq 1 20)
#cd ../

for i in 1000 100 10
	do
		mkdir Island_sampled_$i
		cd Island_sampled_$i
		parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/neutral/Island/{}_0_1.12Loci.directionalSelection.trees --optima ../../../../slim_configs/BC_Map_environments.14x14.txt --output {}_0.003_1.12Loci.directionalSelection_d40n50_i$i --directional --nPops 40 --nInds 50 --intervals $i --bayPass" ::: $(seq 1 20)
		cd ../
		
#		mkdir BC_Map_$i/
#		cd BC_Map_$i/
#		parallel -j1 "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/neutral/Island/{}_0_1.12Loci.directionalSelection.trees --optima ../../../../slim_configs/BC_Map_environments.14x14.txt --output {}_0.003_1.12Loci.directionalSelection_i$i --directional --intervals $i --bayPass" ::: $(seq 1 20)
#		cd ../
	done
	

exit 0


# BC Map 

mkdir BC_Map_sampled/
cd BC_Map_sampled/
parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/neutral/2D_SteppingStone/{}_0.003_1.12Loci.directionalSelection.trees --optima ../../../../slim_configs/BC_Map_environments.14x14.txt --output {}_0.003_1.12Loci.directionalSelection_d40n50_i$i --directional --nPops 40 --nInds 50  --bayPass" ::: $(seq 1 20)
cd ../

#mkdir BC_Map/
#cd BC_Map/
#parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/neutral/2D_SteppingStone/{}_0.003_1.12Loci.directionalSelection.trees --optima ../../../../slim_configs/BC_Map_environments.14x14.txt --output {}_0.003_1.12Loci.directionalSelection_i$i --directional  --bayPass" ::: $(seq 1 20)
#cd ../

for i in 1000 100 10
	do
		mkdir BC_Map_sampled_$i
		cd BC_Map_sampled_$i
		parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/neutral/2D_SteppingStone/{}_0.003_1.12Loci.directionalSelection.trees --optima ../../../../slim_configs/BC_Map_environments.14x14.txt --output {}_0.003_1.12Loci.directionalSelection_d40n50_i$i --directional --nPops 40 --nInds 50 --intervals $i --bayPass" ::: $(seq 1 20)
		cd ../
		
#		mkdir BC_Map_$i/
#		cd BC_Map_$i/
#		parallel -j1 "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/neutral/2D_SteppingStone/{}_0.003_1.12Loci.directionalSelection.trees --optima ../../../../slim_configs/BC_Map_environments.14x14.txt --output {}_0.003_1.12Loci.directionalSelection_i$i --directional --intervals $i --bayPass" ::: $(seq 1 20)
#		cd ../
	done
	
