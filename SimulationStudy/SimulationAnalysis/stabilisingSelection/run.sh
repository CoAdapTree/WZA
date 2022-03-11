
#################
### BC_Map
#################

i=10000

mkdir BC_Map_sampled/

cd BC_Map_sampled/
parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/Vs192/BC_Map/{}_0.5_192.trees --optima ../../../../slim_configs/BC_Map_environments.14x14.txt --output {}_0.5_192_d40n50_i$i  --nPops 40 --nInds 50  --bayPass" ::: $(seq 1 20)
cd ../

mkdir BC_Map/
cd BC_Map/
parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/Vs192/BC_Map/{}_0.5_192.trees --optima ../../../../slim_configs/BC_Map_environments.14x14.txt --output {}_0.5_192_d40n50_i$i   --bayPass" ::: $(seq 1 20)
cd ../

for i in 1000 100 10
	do
		mkdir BC_Map_sampled_$i
		cd BC_Map_sampled_$i
		parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/Vs192/BC_Map/{}_0.5_192.trees --optima ../../../../slim_configs/BC_Map_environments.14x14.txt --output {}_0.5_192_d40n50_i$i  --nPops 40 --nInds 50 --intervals $i --bayPass" ::: $(seq 1 20)
		cd ../
		
#		mkdir BC_Map_$i/
#		cd BC_Map_$i/
#		parallel -j1 "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/Vs192/BC_Map/{}_0.5_192.trees --optima ../../../../slim_configs/BC_Map_environments.14x14.txt --output {}_0.5_192_i$i  --intervals $i --bayPass" ::: $(seq 1 20)
#		cd ../
	done
	
#################
### CLINE
#################

i=10000

mkdir cline_sampled/
cd cline_sampled/
parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/Vs192/cline/{}_0.5_192.trees --optima ../../../../slim_configs/oneDCline_environments.14x14.txt --output {}__0.5_192_d40n50_i$i  --nPops 40 --nInds 50  --bayPass" ::: $(seq 1 20)
cd ../

mkdir cline/
cd cline/
parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/Vs192/cline/{}_0.5_192.trees --optima ../../../../slim_configs/oneDCline_environments.14x14.txt --output {}_0.5_192_i$i   --bayPass" ::: $(seq 1 20)
cd ../

for i in 1000 100 10
	do
		mkdir cline_sampled_$i
		cd cline_sampled_$i
		parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/Vs192/cline/{}_0.5_192.trees --optima ../../../../slim_configs/oneDCline_environments.14x14.txt --output {}_0.5_192_d40n50_i$i  --nPops 40 --nInds 50 --intervals $i --bayPass" ::: $(seq 1 20)
		cd ../
		
#		mkdir cline_$i/
#		cd cline_$i/
#		parallel -j1 "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/Vs192/cline/{}_0.5_192.trees --optima ../../../../slim_configs/oneDCline_environments.14x14.txt --output {}_0.5_192_i$i  --intervals $i --bayPass" ::: $(seq 1 20)
#		cd ../
	done
	
	
	

#################
### TRUNC
#################

i=10000


mkdir trunc_sampled/
cd trunc_sampled/
parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/Vs192/trunc/{}_0.5_192.trees --optima ../../../../slim_configs/trunc_environments.14x14.txt --output {}__0.5_192_d40n50_i$i  --nPops 40 --nInds 50  --bayPass" ::: $(seq 1 20)
cd ../

mkdir trunc/
cd trunc/
parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/Vs192/trunc/{}_0.5_192.trees --optima ../../../../slim_configs/trunc_environments.14x14.txt --output {}_0.5_192_i$i   --bayPass" ::: $(seq 1 20)
cd ../

for i in 1000 100 10
	do
		mkdir trunc_sampled_$i
		cd trunc_sampled_$i
		parallel "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/Vs192/trunc/{}_0.5_192.trees --optima ../../../../slim_configs/trunc_environments.14x14.txt --output {}_0.5_192_d40n50_i$i  --nPops 40 --nInds 50 --intervals $i --bayPass" ::: $(seq 1 20)
		cd ../
		
#		mkdir trunc_$i/
#		cd trunc_$i/
#		parallel -j1 "python ../../../../bin/parseSimulations.py --trees /media/booker/HOWDY/G.2.3/Vs192/trunc/{}_0.5_192.trees --optima ../../../../slim_configs/trunc_environments.14x14.txt --output {}_0.5_192_i$i  --intervals $i --bayPass" ::: $(seq 1 20)
#		cd ../
	done
