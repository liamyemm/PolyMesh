#!/bin/bash
#
# Execute hho file on series of meshes, and calculate outputs
#

###
# Executable file for the considered scheme
executable_name="HHO_Poisson/hho_poisson"

###
# Directories
origin=$(pwd)
if [ ! -f ../directories.sh ]; then
  echo "directories.sh does not exist. Please read the README.txt in the parent folder."
  exit
fi
. ../directories.sh

# Options:
if [[ $1 == "help" ]]; then
	echo -e "\nExecute tests using parameters in data.sh, creates and compile latex file, and calculate rates.\n
Executed without parameters: uses the data in the local data.sh file\n"
	exit;
fi;

# Create/clean output directory
if [ -d $outdir ]; then
	\rm -r $outdir
fi
	mkdir $outdir

###
# LOAD DATA
datafile=$(readlink -f data.sh);
echo -e "Using data file $datafile\n"
. $datafile

#nbmesh=${#mesh[@]}
nbK=${#K[@]}
#for i in `seq 1 $nbmesh`; 
for i in `seq 1 $nbK`; 
do
#  meshfile=${mesh[$i]}
  meshfile=${mesh[1]}
  k=${K[$i]}
  l=$k
#  l=${L[$j]}
#  echo -e "\n*************************\nmesh $i out of $nbmesh: ${mesh[$i]}"
#	 Execute code
  	$executable --mesh $meshfile  --testcase $sol --boundary $bdry --edgedegree $k --celldegree $l --orthonormalise $ortho --use_threads $use_threads --plot plot
#	 Move outputs
	mv results.txt $outdir/results-$i.txt
#done;
done;

# CREATE DATA FILE FOR LATEX
echo -e "MeshSize L2Error H1Error EnergyError NbCells NbEdges NbInternalEdges EdgeDegree" > "outputs/data_rates.dat"
#for i in `seq 1 $nbmesh`; 
for i in `seq 1 $nbK`; 
do	
	MeshSize=$(awk '/MeshSize: / {print $NF}' $outdir/results-$i.txt)
	L2Error=$(awk '/L2Error: / {print $NF}' $outdir/results-$i.txt)
	H1Error=$(awk '/H1Error: / {print $NF}' $outdir/results-$i.txt)
	EnergyError=$(awk '/EnergyError: / {print $NF}' $outdir/results-$i.txt)
	NbCells=$(awk '/NbCells: / {print $NF}' $outdir/results-$i.txt)
	NbEdges=$(awk '/NbEdges: / {print $NF}' $outdir/results-$i.txt)
	NbInternalEdges=$(awk '/NbInternalEdges: / {print $NF}' $outdir/results-$i.txt)
#	DOFs=$(awk '/DOFs: / {print $NF}' $outdir/results-$i-$j.txt)
	EdgeDegree=$(awk '/EdgeDegree: / {print $NF}' $outdir/results-$i.txt)
echo -e "$MeshSize $L2Error $H1Error $EnergyError $NbCells $NbEdges $NbInternalEdges $EdgeDegree" >> "outputs/data_rates.dat"
done;

###
## COMPUTATION OF CONVERGENCE RATES
##
echo -e "\n ----------- Data -----------"
echo "degrees: (edge) k = $k, (cell) l = $l"
#cd ..
#echo "data file: $datafile" >> $outdir/allrates.txt
#echo "Mesh 1: ${mesh[1]}" >> $outdir/allrates.txt
#if [ ! -f ./compute_rates ]; then
#  g++ compute_rates.cpp -o compute_rates
#fi
./compute_rates 8 | tee -a $outdir/allrates.txt
