#!/bin/bash
#
# Execute hho file on series of meshes, and calculate outputs
#

###
# Executable file for the considered scheme
executable_name="HHO_Stokes/hho_stokes"

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

nbradius=${#radius[@]}
nbmesh=${#mesh[@]}
nbK=${#K[@]}
for i in `seq 1 $nbmesh`; 
do
for j in `seq 1 $nbradius`; 
do
for k in `seq 1 $nbK`; 
do
  rad=${radius[$j]}
  meshfile=${mesh[$i]}
#  meshfile=${mesh[1]}
#  k=${K[$i]}
#  l=$k
#  l=${L[$j]}
#  echo -e "\n*************************\nmesh $i out of $nbmesh: ${mesh[$i]}"
#	 Execute code
  	$executable --mesh $meshfile --edgedegree ${K[$k]} --celldegree ${K[$k]} --orthonormalise $ortho --use_threads $use_threads --plot $plotfile --radius $rad
#	 Move outputs
	mv results.txt $outdir/results-$i-$j-$k.txt
done;
done;
done;

# CREATE DATA FILE FOR LATEX
for j in `seq 1 $nbradius`; 
do
for k in `seq 1 $nbK`; 
do
echo -e "MeshSize L2Error H1Error EnergyError PressureError NbCells NbEdges NbInternalEdges DOFs EdgeDegree Radius" > "outputs/data_rates.dat"
for i in `seq 1 $nbmesh`; 
do	
	MeshSize=$(awk '/MeshSize: / {print $NF}' $outdir/results-$i-$j-$k.txt)
	L2Error=$(awk '/L2Error: / {print $NF}' $outdir/results-$i-$j-$k.txt)
	H1Error=$(awk '/H1Error: / {print $NF}' $outdir/results-$i-$j-$k.txt)
	EnergyError=$(awk '/EnergyError: / {print $NF}' $outdir/results-$i-$j-$k.txt)
	PressureError=$(awk '/PressureError: / {print $NF}' $outdir/results-$i-$j-$k.txt)
	NbCells=$(awk '/NbCells: / {print $NF}' $outdir/results-$i-$j-$k.txt)
	NbEdges=$(awk '/NbEdges: / {print $NF}' $outdir/results-$i-$j-$k.txt)
	NbInternalEdges=$(awk '/NbInternalEdges: / {print $NF}' $outdir/results-$i-$j-$k.txt)
	DOFs=$(awk '/DOFs: / {print $NF}' $outdir/results-$i-$j-$k.txt)
	EdgeDegree=$(awk '/EdgeDegree: / {print $NF}' $outdir/results-$i-$j-$k.txt)
	Radius=$(awk '/Radius: / {print $NF}' $outdir/results-$i-$j-$k.txt)
echo -e "$MeshSize $L2Error $H1Error $EnergyError $PressureError $NbCells $NbEdges $NbInternalEdges $DOFs $EdgeDegree $Radius" >> "outputs/data_rates.dat"
done;
#mv outputs/data_rates.dat ~/github/PhD/Papers/XHHO_Stokes/data/nonenriched_4circles_k${K[$k]}.dat
done;
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
  g++ compute_rates.cpp -o compute_rates
#fi
./compute_rates 11 
#| tee -a $outdir/allrates.txt
