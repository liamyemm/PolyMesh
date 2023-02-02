# $origin is loaded in the runseries.sh script that calls directories.sh, and corresponds to the
# path of the folder from where this runseries.sh will be executed

# Source of all the project
root="/home/liam/github/codes/liamyemm/PolyMesh/2D"

# Location for the schemes' executable
executable=$root"/build/schemes/$executable_name"

# Location of mesh files
meshdir=$root"/typ2_meshes"

# Location for all outputs. 
outdir=$origin"/outputs"

# Names for error file and latex file
errorsfile="data_rates.dat"
latexfile="rates.tex"

