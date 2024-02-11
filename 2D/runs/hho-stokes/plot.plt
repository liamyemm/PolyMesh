set size 1,1
set terminal pdf linewidth 0.05 size 30,30
set output "./mesh_plot.pdf"
set xrange [-1:1]
set yrange [-1:1]
unset tics
unset title
unset border
plot "/home/liam/github/codes/liamyemm/PolyMesh/2D/build/schemes/XVEM/mesh_plot.dat" with lines notitle lt rgb "black", 


