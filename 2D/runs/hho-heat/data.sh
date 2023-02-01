#!/bin/bash

# Plot File
plotfile="plot"

# Order
k=2
l=2

R=0.3;
Q=1.0; 
T_ambient=0.0; 
alpha=0.0; 
k1=1.0; 
k2=1E-3;
k3=1E-2;

#K[1]=0
#K[2]=1
#K[3]=2
#K[4]=3
#K[5]=4
#K[6]=5
#K[7]=6
#K[8]=7
#K[9]=8
#K[10]=9
#K[11]=10
#K[12]=11

#L[1]=0
#L[2]=1
#L[3]=2
#L[4]=3
#L[5]=4
#L[6]=5
#L[7]=6
#L[8]=7
#L[9]=8
#L[10]=9

use_threads="true"
ortho="true"

#mesh[1]=$meshdir"/mesh2_1_transformed_2.typ2"
#mesh[2]=$meshdir"/mesh2_2_transformed_2.typ2"
#mesh[3]=$meshdir"/mesh2_3_transformed_2.typ2"
mesh[1]=$meshdir"/mesh2_6_transformed_10_10.typ2"
mesh[2]=$meshdir"/mesh2_6_transformed_10_10_shift_2.typ2"
mesh[3]=$meshdir"/mesh2_6_transformed_10_10_shift_4.typ2"
mesh[4]=$meshdir"/mesh2_6_transformed_10_10_shift_6.typ2"
mesh[5]=$meshdir"/mesh2_6_transformed_10_10_shift_8.typ2"
#mesh[5]=$meshdir"/mesh2_5_transformed_2.typ2"

#mesh[1]=$meshdir"/mesh2_1.typ2"
#mesh[2]=$meshdir"/mesh2_2.typ2"
#mesh[3]=$meshdir"/mesh2_3.typ2"
#mesh[4]=$meshdir"/mesh2_4.typ2"
#mesh[5]=$meshdir"/mesh2_5.typ2"
#mesh[6]=$meshdir"/mesh2_6.typ2"
#mesh[7]=$meshdir"/mesh2_7.typ2"



