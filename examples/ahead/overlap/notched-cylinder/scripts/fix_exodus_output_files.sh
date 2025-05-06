#!/bin/bash

ncdump notched-cylinder-1.e >& n1
ncdump notched-cylinder-2.e >& n2
sed -i "s/num_elem_var = 49 ;/num_elem_var = 48 ;/g" n1
sed -i "s/num_elem_var = 49 ;/num_elem_var = 48 ;/g" n2
grep -wn "stored_energy" n1 | cut -d: -f1 >& a1
grep -wn "stored_energy" n2 | cut -d: -f1 >& a2
X1=$(< a1)
Y1=$(($X1+1))
X2=$(< a2)
Y2=$(($X2+1))
sed -i ''"$Y1"'d' n1
sed -i ''"$Y2"'d' n2
sed -i 's/"stored_energy",/"stored_energy" ;/g' n1
sed -i 's/"stored_energy",/"stored_energy" ;/g' n2
ncgen n1 -o notched-cylinder-1.e
ncgen n2 -o notched-cylinder-2.e
rm a1 a2 n1 n2
