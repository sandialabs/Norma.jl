#!/bin/bash

ncdump holder-0.e >& h0
ncdump holder-1.e >& h1
ncdump gauge.e >& g
sed -i "s/num_elem_var = 49 ;/num_elem_var = 48 ;/g" h0
sed -i "s/num_elem_var = 49 ;/num_elem_var = 48 ;/g" h1
sed -i "s/num_elem_var = 49 ;/num_elem_var = 48 ;/g" g
grep -wn "stored_energy" h0 | cut -d: -f1 >& a1
grep -wn "stored_energy" h1 | cut -d: -f1 >& a2
grep -wn "stored_energy" g | cut -d: -f1 >& a3
X1=$(< a1)
Y1=$(($X1+1))
X2=$(< a2)
Y2=$(($X2+1))
X3=$(< a3)
Y3=$(($X3+1))
sed -i ''"$Y1"'d' h0
sed -i ''"$Y2"'d' h1
sed -i ''"$Y3"'d' g
sed -i 's/"stored_energy",/"stored_energy" ;/g' h0
sed -i 's/"stored_energy",/"stored_energy" ;/g' h1
sed -i 's/"stored_energy",/"stored_energy" ;/g' g
ncgen h0 -o holder-0.e
ncgen h1 -o holder-1.e
ncgen g -o gauge.e
rm a1 a2 a3 h0 h1 g
