#/bin/bash


grep "Performed" out >& schwarz-iters.txt
sed -i "s/Performed//g" schwarz-iters.txt
sed -i "s/Schwarz Iterations//g" schwarz-iters.txt
awk '{print substr($0, 6)}' schwarz-iters.txt >& a
mv a schwarz-iters.txt

