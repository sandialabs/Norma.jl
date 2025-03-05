#/bin/bash

grep "Performed" out >& schwarz-iters.txt 
sed -i "s/Performed//g" schwarz-iters.txt
sed -i "s/Schwarz iterations//g" schwarz-iters.txt
