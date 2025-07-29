#!/bin/bash

echo "Starting run M=30, reg=1e-3 run..."
cd M20_reg1em3
julia --project=@. /home/ikalash/Norma.jl-Eric/src/Norma.jl laser-weld.yaml | tee out
echo "...done."
echo "Starting run M=30, reg=1e-4 run..."
cd ../M20_reg1em4
julia --project=@. /home/ikalash/Norma.jl-Eric/src/Norma.jl laser-weld.yaml | tee out
echo "...done."
echo "Starting run M=30, reg=1e-5 run..."
cd ../M20_reg1em5
julia --project=@. /home/ikalash/Norma.jl-Eric/src/Norma.jl laser-weld.yaml | tee out
echo "...done."
echo "Starting run M=30, reg=1e-6 run..."
cd ../M20_reg1em6
julia --project=@. /home/ikalash/Norma.jl-Eric/src/Norma.jl laser-weld.yaml | tee out
echo "...done."
echo "Starting run M=30, reg=1e-7 run..."
cd M20_reg1em7
julia --project=@. /home/ikalash/Norma.jl-Eric/src/Norma.jl laser-weld.yaml | tee out
echo "...done."
echo "Starting run M=30, reg=1e-8 run..."
cd ../M20_reg1em8
julia --project=@. /home/ikalash/Norma.jl-Eric/src/Norma.jl laser-weld.yaml | tee out
echo "...done."
echo "Starting run M=30, reg=1e-9 run..."
cd ../M20_reg1em9
julia --project=@. /home/ikalash/Norma.jl-Eric/src/Norma.jl laser-weld.yaml | tee out
echo "...done."
echo "Starting run M=30, reg=1e-10 run..."
cd ../M20_reg1em10
julia --project=@. /home/ikalash/Norma.jl-Eric/src/Norma.jl laser-weld.yaml | tee out
echo "...done."
cd ../
