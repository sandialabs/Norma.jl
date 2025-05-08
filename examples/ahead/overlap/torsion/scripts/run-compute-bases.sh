#!/bin/bash

#echo "Starting run M=30, reg=1e-1 run..." 
#cd reg1em1
#python make_cubic_opinf_models.py
#echo "...done." 
echo "Starting run M=30, reg=1e-2 run..." 
cd reg1em2
python make_cubic_opinf_models.py
echo "...done." 
echo "Starting run M=30, reg=1e-3 run..." 
cd ../reg1em3
python make_cubic_opinf_models.py
echo "...done." 
echo "Starting run M=30, reg=1e-4 run..." 
cd ../reg1em4
python make_cubic_opinf_models.py
echo "...done." 
echo "Starting run M=30, reg=1e-5 run..." 
cd ../reg1em5
python make_cubic_opinf_models.py
echo "...done." 
echo "Starting run M=30, reg=1e-6 run..." 
cd ../reg1em6
python make_cubic_opinf_models.py
echo "...done." 
echo "Starting run M=30, reg=1e-7 run..." 
cd ../reg1em7
python make_cubic_opinf_models.py
echo "...done." 
echo "Starting run M=30, reg=1e-8 run..." 
cd ../reg1em8
python make_cubic_opinf_models.py
echo "...done." 
echo "Starting run M=30, reg=1e-9 run..." 
cd ../reg1em9
python make_cubic_opinf_models.py
echo "...done." 
echo "Starting run M=30, reg=1e-10 run..." 
cd ../reg1em10
python make_cubic_opinf_models.py
echo "...done."
echo "Starting run M=30, reg=1e-11 run..." 
cd ../reg1em11
python make_cubic_opinf_models.py
echo "...done."
cd ../
