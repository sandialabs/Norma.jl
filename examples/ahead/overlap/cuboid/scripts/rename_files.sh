#!/bin/bash
for i in {0..9}; do
  name_old=cuboid-1-disp-000$i.csv
  name_new=01-disp-000$i.csv
  mv $name_old $name_new 
  name_old=cuboid-2-disp-000$i.csv
  name_new=02-disp-000$i.csv
  mv $name_old $name_new 
  name_old=cuboid-1-velo-000$i.csv
  name_new=01-velo-000$i.csv
  mv $name_old $name_new 
  name_old=cuboid-2-velo-000$i.csv
  name_new=02-velo-000$i.csv
  mv $name_old $name_new 
  name_old=cuboid-1-acce-000$i.csv
  name_new=01-acce-000$i.csv
  mv $name_old $name_new 
  name_old=cuboid-2-acce-000$i.csv
  name_new=02-acce-000$i.csv
  mv $name_old $name_new 
  name_old=cuboid-1-time-000$i.csv
  name_new=01-time-000$i.csv
  mv $name_old $name_new 
  name_old=cuboid-2-time-000$i.csv
  name_new=02-time-000$i.csv
  mv $name_old $name_new 
  name_old=cuboid-1-free_dofs-000$i.csv
  name_new=01-free_dofs-000$i.csv
  mv $name_old $name_new 
  name_old=cuboid-2-free_dofs-000$i.csv
  name_new=02-free_dofs-000$i.csv
  mv $name_old $name_new 
  name_old=cuboid-1-nsx--x-velo-000$i.csv
  name_new=01-nsx--x-velo-000$i.csv
  mv $name_old $name_new 
  name_old=cuboid-1-nsx--x-disp-000$i.csv
  name_new=01-nsx--x-disp-000$i.csv
  mv $name_old $name_new 
  name_old=cuboid-1-nsx--x-acce-000$i.csv
  name_new=01-nsx--x-acce-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsy--y-velo-000$i.csv
  name_new=01-nsy--y-velo-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsy--y-curr-000$i.csv
  name_new=01-nsy--y-curr-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsy--y-acce-000$i.csv
  name_new=01-nsy--y-acce-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsz--z-curr-000$i.csv
  name_new=01-nsz--z-curr-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsy--y-disp-000$i.csv
  name_new=01-nsy--y-disp-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsz--z-velo-000$i.csv
  name_new=01-nsz--z-velo-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsz--z-disp-000$i.csv
  name_new=01-nsz--z-disp-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsz--z-acce-000$i.csv
  name_new=01-nsz--z-acce-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsx--x-velo-000$i.csv
  name_new=02-nsx--x-velo-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsx--x-disp-000$i.csv
  name_new=02-nsx--x-disp-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsx--x-acce-000$i.csv
  name_new=02-nsx--x-acce-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsy--y-velo-000$i.csv
  name_new=02-nsy--y-velo-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsy--y-curr-000$i.csv
  name_new=02-nsy--y-curr-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsy--y-acce-000$i.csv
  name_new=02-nsy--y-acce-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsz+-z-velo-000$i.csv
  name_new=02-nsz+-z-velo-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsz+-z-curr-000$i.csv
  name_new=02-nsz+-z-curr-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsy--y-disp-000$i.csv
  name_new=02-nsy--y-disp-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsz+-z-disp-000$i.csv
  name_new=02-nsz+-z-disp-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsz+-z-acce-000$i.csv
  name_old=02-nsz+-z-acce-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-ssz+-curr-000$i.csv
  name_new=01-ssz+-curr-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-ssz+-velo-000$i.csv
  name_new=01-ssz+-velo-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-ssz+-acce-000$i.csv
  name_new=01-ssz+-acce-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-ssz+-disp-000$i.csv
  name_new=01-ssz+-disp-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-ssz--curr-000$i.csv
  name_new=02-ssz--curr-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-ssz--velo-000$i.csv
  name_new=02-ssz--velo-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-ssz--disp-000$i.csv
  name_new=02-ssz--disp-000$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-ssz--acce-000$i.csv
  name_new=02-ssz--acce-000$i.csv
  mv $name_old $name_new
done
for i in {10..99}; do
  name_old=cuboid-1-disp-00$i.csv
  name_new=01-disp-00$i.csv
  mv $name_old $name_new 
  name_old=cuboid-2-disp-00$i.csv
  name_new=02-disp-00$i.csv
  mv $name_old $name_new 
  name_old=cuboid-1-velo-00$i.csv
  name_new=01-velo-00$i.csv
  mv $name_old $name_new 
  name_old=cuboid-2-velo-00$i.csv
  name_new=02-velo-00$i.csv
  mv $name_old $name_new 
  name_old=cuboid-1-acce-00$i.csv
  name_new=01-acce-00$i.csv
  mv $name_old $name_new 
  name_old=cuboid-2-acce-00$i.csv
  name_new=02-acce-00$i.csv
  mv $name_old $name_new 
  name_old=cuboid-1-time-00$i.csv
  name_new=01-time-00$i.csv
  mv $name_old $name_new 
  name_old=cuboid-2-time-00$i.csv
  name_new=02-time-00$i.csv
  mv $name_old $name_new 
  name_old=cuboid-1-free_dofs-00$i.csv
  name_new=01-free_dofs-00$i.csv
  mv $name_old $name_new 
  name_old=cuboid-2-free_dofs-00$i.csv
  name_new=02-free_dofs-00$i.csv
  mv $name_old $name_new 
  name_old=cuboid-1-nsx--x-velo-00$i.csv
  name_new=01-nsx--x-velo-00$i.csv
  mv $name_old $name_new 
  name_old=cuboid-1-nsx--x-disp-00$i.csv
  name_new=01-nsx--x-disp-00$i.csv
  mv $name_old $name_new 
  name_old=cuboid-1-nsx--x-acce-00$i.csv
  name_new=01-nsx--x-acce-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsy--y-velo-00$i.csv
  name_new=01-nsy--y-velo-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsy--y-curr-00$i.csv
  name_new=01-nsy--y-curr-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsy--y-acce-00$i.csv
  name_new=01-nsy--y-acce-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsz--z-curr-00$i.csv
  name_new=01-nsz--z-curr-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsy--y-disp-00$i.csv
  name_new=01-nsy--y-disp-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsz--z-velo-00$i.csv
  name_new=01-nsz--z-velo-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsz--z-disp-00$i.csv
  name_new=01-nsz--z-disp-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsz--z-acce-00$i.csv
  name_new=01-nsz--z-acce-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsx--x-velo-00$i.csv
  name_new=02-nsx--x-velo-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsx--x-disp-00$i.csv
  name_new=02-nsx--x-disp-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsx--x-acce-00$i.csv
  name_new=02-nsx--x-acce-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsy--y-velo-00$i.csv
  name_new=02-nsy--y-velo-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsy--y-curr-00$i.csv
  name_new=02-nsy--y-curr-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsy--y-acce-00$i.csv
  name_new=02-nsy--y-acce-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsz+-z-velo-00$i.csv
  name_new=02-nsz+-z-velo-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsz+-z-curr-00$i.csv
  name_new=02-nsz+-z-curr-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsy--y-disp-00$i.csv
  name_new=02-nsy--y-disp-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsz+-z-disp-00$i.csv
  name_new=02-nsz+-z-disp-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsz+-z-acce-00$i.csv
  name_old=02-nsz+-z-acce-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-ssz+-curr-00$i.csv
  name_new=01-ssz+-curr-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-ssz+-velo-00$i.csv
  name_new=01-ssz+-velo-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-ssz+-acce-00$i.csv
  name_new=01-ssz+-acce-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-ssz+-disp-00$i.csv
  name_new=01-ssz+-disp-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-ssz--curr-00$i.csv
  name_new=02-ssz--curr-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-ssz--velo-00$i.csv
  name_new=02-ssz--velo-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-ssz--disp-00$i.csv
  name_new=02-ssz--disp-00$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-ssz--acce-00$i.csv
  name_new=02-ssz--acce-00$i.csv
  mv $name_old $name_new
done
for i in {100..999}; do
  name_old=cuboid-1-disp-0$i.csv
  name_new=01-disp-0$i.csv
  mv $name_old $name_new 
  name_old=cuboid-2-disp-0$i.csv
  name_new=02-disp-0$i.csv
  mv $name_old $name_new 
  name_old=cuboid-1-velo-0$i.csv
  name_new=01-velo-0$i.csv
  mv $name_old $name_new 
  name_old=cuboid-2-velo-0$i.csv
  name_new=02-velo-0$i.csv
  mv $name_old $name_new 
  name_old=cuboid-1-acce-0$i.csv
  name_new=01-acce-0$i.csv
  mv $name_old $name_new 
  name_old=cuboid-2-acce-0$i.csv
  name_new=02-acce-0$i.csv
  mv $name_old $name_new 
  name_old=cuboid-1-time-0$i.csv
  name_new=01-time-0$i.csv
  mv $name_old $name_new 
  name_old=cuboid-2-time-0$i.csv
  name_new=02-time-0$i.csv
  mv $name_old $name_new 
  name_old=cuboid-1-free_dofs-0$i.csv
  name_new=01-free_dofs-0$i.csv
  mv $name_old $name_new 
  name_old=cuboid-2-free_dofs-0$i.csv
  name_new=02-free_dofs-0$i.csv
  mv $name_old $name_new 
  name_old=cuboid-1-nsx--x-velo-0$i.csv
  name_new=01-nsx--x-velo-0$i.csv
  mv $name_old $name_new 
  name_old=cuboid-1-nsx--x-disp-0$i.csv
  name_new=01-nsx--x-disp-0$i.csv
  mv $name_old $name_new 
  name_old=cuboid-1-nsx--x-acce-0$i.csv
  name_new=01-nsx--x-acce-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsy--y-velo-0$i.csv
  name_new=01-nsy--y-velo-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsy--y-curr-0$i.csv
  name_new=01-nsy--y-curr-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsy--y-acce-0$i.csv
  name_new=01-nsy--y-acce-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsz--z-curr-0$i.csv
  name_new=01-nsz--z-curr-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsy--y-disp-0$i.csv
  name_new=01-nsy--y-disp-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsz--z-velo-0$i.csv
  name_new=01-nsz--z-velo-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsz--z-disp-0$i.csv
  name_new=01-nsz--z-disp-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-nsz--z-acce-0$i.csv
  name_new=01-nsz--z-acce-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsx--x-velo-0$i.csv
  name_new=02-nsx--x-velo-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsx--x-disp-0$i.csv
  name_new=02-nsx--x-disp-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsx--x-acce-0$i.csv
  name_new=02-nsx--x-acce-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsy--y-velo-0$i.csv
  name_new=02-nsy--y-velo-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsy--y-curr-0$i.csv
  name_new=02-nsy--y-curr-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsy--y-acce-0$i.csv
  name_new=02-nsy--y-acce-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsz+-z-velo-0$i.csv
  name_new=02-nsz+-z-velo-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsz+-z-curr-0$i.csv
  name_new=02-nsz+-z-curr-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsy--y-disp-0$i.csv
  name_new=02-nsy--y-disp-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsz+-z-disp-0$i.csv
  name_new=02-nsz+-z-disp-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-nsz+-z-acce-0$i.csv
  name_old=02-nsz+-z-acce-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-ssz+-curr-0$i.csv
  name_new=01-ssz+-curr-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-ssz+-velo-0$i.csv
  name_new=01-ssz+-velo-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-ssz+-acce-0$i.csv
  name_new=01-ssz+-acce-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-1-ssz+-disp-0$i.csv
  name_new=01-ssz+-disp-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-ssz--curr-0$i.csv
  name_new=02-ssz--curr-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-ssz--velo-0$i.csv
  name_new=02-ssz--velo-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-ssz--disp-0$i.csv
  name_new=02-ssz--disp-0$i.csv
  mv $name_old $name_new
  name_old=cuboid-2-ssz--acce-0$i.csv
  name_new=02-ssz--acce-0$i.csv
  mv $name_old $name_new
done
