#!/bin/bash

gfortran -fPIC -c ./modules/xraylib.F90 -L/usr/lib64 -lxrlf03a -lxrl -J ./
gfortran -c -fPIC ./modules/types.f90 -o types.o -J./
gfortran -c -fPIC ./modules/xrldata.f90 -o xrldata.o -J./
gfortran -c -fPIC ./modules/cfgdata.f90 -o cfgdata.o -J./
gfortran -c -fPIC ./modules/math.f90 -o math.o -J./
gfortran -c -fPIC ./modules/anode.f90 -o anode.o -J./
gfortran -shared -fPIC types.o xrldata.o cfgdata.o math.o anode.o wrap.f90 -o libanode.so
