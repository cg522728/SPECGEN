#!/bin/bash
rm *.pyf
rm *.o
rm *.mod
rm *.c

gfortran -fPIC -c ../src/modules/xraylib.F90 -L/usr/lib64 -lxrlf03a -lxrl -J ./

f2py -m integrate -h integrate.pyf ../src/modules/math.f90 only: integrate func :
f2py --f90flags=-lxrlf03 -lxrl -lxrlf03 -lxrl -c integrate.pyf ../src/modules/math.f90 -I../OUTPUT/modules -I../OUTPUT/modules

f2py -m tube_atten -h tube_atten.pyf ../src/modules/anode.f90 only: tube_atten :
f2py --f90flags='-lxrlf03 -lxrl' -lxrlf03 -lxrl -c xraylib.o tube_atten.pyf ../src/modules/anode.f90 ../src/modules/xrldata.f90 ../src/modules/cfgdata.f90 ../src/modules/math.f90 -I../OUTPUT/modules

#f2py -m manode -h manode.pyf ../src/modules/anode.f90
#f2py --f90flags='-lxrlf03 -lxrl' -lxrlf03 -lxrl -c xraylib.o manode.pyf ../src/modules/anode.f90 ../src/modules/xrldata.f90 ../src/modules/cfgdata.f90 ../src/modules/math.f90 ../src/modules/types.f90 -I../OUTPUT/modules

