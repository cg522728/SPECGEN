#!/bin/bash
f2py -m mathchg -h mathchg.pyf ../src/modules/math.f90
f2py -c mathchg.pyf ../src/modules/math.f90 -I../OUTPUT/modules
