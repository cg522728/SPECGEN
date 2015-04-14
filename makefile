odir = modules
srcdir = src
_obj = \
	xraylib.o\
	cfgdata.o\
	xrldata.o\
	constants.o\
	anode.o\
	pella.o\
	seccomp.o\
	micromatter.o
obj = $(patsubst %,$(odir)/%,$(_obj))
all:
	gfortran -fopenmp -funderscoring -cpp -Wall -L/usr/lib64 -c src/modules/xraylib.f90 -o modules/xraylib.o -lxrlf03 -lxrl -lgomp
	gfortran -fopenmp -funderscoring -cpp -Wall -L/usr/lib64 -c src/modules/cfgdata.f90 -o modules/cfgdata.o -lxrlf03 -lxrl -lgomp
	gfortran -fopenmp -funderscoring -cpp -Wall -L/usr/lib64 -c src/modules/xrldata.f90 -o modules/xrldata.o -lxrlf03 -lxrl -lgomp
	gfortran -fopenmp -funderscoring -cpp -Wall -L/usr/lib64 -c src/modules/constants.f90 -o modules/constants.o -lxrlf03 -lxrl -lgomp
	gfortran -fopenmp -funderscoring -cpp -Wall -L/usr/lib64 -c src/modules/anode.f90 -o modules/anode.o -lxrlf03 -lxrl -lgomp
	gfortran -fopenmp -funderscoring -cpp -Wall -L/usr/lib64 -c src/modules/pella.f90 -o modules/pella.o -lxrlf03 -lxrl -lgomp
	gfortran -fopenmp -funderscoring -cpp -Wall -L/usr/lib64 -c src/modules/seccomp.f90 -o modules/seccomp.o -lxrlf03 -lxrl -lgomp
	gfortran -fopenmp -funderscoring -cpp -Wall -L/usr/lib64 -c src/modules/micromatter.f90 -o modules/micromatter.o -lxrlf03 -lxrl -lgomp
	gfortran -fopenmp -funderscoring -cpp -Wall -L/usr/lib64 -o micromatter src/micromatter.f90 $(obj) -lxrlf03 -lxrl -lgomp
	gfortran -fopenmp -funderscoring -cpp -Wall -L/usr/lib64 -o tspec src/tspec.f90 $(obj) -lxrlf03 -lxrl -lgomp
	gfortran -fopenmp -funderscoring -cpp -Wall -L/usr/lib64 -o stspec src/stspec.f90 $(obj) -lxrlf03 -lxrl -lgomp
clean:
	rm -fr modules/*.o
	rm -fr *.mod