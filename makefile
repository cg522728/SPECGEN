CC = gfortran
CDFLAGS = -funderscoring -O0 -g3 -pg -p -Wall -fopenmp -fmessage-length=0 -cpp
COFLAGS = -funderscoring -O3 -Wall -fopenmp -fmessage-length=0 -cpp
CFLAGS = $(COFLAGS)
LIBS = -lxrlf03 -lxrl -lgomp
_SOURCES =\
	types\
	cfgdata\
	xrldata\
	math\
	anode\
	seccomp\
	micromatter
_MAINS = \
	tspec\
	stspec\
	mmsens\
	deteff
_OBJECTS = $(_SOURCES:%=%.o)
_EXEC = $(_MAINS)
SDIR = src
WDIR = OUTPUT
MODDIR = modules
SOURCES = $(patsubst %,$(MODDIR:%=$(SDIR)/%)/%,$(_SOURCES:%=%.f90))
OBJECTS = $(patsubst %,$(MODDIR:%=$(WDIR)/%)/%,$(_OBJECTS))
MAINS = $(patsubst %,$(SDIR)/%,$(_MAINS:%=%.f90))
EXEC = $(patsubst %,$(WDIR)/%,$(_EXEC))

all: executables

.PHONY: all

.PHONY: modules $(_SOURCES)

.PHONY: executables $(_MAINS)

executables :	$(_MAINS)

$(_MAINS): modules
	$(CC) $(CFLAGS) -o $(patsubst %,$(WDIR)/%,$@) $(patsubst %,$(SDIR)/%,$@).f90 $(OBJECTS) -I $(MODDIR:%=$(WDIR)/%) $(LIBS)

modules : $(_SOURCES)

$(_SOURCES): xraylib.o
		$(CC) $(CFLAGS) -c $(patsubst %,$(MODDIR:%=$(SDIR)/%)/%,$@).f90 -o $(patsubst %,$(MODDIR:%=$(WDIR)/%)/%,$@).o -J $(MODDIR:%=$(WDIR)/%) $(LIBS)
		
xraylib.o:
	$(CC) -funderscoring -O0 -g3 -pg -p -fopenmp -fmessage-length=0 -cpp -L/usr/lib64 -c $(SDIR)/$(MODDIR)/xraylib.F90 -o $(WDIR)/$(MODDIR)/xraylib.o -J $(MODDIR:%=$(WDIR)/%) $(LIBS)
clean:
	rm -fr modules/*.o
	rm -fr *.mod
	rm -fr mmsens deteff
