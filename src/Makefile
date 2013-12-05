# SPARTA multiple-machine Makefile

SHELL = /bin/bash

# Definitions

ROOT =	spa
EXE =	$(ROOT)_$@
SRC =	$(wildcard *.cpp)
INC =	$(wildcard *.h)
OBJ = 	$(SRC:.cpp=.o)

# List of all targets

help:
	@echo ''
	@echo 'make clean-all           delete all object files'
	@echo 'make clean-machine       delete object files for one machine'
	@echo 'make tar                 lmp_src.tar.gz of src dir and packages'
	@echo 'make makelib             create Makefile.lib for static library build'
	@echo 'make makeshlib           create Makefile.shlib for shared library build'
	@echo 'make makelist            create Makefile.list used by old makes'
	@echo 'make -f Makefile.lib machine      build SPARTA as static library for machine'
	@echo 'make -f Makefile.shlib machine    build SPARTA as shared library for machine'
	@echo 'make -f Makefile.list machine     build SPARTA from explicit list of files'
	@echo 'make stubs               build dummy MPI library in STUBS'
	@echo 'make install-python      install SPARTA wrapper in Python'
	@echo ''
	@echo 'make machine             build SPARTA where machine is one of:'
	@echo ''
	@files="`ls MAKE/Makefile.*`"; \
	  for file in $$files; do head -1 $$file; done
	@echo ''

# Build the code

.DEFAULT:
	@test -f MAKE/Makefile.$@
	@if [ ! -d Obj_$@ ]; then mkdir Obj_$@; fi
	@$(SHELL) Make.sh style
	@cp MAKE/Makefile.$@ Obj_$@/Makefile
	@cd Obj_$@; \
	$(MAKE) $(MFLAGS) "OBJ = $(OBJ)" "INC = $(INC)" "SHFLAGS =" \
	  "EXE = ../$(EXE)" ../$(EXE)

# Remove machine-specific object files

clean:
	@echo 'make clean-all           delete all object files'
	@echo 'make clean-machine       delete object files for one machine'

clean-all:
	rm -rf Obj_*

clean-%:
	rm -rf Obj_$(@:clean-%=%)

# Create a tarball of src dir and packages

tar:
	@cd STUBS; make clean
	@cd ..; tar cvzf src/$(ROOT)_src.tar.gz \
	  src/Make* src/MAKE src/*.cpp src/*.h src/STUBS \
          --exclude=*/.svn
	@cd STUBS; make
	@echo "Created $(ROOT)_src.tar.gz"

# Make MPI STUBS lib

stubs:
	@cd STUBS; make clean; make

# Create Makefile.lib, Makefile.shlib, and Makefile.list

makelib:
	@$(SHELL) Make.sh style
	@$(SHELL) Make.sh Makefile.lib

makeshlib:
	@$(SHELL) Make.sh style
	@$(SHELL) Make.sh Makefile.shlib

makelist:
	@$(SHELL) Make.sh style
	@$(SHELL) Make.sh Makefile.list

# install LAMMPS shared lib and Python wrapper for Python usage

install-python:
	@python ../python/install.py
