#
#  Makefile.in: 
#
RM = /bin/rm -f

PROGRAM = Create_Sparta_Input
TESTNAME = $(PROGRAM)

DATA_FILES =  air.xml airSM.xml problem_setup.sh
BLESS_FILES = output_blessed.xml

TEST_FLAGS = runGood runBad runInProgress.txt

TEST_RESULTS = out.txt

TEST_INTERMEDIATES =  diff_out.txt vcs_equilibrate_res.csv err_out.txt diff_summary.txt

OBJ_FILES = $(PROGRAM).o
#
##########################################################################
#           Setup The Zuzax Directory locations
#########################################################################

all:  $(PROGRAM)
#
# Point to the Sparta MAKE include file for serial_debug. This must include the path to a Zuzax
# installation. Include the whole Makefile (or just include the directory location!)
#
MSERIALDEBUG_MAK = ../../src/MAKE/Makefile.serial_debug
include $(MSERIALDEBUG_MAK)
#
#CANTERA_INCROOT=/ascldap/users/hkmoffa/arch/linux64_gcc720_rh7/Zzuzax_master_dbg/include
#
#include $(CANTERA_INCROOT)/zuzax/Zuzax.mak
#
#
# Set up the various local and global defines
#
LDEFINES = -DDEBUG  -DDEBUG_TEST
#
# Set up the global include directory macro
#
INCLUDE_DIRS = -I../../src $(CANTERA_TOTAL_INCLUDES) 
#
# the C++ compiler
#
CXX = /projects/sierra/linux_rh7/SDK/compilers/gcc/7.2.0-RHEL6/bin/g++

debug_mode = 0
ifeq ($(debug_mode), 1)
   DEBUG_FLAG=-DDEBUG_MODE
else
   DEBUG_FLAG=
endif

CXXFLAGS = -std=c++11 -g -Wall -DDEBUG_MODE $(DEBUG_FLAG)
LCXX_FLAGS = -std=c++11 -g -Wall -DDEBUG_MODE
LCXX_END_LIBS =  -lrt -lgfortran -lm -lstdc++ 
CXXLDR = $(CXX)
CFLAGS=-std=c++11 -g -Wall -DDEBUG_MODE
#
# Both the C and C++ compiler predefines
#
CCFLAGS =   $(LDEFINES) $(INCLUDE_DIRS) 
######################################################################

all: .dependsT $(PROGRAM)

######################################################################

.SUFFIXES:  .d .cpp

.cpp.o:
	$(CXX) $(CXXFLAGS) $(CCFLAGS) -c $<
# 
# These classes compile against the Zuzax kernel
# format
#
SRCDIR=../../src

###########################################################################
# section to set up depends
###########################################################################
DEPENDS = $(PROGRAM).d

.cpp.d:
	/projects/sierra/linux_rh7/SDK/mpi/openmpi/1.10.2-gcc-7.2.0-RHEL6/bin/mpiCC -MM $(CXXFLAGS) $(CCFLAGS) -c $*.cpp > $*.d

.dependsT: $(DEPENDS)
	cat *.d > .depends
	touch .dependsT

depends:
	$(RM) *.d .depends
	$(MAKE) .dependsT

$(PROGRAM).d: $(PROGRAM).cpp Makefile

############################################################################
############################################################################
$(OBJ_FILES) : Makefile

$(PROGRAM): $(OBJ_FILES)  $(CANTERA_CORE_LIBS_DEP) 
	$(CXXLDR) $(LCXX_FLAGS) $(FLIB_BSD) -o $(PROGRAM) \
                  $(OBJ_FILES) \
                   $(CANTERA_TOTAL_LIBS) $(LCXX_END_LIBS) 



############################################################################
#           Section to set up the clean steps
############################################################################

clean:
	$(RM) $(OBJ_FILES) core *.trace $(PROGRAM) $(DEPENDS) .depends .dependsT $(TEST_FLAGS) \
	      *.trace .nsf*

test_clean:
	$(RM) $(TEST_RESULTS) $(TEST_INTERMEDIATES) $(TEST_FLAGS)

dist_clean:
	$(MAKE) -s test_clean
	$(MAKE) -s clean

###########################################################################
#          INSTALL RULES
###########################################################################

ifeq ($(wildcard .depends), .depends)
include .depends
endif

install:

