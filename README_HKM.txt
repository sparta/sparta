This is a clone of the gitlab project spartaZZ


The zuzax_tools directory contains additions to sparta involved with the coupling of the
two codes.


Coupling of Zuzax and Sparta
--------------------------------------------

Zuzax is an external optional library. Zuzax may be linked into Sparta by
Defining the Makefile variable, ZUZAX_HOME, which is the top of the Zuzax install directory.
Then, Sparta's Makefile will include the Makefile include, Zuzax.mak,
which defines the additions to compilation environment necessary to compile and
link against the Zuzax libraries.

For example, below are the Makefile additions necessary:

ZUZAX_HOME = /home/hkmoffa/arch/linux64_gcc720_rh7/Zzuzax_master_dbg
include $(ZUZAX_HOME)/include/zuzax/Cantera.mak

ZUZAX_DEF = -DUSE_ZUZAX
ZSURF_DEF = -DUSE_ZSURF
ZUZAX_INC = -I$(ZUZAX_HOME)/include $(ZUZAX_DEF) $(ZSURF_DEF)
ZUZAX_PATH = -L$(ZUZAX_HOME)/lib
ZUZAX_LIB = $(CANTERA_TOTAL_LIBS)
   

All calls to Zuzax within Sparta are enclosed within one of two ifdef blocks.
The USE_ZUZAX ifdef encloses general calls having to do with the calculation of
thermodynamics routines and the ability to get species data from the Zuzax databases.

The USE_ZSURF ifdef is used to calculate surface reactions whose surfaces are controlled
by Zuzax.

