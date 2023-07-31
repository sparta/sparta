## src

* Files related to a new compute momentum command, and two new fixes: `fix_emit_surf_timeavg` and `fix_emit_surf_mflow`.
* The two fixes are based on `fix_emit_surf` of the SPARTA vanilla implementation, which is included in the path for easy comparision.
* MAKE/Makefile.kokkos_mpi_only: `SPARTA_INC = -DSPARTA_GZIP –DSPARTA_BIGBIG`. The BIGBIG flag is added.

## tests

* Two test directories `surf_mflow_test` and `surf_timeavg_test` for testing the two new fixes
* Input files and submission (to PBS scheduler) files included

## tools

* minor bugfixes made in `log2txt.py` and  `pizza/olog.py`


============================================================================================================================================
This is the SPARTA software package.

SPARTA stands for Stochastic PArallel Rarefied-gas Time-accurate
Analyzer.

Copyright (2014) Sandia Corporation.  Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software.  This software is distributed under
the GNU General Public License.

----------------------------------------------------------------------

SPARTA is a Direct Simulation Monte Carlo (DSMC) code designed to run
efficiently on parallel computers.  It was developed at Sandia
National Laboratories, a US Department of Energy facility, with
funding from the DOE.  It is an open-source code, distributed freely
under the terms of the GNU Public License (GPL).

The primary authors of the code are Steve Plimpton and Michael Gallis,
who can be emailed at sjplimp@sandia.gov and magalli@sandia.gov.  The
SPARTA web site at http://sparta.sandia.gov has more information about
the code and its uses.

The SPARTA distribution includes the following files and directories:

README			   this file
LICENSE			   the GNU General Public License (GPL)
bench                      benchmark problems
data                       files with species/reaction params, surface files
doc                        documentation
examples                   simple test problems
lib                        additional library files
python                     Python wrapper on SPARTA as a library
src                        source files
tools                      pre- and post-processing tools

Point your browser at any of these files to get started:

doc/Manual.html	           the SPARTA manual
doc/Section_intro.html	   hi-level introduction to SPARTA
doc/Section_start.html	   how to build and use SPARTA
doc/Developer.pdf          SPARTA developer guide
