/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifndef SPARTA_MASK_H
#define SPARTA_MASK_H

// data masks

#define EMPTY_MASK     0x00000000
#define ALL_MASK       0xffffffff

// particles

#define PARTICLE_MASK  0x00000001
#define SPECIES_MASK   0x00000002
#define CUSTOM_MASK    0x00004096

// grid

#define CELL_MASK      0x00000004
#define CINFO_MASK     0x00000008
#define PCELL_MASK     0x00000016
#define SINFO_MASK     0x00000032
#define PLEVEL_MASK    0x00000064

// collide

#define VREMAX_MASK    0x00000128
#define REMAIN_MASK    0x00000256

// surf

#define PT_MASK          0x00000512
#define LINE_MASK        0x00001024
#define TRI_MASK         0x00002048
#define SURF_CUSTOM_MASK 0x00008192


#endif
