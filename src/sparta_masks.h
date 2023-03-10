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
#define CUSTOM_MASK    0x00000004

// grid

#define CELL_MASK      0x00000008
#define CINFO_MASK     0x00000010
#define PCELL_MASK     0x00000020
#define SINFO_MASK     0x00000040
#define PLEVEL_MASK    0x00000080

// collide

#define VREMAX_MASK    0x00000100
#define REMAIN_MASK    0x00000200

// surf

#define PT_MASK          0x00000400
#define LINE_MASK        0x00000800
#define TRI_MASK         0x00001000
#define SURF_CUSTOM_MASK 0x00002000


#endif
