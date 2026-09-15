/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "variable_kokkos.h"
#include "particle_kokkos.h"
#include "grid_kokkos.h"
#include "surf_kokkos.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

// the 2 enums below are copies of the enums in variable.cpp
// they must be kept consistent with that file

enum{PARTICLE_CUSTOM,GRID_CUSTOM,SURF_CUSTOM};

enum{DONE,ADD,SUBTRACT,MULTIPLY,DIVIDE,CARAT,MODULO,UNARY,
     NOT,EQ,NE,LT,LE,GT,GE,AND,OR,
     SQRT,EXP,LN,LOG,ABS,SIN,COS,TAN,ASIN,ACOS,ATAN,ATAN2,ERF,
     RANDOM,NORMAL,CEIL,FLOOR,ROUND,RAMP,STAGGER,LOGFREQ,STRIDE,
     VDISPLACE,SWIGGLE,CWIGGLE,
     PYWRAPPER,
     VALUE,ARRAY,ARRAYINT,PARTARRAYDOUBLE,PARTARRAYINT,SPECARRAY,PARTGRIDARRAY};

/* ---------------------------------------------------------------------- */

VariableKokkos::VariableKokkos(SPARTA *sparta) : Variable(sparta)
{

}

/* ----------------------------------------------------------------------
   sync a custom particle, grid, or surf attribute to the host
   called from Variable::evaluate() when a p_name, g_name, or s_name
     custom attribute is encountered in a formula
------------------------------------------------------------------------- */

void VariableKokkos::custom_sync(int cwhich)
{
  if (cwhich == PARTICLE_CUSTOM)
    ((ParticleKokkos*)particle)->sync(Host,CUSTOM_MASK);
  else if (cwhich == GRID_CUSTOM)
    ((GridKokkos*)grid)->sync(Host,CUSTOM_MASK);
  else if (cwhich == SURF_CUSTOM)
    ((SurfKokkos*)surf)->sync(Host,CUSTOM_MASK);
}

/* ----------------------------------------------------------------------
   evaluate a particle-style variable parse tree for particle I
     or a grid-style variable parse tree for grid cell I
     or a surf-style variable parse tree for surf element I
   tree was created by one-time parsing of formula string via evaulate()
   customize by adding a function:
     sqrt(),exp(),ln(),log(),sin(),cos(),tan(),asin(),acos(),atan(),
     atan2(y,x),random(x,y),normal(x,y),ceil(),floor(),round(),
     ramp(x,y),stagger(x,y),logfreq(x,y,z),stride(x,y,z),
     vdisplace(x,y),swiggle(x,y,z),cwiggle(x,y,z)
---------------------------------------------------------------------- */

double VariableKokkos::eval_tree(Tree *tree, int i)
{
  ParticleKokkos* particle_kk = ((ParticleKokkos*)particle);

  if (tree->type == SPECARRAY)
    particle_kk->sync(Host,PARTICLE_MASK);

  if (tree->type == PARTGRIDARRAY)
    particle_kk->sync(Host,PARTICLE_MASK);

  return Variable::eval_tree(tree,i);
}

/* ----------------------------------------------------------------------
   process a particle vector in formula
   push result onto tree
   word = particle vector
   customize by adding a particle vector:
     id,x,y,z,vx,vy,vz,type,mass,q,mu
------------------------------------------------------------------------- */

void VariableKokkos::particle_vector(char *word, Tree **tree,
                               Tree **treestack, int &ntreestack)
{
  ParticleKokkos* particle_kk = ((ParticleKokkos*)particle);

  if (strcmp(word,"x") == 0)
    particle_kk->sync(Host,PARTICLE_MASK);
  else if (strcmp(word,"y") == 0)
    particle_kk->sync(Host,PARTICLE_MASK);
  else if (strcmp(word,"z") == 0)
    particle_kk->sync(Host,PARTICLE_MASK);
  else if (strcmp(word,"vx") == 0)
    particle_kk->sync(Host,PARTICLE_MASK);
  else if (strcmp(word,"vy") == 0)
    particle_kk->sync(Host,PARTICLE_MASK);
  else if (strcmp(word,"vz") == 0)
    particle_kk->sync(Host,PARTICLE_MASK);

  else if (strcmp(word,"id") == 0)
    particle_kk->sync(Host,PARTICLE_MASK);
  else if (strcmp(word,"type") == 0)
    particle_kk->sync(Host,PARTICLE_MASK);
  else if (strcmp(word,"mass") == 0)
    particle_kk->sync(Host,SPECIES_MASK);
  else if (strcmp(word,"q") == 0)
    particle_kk->sync(Host,SPECIES_MASK);
  else if (strcmp(word,"mu") == 0)
    particle_kk->sync(Host,SPECIES_MASK);

  Variable::particle_vector(word,tree,treestack,ntreestack);
}

/* ----------------------------------------------------------------------
   process a grid vector in formula
   push result onto tree
   word = grid vector
   customize by adding a grid vector:
     cxlo,cxhi,cylo,cyhi,czlo,czhi
------------------------------------------------------------------------- */

void VariableKokkos::grid_vector(char *word, Tree **tree,
                           Tree **treestack, int &ntreestack)
{
  GridKokkos* grid_kk = ((GridKokkos*)grid);
  grid_kk->sync(Host,CELL_MASK);

  Variable::grid_vector(word,tree,treestack,ntreestack);
}
