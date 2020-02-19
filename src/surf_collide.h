/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifndef SPARTA_SURF_COLLIDE_H
#define SPARTA_SURF_COLLIDE_H

#include "pointers.h"
#include "particle.h"

namespace SPARTA_NS {

class SurfCollide : protected Pointers {
 public:
  char *id;
  char *style;
 
  int dynamicflag;          // 1 if any param is dynamically updated
  int allowreact;           // 1 if allows for surface reactions
  int transparent;          // 1 if transparent collision model
  int vector_flag;          // 0/1 if compute_vector() function exists
  int size_vector;          // length of global vector

  SurfCollide(class SPARTA *, int, char **);
  SurfCollide(class SPARTA *sparta) : Pointers(sparta) {}
  virtual ~SurfCollide();
  virtual void init();

  //! Main member function for SurfCollide class where one particle 
  //! collides with a surface.
  /*!
   *  @param[in,out]  ipart       Reference to the pointer to the incoming particle
   *  @param[in]      norm        value of the surface normal
   *  @param[in,out]  dtremain    Remaining time
   *  @param[in]      isr         Index of the surface collision model
   *  @param[in,out]  surfaceState   Pointer to the surface state
   */
  virtual Particle::OnePart *collide(Particle::OnePart *& ipart, double * norm, 
                                     double & dtremain, int isr, void * surfaceState = nullptr) = 0;

  virtual void dynamic() {}
  void tally_update();
  double compute_vector(int i);

  int copy,copymode;

 protected:
  //! Number of collisions during the current time step
  int nsingle;
  //! Number of total collisions
  int ntotal;
  double one[2],all[2];
};

}

#endif

/* ERROR/WARNING messages:

E: Surf_collide ID must be alphanumeric or underscore characters

Self-explanatory.

*/
