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

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "update.h"
#include "math_const.h"
#include "particle.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "domain.h"
#include "comm.h"
#include "collide.h"
#include "grid.h"
#include "surf.h"
#include "surf_collide.h"
#include "output.h"
#include "geometry.h"
#include "random_mars.h"
#include "timer.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};         // same as Domain
enum{PERIODIC,OUTFLOW,REFLECT,SURFACE,AXISYM};  // same as Domain
enum{OUTSIDE,INSIDE,ONSURF2OUT,ONSURF2IN};      // several files
enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT};   // several files
enum{NCHILD,NPARENT,NUNKNOWN,NPBCHILD,NPBPARENT,NPBUNKNOWN,NBOUND};  // Grid

#define MAXSTUCK 20

// either set ID or PROC/INDEX, set other to -1

//#define MOVE_DEBUG 1              // un-comment to debug one particle
#define MOVE_DEBUG_ID 1010784267  // particle ID
#define MOVE_DEBUG_PROC -1        // owning proc
#define MOVE_DEBUG_INDEX 13       // particle index on owning proc
#define MOVE_DEBUG_STEP 207       // timestep

/* ---------------------------------------------------------------------- */

Update::Update(SPARTA *sparta) : Pointers(sparta)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  ntimestep = 0;
  firststep = laststep = 0;
  beginstep = endstep = 0;
  runflag = 0;

  unit_style = NULL;
  set_units("si");

  fnum = 1.0;
  nrho = 1.0;
  vstream[0] = vstream[1] = vstream[2] = 0.0;
  temp_thermal = 273.15;
  gravity[0] = gravity[1] = gravity[2] = 0.0;

  maxmigrate = 0;
  mlist = NULL;

  nslist_compute = nblist_compute = 0;
  slist_compute = blist_compute = NULL;
  slist_active = blist_active = NULL;

  ranmaster = new RanMars(sparta);
}

/* ---------------------------------------------------------------------- */

Update::~Update()
{
  delete [] unit_style;
  memory->destroy(mlist);
  delete [] slist_compute;
  delete [] blist_compute;
  delete [] slist_active;
  delete [] blist_active;
  delete ranmaster;
}

/* ---------------------------------------------------------------------- */

void Update::set_units(const char *style)
{
  // physical constants from:
  // http://physics.nist.gov/cuu/Constants/Table/allascii.txt
  
  if (strcmp(style,"cgs") == 0) {
    boltz = 1.3806488e-16;
    mvv2e = 1.0;
    dt = 1.0;

  } else if (strcmp(style,"si") == 0) {
    boltz = 1.3806488e-23;
    mvv2e = 1.0;
    dt = 1.0;
    
  } else error->all(FLERR,"Illegal units command");

  delete [] unit_style;
  int n = strlen(style) + 1;
  unit_style = new char[n];
  strcpy(unit_style,style);
}

/* ---------------------------------------------------------------------- */

void Update::init()
{
  // init the Update class if performing a run, else just return
  // only set first_update if a run is being performed

  if (runflag == 0) return;
  first_update = 1;

  // choose the appropriate move method

  if (domain->dimension == 3) {
    if (surf->exist) moveptr = &Update::move<3,1>;
    else moveptr = &Update::move<3,0>;
  } else if (domain->dimension == 2) {
    if (surf->exist) moveptr = &Update::move<2,1>;
    else moveptr = &Update::move<2,0>;
  }

  // check gravity vector

  if (domain->dimension == 2 && gravity[2] != 0.0)
    error->all(FLERR,"Gravity in z not allowed for 2d");

  // moveperturb method is set if particle motion is perturbed

  moveperturb = NULL;
  if (domain->axisymmetric) moveperturb = &Update::axisymmetric;
  if (gravity[0] != 0.0 || gravity[1] != 0.0 || gravity[2] != 0.0) {
    if (domain->axisymmetric && gravity[1] != 0.0) 
      error->all(FLERR,"Gravity in y not allowed for axi-symmetry");
    if (domain->dimension == 3) moveperturb = &Update::gravity3d;
    if (domain->dimension == 2) moveperturb = &Update::gravity2d;
  }

  if (moveperturb) perturbflag = 1;
  else perturbflag = 0;
}

/* ---------------------------------------------------------------------- */

void Update::setup()
{
  // initialize counters in case stats outputs them
  // initialize running stats before each run

  ntouch_one = ncomm_one = 0;
  nboundary_one = nexit_one = 0;
  nscheck_one = nscollide_one = 0;

  niterate_running = 0;
  nmove_running = ntouch_running = ncomm_running = 0;
  nboundary_running = nexit_running = 0;
  nscheck_running = nscollide_running = 0;
  nstuck = 0;

  bounce_tally = bounce_setup();
  modify->setup();
  output->setup(1);
}

/* ---------------------------------------------------------------------- */

void Update::run(int nsteps)
{
  int n_start_of_step = modify->n_start_of_step;
  int n_end_of_step = modify->n_end_of_step;
  int dynamic = 0;

  // cellweightflag = 1 if grid-based particle weighting is ON

  int cellweightflag = 0;
  if (grid->cellweightflag) cellweightflag = 1;

  // loop over timesteps

  for (int i = 0; i < nsteps; i++) {

    ntimestep++;
    if (bounce_tally) bounce_set(ntimestep);

    // start of step fixes

    if (n_start_of_step) modify->start_of_step();

    //if (dynamic) domain->dynamic();

    // move particles

    timer->stamp();
    if (cellweightflag) particle->pre_weight();
    (this->*moveptr)();
    timer->stamp(TIME_MOVE);

    // communicate particles

    timer->stamp();
    comm->migrate_particles(nmigrate,mlist);
    if (cellweightflag) particle->post_weight();
    timer->stamp(TIME_COMM);

    if (collide) {
      timer->stamp();
      particle->sort();
      timer->stamp(TIME_SORT);

      timer->stamp();
      collide->collisions();
      timer->stamp(TIME_COLLIDE);
    }

    // diagnostic fixes

    if (n_end_of_step) modify->end_of_step();

    // all output

    if (ntimestep == output->next) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(TIME_OUTPUT);
    }
  }
}

/* ----------------------------------------------------------------------
   advect particles thru grid
   DIM = 2/3 for 2d/3d
   SURF = 0/1 for no surfs or surfs
   use multiple iterations of move/comm if necessary
------------------------------------------------------------------------- */

template < int DIM, int SURF > void Update::move()
{
  bool hitflag;
  int m,icell,icell_original,nmask,outface,bflag,nflag,pflag,pexit;
  int side,minside,minsurf,nsurf,cflag,isurf,exclude,stuck_iterate;
  int pstart,pstop,entryexit,any_entryexit;
  int *csurfs;
  cellint *neigh;
  double dtremain,frac,newfrac,param,minparam;
  double xnew[3],vold[3],xc[3],xhold[3],minxc[3];
  double *x,*v,*lo,*hi;
  Surf::Tri *tri;
  Surf::Line *line;

  // for 2d only
  // xnew,xc passed to geometry routines which use or set 3rd component

  if (DIM == 2) xnew[2] = xc[2] = 0.0;

  // extend migration list if necessary

  int nlocal = particle->nlocal;
  int maxlocal = particle->maxlocal;

  if (nlocal > maxmigrate) {
    maxmigrate = maxlocal;
    memory->destroy(mlist);
    memory->create(mlist,maxmigrate,"particle:mlist");
  }

  // counters

  niterate = 0;
  ntouch_one = ncomm_one = 0;
  nboundary_one = nexit_one = 0;
  nscheck_one = nscollide_one = 0;

  // move/migrate iterations

  Grid::ChildCell *cells = grid->cells;
  Surf::Tri *tris = surf->tris;
  Surf::Line *lines = surf->lines;
  Surf::Point *pts = surf->pts;
  double dt = update->dt;
  int notfirst = 0;

  while (1) {

    // loop over particles
    // first iteration = all my particles
    // subsequent iterations = received particles

    niterate++;
    Particle::OnePart *particles = particle->particles;
    nmigrate = 0;
    entryexit = 0;

    if (notfirst == 0) {
      notfirst = 1;
      pstart = 0;
      pstop = nlocal;
    }

    for (int i = pstart; i < pstop; i++) {
      pflag = particles[i].flag;

      // received from another proc and move is done
      // if first iteration, PDONE is from a previous step,
      //   set pflag to PKEEP so move the particle on this step
      // else do nothing

      if (pflag == PDONE) {
        pflag = particles[i].flag = PKEEP;
        if (niterate > 1) continue;
      }

      x = particles[i].x;
      v = particles[i].v;
      pexit = 0;

      if (pflag == PKEEP) {
        dtremain = dt;
        xnew[0] = x[0] + dtremain*v[0];
        xnew[1] = x[1] + dtremain*v[1];
        if (DIM == 3) xnew[2] = x[2] + dtremain*v[2];
        if (perturbflag) (this->*moveperturb)(dtremain,xnew,v);
      } else if (pflag == PINSERT) {
        dtremain = particles[i].dtremain;
        xnew[0] = x[0] + dtremain*v[0];
        xnew[1] = x[1] + dtremain*v[1];
        if (DIM == 3) xnew[2] = x[2] + dtremain*v[2];
        if (perturbflag) (this->*moveperturb)(dtremain,xnew,v);
      } else if (pflag == PENTRY) {
        icell = particles[i].icell;
        if (cells[icell].nsplit > 1) {
          if (DIM == 3 && SURF) {
            icell = split3d(icell,x);
          }
          if (DIM == 2 && SURF) {
            icell = split2d(icell,x);
          }
          particles[i].icell = icell;
        }
        dtremain = particles[i].dtremain;
        xnew[0] = x[0] + dtremain*v[0];
        xnew[1] = x[1] + dtremain*v[1];
        if (DIM == 3) xnew[2] = x[2] + dtremain*v[2];
      } else if (pflag == PEXIT) {
        pexit = 1;
        dtremain = particles[i].dtremain;
        xnew[0] = x[0];
        xnew[1] = x[1];
        if (DIM == 3) xnew[2] = x[2];
      }

      particles[i].flag = PKEEP;
      icell = particles[i].icell;
      lo = cells[icell].lo;
      hi = cells[icell].hi;
      neigh = cells[icell].neigh;
      nmask = cells[icell].nmask;
      exclude = -1;
      stuck_iterate = 0;
      ntouch_one++;

      // advect particle from cell to cell until move is done

      while (1) {

#ifdef MOVE_DEBUG
        if (DIM == 3) {
          if (ntimestep == MOVE_DEBUG_STEP && 
              (MOVE_DEBUG_ID == particles[i].id ||
               (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX))) 
            printf("PARTICLE %d %ld: %d %d: %d: x %g %g %g: xnew %g %g %g: %d " 
                   CELLINT_FORMAT ": lo %g %g %g: hi %g %g %g\n",
                   me,update->ntimestep,i,particles[i].id,
                   cells[icell].nsurf,
                   x[0],x[1],x[2],xnew[0],xnew[1],xnew[2],
                   icell,cells[icell].id,
                   lo[0],lo[1],lo[2],hi[0],hi[1],hi[2]);
        }
        if (DIM == 2) {
          if (ntimestep == MOVE_DEBUG_STEP && 
              (MOVE_DEBUG_ID == particles[i].id ||
               (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX))) 
            printf("PARTICLE %d %ld: %d %d: %d: x %g %g: xnew %g %g: %d "
                   CELLINT_FORMAT ": lo %g %g: hi %g %g\n",
                   me,update->ntimestep,i,particles[i].id,
                   cells[icell].nsurf,
                   x[0],x[1],xnew[0],xnew[1],
                   icell,cells[icell].id,
                   lo[0],lo[1],hi[0],hi[1]);
        }
#endif

        // check if particle crosses any cell face
        // frac = fraction of move completed before hitting cell face
        // check pexit to avoid divide by 0.0
        // this section should be as efficient as possible,
        // since most particles won't do anything else

        outface = INTERIOR;
        frac = 1.0;
        
        if (xnew[0] < lo[0]) {
          frac = (lo[0]-x[0]) / (xnew[0]-x[0]);
          outface = XLO;
        } else if (xnew[0] >= hi[0]) {
          if (pexit) frac = 0.0;
          else frac = (hi[0]-x[0]) / (xnew[0]-x[0]);
          outface = XHI;
        }
        
        if (xnew[1] < lo[1]) {
          newfrac = (lo[1]-x[1]) / (xnew[1]-x[1]);
          if (newfrac < frac) {
            frac = newfrac;
            outface = YLO;
          }
        } else if (xnew[1] >= hi[1]) {
          if (pexit) newfrac = 0.0;
          else newfrac = (hi[1]-x[1]) / (xnew[1]-x[1]);
          if (newfrac < frac) {
            frac = newfrac;
            outface = YHI;
          }
        }
        
        if (DIM == 3) {
          if (xnew[2] < lo[2]) {
            newfrac = (lo[2]-x[2]) / (xnew[2]-x[2]);
            if (newfrac < frac) {
              frac = newfrac;
              outface = ZLO;
            }
          } else if (xnew[2] >= hi[2]) {
            if (pexit) newfrac = 0.0;
            else newfrac = (hi[2]-x[2]) / (xnew[2]-x[2]);
            if (newfrac < frac) {
              frac = newfrac;
              outface = ZHI;
            }
          }
        }

#ifdef MOVE_DEBUG
        if (ntimestep == MOVE_DEBUG_STEP && 
            (MOVE_DEBUG_ID == particles[i].id ||
             (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX))) {
          if (outface != INTERIOR)
            printf("  OUT CHECK %d out: %d %d\n",
                   outface,grid->neigh_decode(nmask,outface),neigh[outface]);
          else 
            printf("  INT CHECK %d %d\n",outface,INTERIOR);
        }
#endif

        // START of code specific to surfaces

        if (SURF) {

          nsurf = cells[icell].nsurf;
          if (nsurf) {

            // particle crosses cell face, reset xnew exactly on face of cell
            // so surface check occurs only for particle path within grid cell
            // save xnew so can restore later if needed

            if (outface != INTERIOR) {
              xhold[0] = xnew[0];
              xhold[1] = xnew[1];
              if (DIM == 3) xhold[2] = xnew[2];
              
              xnew[0] = x[0] + frac*(xnew[0]-x[0]);
              xnew[1] = x[1] + frac*(xnew[1]-x[1]);
              if (DIM == 3) xnew[2] = x[2] + frac*(xnew[2]-x[2]);
              
              if (outface == XLO) xnew[0] = lo[0];
              else if (outface == XHI) xnew[0] = hi[0]; 
              else if (outface == YLO) xnew[1] = lo[1];
              else if (outface == YHI) xnew[1] = hi[1];
              else if (outface == ZLO) xnew[2] = lo[2];
              else if (outface == ZHI) xnew[2] = hi[2]; 
            }

            // check for collisions with triangles or lines in cell
            // find 1st surface hit via minparam
            // not considered collision if particles starts on surf, moving out
            // not considered collision if 2 params are tied and one INSIDE surf
            // if collision occurs, perform collision with surface model
            // reset x,v,xnew,dtremain and continue particle trajectory
          
            cflag = 0;
            minparam = 2.0;
            csurfs = cells[icell].csurfs;
            for (m = 0; m < nsurf; m++) {
              isurf = csurfs[m];
              if (isurf == exclude) continue;
              if (DIM == 3) {
                tri = &tris[isurf];
                hitflag = Geometry::
                  line_tri_intersect(x,xnew,
                                     pts[tri->p1].x,pts[tri->p2].x,
                                     pts[tri->p3].x,tri->norm,xc,param,side);
              }
              if (DIM == 2) {line = &lines[isurf];
                line = &lines[isurf];
                hitflag = Geometry::
                  line_line_intersect(x,xnew,
                                      pts[line->p1].x,pts[line->p2].x,
                                      line->norm,xc,param,side);
              }
              
#ifdef MOVE_DEBUG
              if (DIM == 3) {
                if (hitflag && ntimestep == MOVE_DEBUG_STEP && 
                    (MOVE_DEBUG_ID == particles[i].id ||
                     (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
                  printf("SURF COLLIDE: %d %d %d %d: "
                         "P1 %g %g %g: P2 %g %g %g: "
                         "T1 %g %g %g: T2 %g %g %g: T3 %g %g %g: "
                         "TN %g %g %g: XC %g %g %g: "
                         "Param %g: Side %d: DTR %g\n",
                         MOVE_DEBUG_INDEX,icell,nsurf,isurf,
                         x[0],x[1],x[2],xnew[0],xnew[1],xnew[2],
                         pts[tri->p1].x[0],pts[tri->p1].x[1],pts[tri->p1].x[2],
                         pts[tri->p2].x[0],pts[tri->p2].x[1],pts[tri->p2].x[2],
                         pts[tri->p3].x[0],pts[tri->p3].x[1],pts[tri->p3].x[2],
                         tri->norm[0],tri->norm[1],tri->norm[2],
                         xc[0],xc[1],xc[2],param,side,
                         dtremain);
              }
              if (DIM == 2) {
                if (hitflag && ntimestep == MOVE_DEBUG_STEP && 
                    (MOVE_DEBUG_ID == particles[i].id ||
                     (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
                  printf("SURF COLLIDE: %d %d %d %d: P1 %g %g: P2 %g %g: "
                         "L1 %g %g: L2 %g %g: LN %g %g: XC %g %g: "
                         "Param %g: Side %d: DTR %g\n",
                         MOVE_DEBUG_INDEX,icell,nsurf,isurf,
                         x[0],x[1],xnew[0],xnew[1],
                         pts[line->p1].x[0],pts[line->p1].x[1],
                         pts[line->p2].x[0],pts[line->p2].x[1],
                         line->norm[0],line->norm[1],
                         xc[0],xc[1],param,side,
                         dtremain);
              }
#endif
              
              if (hitflag && side != ONSURF2OUT && param <= minparam) {
                if (param == minparam && side == INSIDE) continue;
                cflag = 1;
                minparam = param;
                minside = side;
                minsurf = isurf;
                minxc[0] = xc[0];
                minxc[1] = xc[1];
                if (DIM == 3) minxc[2] = xc[2];
              }
            }
            nscheck_one += nsurf;
            
            if (cflag) {
              if (minside == INSIDE) {
                char str[128];
                sprintf(str,
                        "Particle %d on proc %d hit inside of "
                        "surf %d on step " BIGINT_FORMAT,
                        i,me,minsurf,update->ntimestep);
                error->one(FLERR,str);
              }
              x[0] = minxc[0];
              x[1] = minxc[1];
              if (DIM == 3) x[2] = minxc[2];
              if (DIM == 3) tri = &tris[minsurf];
              if (DIM == 2) line = &lines[minsurf];
              
              if (nsurf_tally) {
                vold[0] = v[0];
                vold[1] = v[1];
                if (DIM == 3) vold[2] = v[2];
                surf->sc[tri->isc]->collide(&particles[i],tri->norm);
                for (m = 0; m < nsurf_tally; m++)
                  slist_active[m]->surf_tally(minsurf,vold,&particles[i]);
              } else {
                if (DIM == 3)
                  surf->sc[tri->isc]->collide(&particles[i],tri->norm);
                if (DIM == 2)
                  surf->sc[line->isc]->collide(&particles[i],line->norm);
              }
              
              // nstuck = consective iterations particle is immobile

              if (minparam == 0.0) stuck_iterate++;
              else stuck_iterate = 0;

              // decrement dtremain and reset post-bounce xnew

              dtremain *= 1.0 - minparam*frac;
              xnew[0] = x[0] + dtremain*v[0];
              xnew[1] = x[1] + dtremain*v[1];
              if (DIM == 3) xnew[2] = x[2] + dtremain*v[2];
              exclude = minsurf;
              nscollide_one++;
              
#ifdef MOVE_DEBUG
              if (DIM == 3) {
                if (ntimestep == MOVE_DEBUG_STEP && 
                    (MOVE_DEBUG_ID == particles[i].id ||
                     (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
                  printf("POST COLLISION %d: %g %g %g: %g %g %g: %g %g %g\n",
                         MOVE_DEBUG_INDEX,
                         x[0],x[1],x[2],xnew[0],xnew[1],xnew[2],
                         minparam,frac,dtremain);
              }
              if (DIM == 2) {
                if (ntimestep == MOVE_DEBUG_STEP && 
                    (MOVE_DEBUG_ID == particles[i].id ||
                     (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
                  printf("POST COLLISION %d: %g %g: %g %g: %g %g %g\n",
                         MOVE_DEBUG_INDEX,
                         x[0],x[1],xnew[0],xnew[1],
                         minparam,frac,dtremain);
              }
#endif

              // if particle stuck on surfs, delete it, else continue

              if (stuck_iterate < MAXSTUCK) continue;
              particles[i].flag = PDISCARD;
              nstuck++;
            }

            // no collision, so restore saved xnew if changed it above
            
            if (outface != INTERIOR) {
              xnew[0] = xhold[0];
              xnew[1] = xhold[1];
              if (DIM == 3) xnew[2] = xhold[2];
            }
          }
        } 

        // END of code specific to surfaces

        // no cell crossing and no surface collision
        // set final particle position to xnew
        // if migrating to another proc, 
        //   flag as PDONE so new proc won't move it more on this step
        
        if (outface == INTERIOR) {
          x[0] = xnew[0];
          x[1] = xnew[1];
          if (DIM == 3) x[2] = xnew[2];
          if (cells[icell].proc != me) particles[i].flag = PDONE;
          break;
        }
          
        // particle crosses cell face
        // decrement dtremain in case particle is passed to another proc
        // reset particle position exactly on cell face
          
        exclude = -1;
        dtremain *= 1.0-frac;
          
        x[0] += frac * (xnew[0]-x[0]);
        x[1] += frac * (xnew[1]-x[1]);
        if (DIM == 3) x[2] += frac * (xnew[2]-x[2]);
          
        if (outface == XLO) x[0] = lo[0];
        else if (outface == XHI) x[0] = hi[0]; 
        else if (outface == YLO) x[1] = lo[1];
        else if (outface == YHI) x[1] = hi[1];
        else if (outface == ZLO) x[2] = lo[2];
        else if (outface == ZHI) x[2] = hi[2]; 

        // nflag = type of neighbor cell: child, parent, unknown, boundary
        // if parent, use id_find_child to identify child cell
        //   result of id_find_child could be unknown:
        //     particle is hitting face of a ghost child cell which extends
        //     beyond my ghost halo, cell on other side of face is a parent,
        //     it's child which the particle is in is entirely beyond my halo
        // if new cell is child and surfs exist, check if a split cell

        nflag = grid->neigh_decode(nmask,outface);
        icell_original = icell;

        if (nflag == NCHILD) {
          icell = neigh[outface];
          if (DIM == 3 && SURF) {
            if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
              icell = split3d(icell,x);
          }
          if (DIM == 2 && SURF) {
            if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
              icell = split2d(icell,x);
          }
        } else if (nflag == NPARENT) {
          icell = grid->id_find_child(neigh[outface],x);
          if (icell >= 0) {
            if (DIM == 3 && SURF) {
              if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
                icell = split3d(icell,x);
            }
            if (DIM == 2 && SURF) {
              if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
                icell = split2d(icell,x);
            }
          }
        } else if (nflag == NUNKNOWN) icell = -1;
        
        // neighbor cell is global boundary
        // tally boundary stats if requested
        // collide() updates x,v,xnew as needed due to boundary interaction
        // OUTFLOW: exit with particle flag = PDISCARD
        // PERIODIC: new cell via same logic as above for child/parent/unknown
        // other = reflected particle stays in same grid cell

        else {
          if (nboundary_tally) {
            vold[0] = v[0]; 
            vold[1] = v[1]; 
            if (DIM == 3) vold[2] = v[2];
            bflag = domain->collide(&particles[i],outface,icell,xnew);
            for (m = 0; m < nboundary_tally; m++)
              blist_active[m]->
                boundary_tally(outface,bflag,vold,&particles[i]);
          } else
            bflag = domain->collide(&particles[i],outface,icell,xnew);

          if (bflag == OUTFLOW) {
            particles[i].flag = PDISCARD;
            nexit_one++;
            break;

          } else if (bflag == PERIODIC) {
            if (nflag == NPBCHILD) {
              icell = neigh[outface];
              if (DIM == 3 && SURF) {
                if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
                  icell = split3d(icell,x);
              }
              if (DIM == 2 && SURF) {
                if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
                  icell = split2d(icell,x);
              }
            } else if (nflag == NPBPARENT) {
              icell = grid->id_find_child(neigh[outface],x);
              if (icell >= 0) {
                if (DIM == 3 && SURF) {
                  if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
                    icell = split3d(icell,x);
                }
                if (DIM == 2 && SURF) {
                  if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
                    icell = split2d(icell,x);
                }
              } else domain->uncollide(outface,x);
            } else if (nflag == NPBUNKNOWN) {
              icell = -1;
              domain->uncollide(outface,x);
            }

          } else {
            nboundary_one++;
            ntouch_one--;    // decrement here since will increment below
          }
        }

        // neighbor cell is unknown
        // reset icell to original icell which must be a ghost cell
        // exit with particle flag = PEXIT, so receiver can identify neighbor

        if (icell < 0) {
          icell = icell_original;
          particles[i].flag = PEXIT;
          particles[i].dtremain = dtremain;
          entryexit = 1;
          break;
        }

        // if nsurf < 0, new cell is EMPTY ghost
        // exit with particle flag = PENTRY, so receiver can continue move
        
        if (cells[icell].nsurf < 0) {
          particles[i].flag = PENTRY;
          particles[i].dtremain = dtremain;
          entryexit = 1;
          break;
        }

        // move particle into new grid cell for next stage of move

        lo = cells[icell].lo;
        hi = cells[icell].hi;
        neigh = cells[icell].neigh;
        nmask = cells[icell].nmask;
        ntouch_one++;
      }

#ifdef MOVE_DEBUG
      if (ntimestep == MOVE_DEBUG_STEP && 
          (MOVE_DEBUG_ID == particles[i].id ||
           (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
        printf("MOVE DONE %d %d %d: %g %g %g: %g\n",
               MOVE_DEBUG_INDEX,particles[i].flag,icell,
               x[0],x[1],x[2],particles[i].dtremain);
#endif

      // move is complete, or as much as can be done on this proc
      // update particle's grid cell
      // if particle flag set, add particle to migrate list
      // if discarding, migration will delete particle
    
      particles[i].icell = icell;
      
      if (particles[i].flag != PKEEP) {
        mlist[nmigrate++] = i;
        if (particles[i].flag != PDISCARD) {
          if (cells[icell].proc == me)
            error->one(FLERR,"Sending particle to self");
          ncomm_one++;
        }
      }

      // END of while loop advecting one particle

    }

    // END of pstart/pstop loop advecting all particles
    
    // if gridcut >= 0.0, check if another iteration of move is required
    // only the case if some particle flag = PENTRY/PEXIT
    //   in which case perform particle migration
    // if not, move is done and final particle comm will occur in run()
    // if iterating, reset pstart/pstop and extend migration list if necessary

    if (grid->cutoff < 0.0) break;

    MPI_Allreduce(&entryexit,&any_entryexit,1,MPI_INT,MPI_MAX,world);
    if (any_entryexit) {
      pstart = comm->migrate_particles(nmigrate,mlist);
      pstop = particle->nlocal;
      if (pstop-pstart > maxmigrate) {
        maxmigrate = pstop-pstart;
        memory->destroy(mlist);
        memory->create(mlist,maxmigrate,"particle:mlist");
      }
    } else break;

    // END of single move/migrate iteration

  }

  // END of all move/migrate iterations

  // accumulate running totals

  niterate_running += niterate;
  nmove_running += nlocal;
  ntouch_running += ntouch_one;
  ncomm_running += ncomm_one;
  nboundary_running += nboundary_one;
  nexit_running += nexit_one;
  nscheck_running += nscheck_one;
  nscollide_running += nscollide_one;
}

/* ----------------------------------------------------------------------
   particle is entering split parent icell at x
   determine which split child cell it is in
   return ccell = index of sub-cell in ChildCell
------------------------------------------------------------------------- */

int Update::split3d(int icell, double *x)
{
  int m,cflag,isurf,hitflag,side,minside,minsurfindex;
  double param,minparam;
  double xc[3];
  Surf::Tri *tri;

  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;
  Surf::Tri *tris = surf->tris;
  Surf::Point *pts = surf->pts;

  // check for collisions with lines in cell
  // find 1st surface hit via minparam
  // only consider tris that are mapped via csplits to a split cell
  //   unmapped tris only touch cell surf at xnew
  //   another mapped tri should include same xnew
  // not considered a collision if particles starts on surf, moving out
  // not considered a collision if 2 params are tied and one is INSIDE surf

  int nsurf = cells[icell].nsurf;
  int *csurfs = cells[icell].csurfs;
  int isplit = cells[icell].isplit;
  int *csplits = sinfo[isplit].csplits;
  double *xnew = sinfo[isplit].xsplit;

  cflag = 0;
  minparam = 2.0;
  for (m = 0; m < nsurf; m++) {
    if (csplits[m] < 0) continue;
    isurf = csurfs[m];
    tri = &tris[isurf];
    hitflag = Geometry::
      line_tri_intersect(x,xnew,
                         pts[tri->p1].x,pts[tri->p2].x,pts[tri->p3].x,
                         tri->norm,xc,param,side);
    
    if (hitflag && side != INSIDE && param < minparam) {
      cflag = 1;
      minparam = param;
      minside = side;
      minsurfindex = m;
    }
  }

  if (!cflag) return sinfo[isplit].csubs[sinfo[isplit].xsub];
  int index = csplits[minsurfindex];
  return sinfo[isplit].csubs[index];
}

/* ----------------------------------------------------------------------
   particle is entering split ICELL at X
   determine which split sub-cell it is in
   return ccell = index of sub-cell in ChildCell
------------------------------------------------------------------------- */

int Update::split2d(int icell, double *x)
{
  int m,cflag,isurf,hitflag,side,minside,minsurfindex;
  double param,minparam;
  double xc[3];
  Surf::Line *line;

  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;
  Surf::Line *lines = surf->lines;
  Surf::Point *pts = surf->pts;

  // check for collisions with lines in cell
  // find 1st surface hit via minparam
  // only consider lines that are mapped via csplits to a split cell
  //   unmapped lines only touch cell surf at xnew
  //   another mapped line should include same xnew
  // not considered a collision if particle starts on surf, moving out
  // not considered a collision if 2 params are tied and one is INSIDE surf

  int nsurf = cells[icell].nsurf;
  int *csurfs = cells[icell].csurfs;
  int isplit = cells[icell].isplit;
  int *csplits = sinfo[isplit].csplits;
  double *xnew = sinfo[isplit].xsplit;

  cflag = 0;
  minparam = 2.0;
  for (m = 0; m < nsurf; m++) {
    if (csplits[m] < 0) continue;
    isurf = csurfs[m];
    line = &lines[isurf];
    hitflag = Geometry::
      line_line_intersect(x,xnew,
                          pts[line->p1].x,pts[line->p2].x,line->norm,
                          xc,param,side);
    
    if (hitflag && side != INSIDE && param < minparam) {
      cflag = 1;
      minparam = param;
      minside = side;
      minsurfindex = m;
    }
  }

  if (!cflag) return sinfo[isplit].csubs[sinfo[isplit].xsub];
  int index = csplits[minsurfindex];
  return sinfo[isplit].csubs[index];
}

/* ----------------------------------------------------------------------
   setup lists of all computes that tally surface and boundary bounce info
   return 1 if there are any, 0 if not
------------------------------------------------------------------------- */

int Update::bounce_setup()
{
  delete [] slist_compute;
  delete [] blist_compute;
  delete [] slist_active;
  delete [] blist_active;

  slist_compute = blist_compute = NULL;
  nslist_compute = nblist_compute = 0;
  for (int i = 0; i < modify->ncompute; i++) {
    if (modify->compute[i]->surf_tally_flag) nslist_compute++;
    if (modify->compute[i]->boundary_tally_flag) nblist_compute++;
  }

  if (nslist_compute) slist_compute = new Compute*[nslist_compute];
  if (nblist_compute) blist_compute = new Compute*[nblist_compute];
  if (nslist_compute) slist_active = new Compute*[nslist_compute];
  if (nblist_compute) blist_active = new Compute*[nblist_compute];

  nslist_compute = nblist_compute = 0;
  for (int i = 0; i < modify->ncompute; i++) {
    if (modify->compute[i]->surf_tally_flag)
      slist_compute[nslist_compute++] = modify->compute[i];
    if (modify->compute[i]->boundary_tally_flag)
      blist_compute[nblist_compute++] = modify->compute[i];
  }

  if (nslist_compute || nblist_compute) return 1;
  nsurf_tally = nboundary_tally = 0;
  return 0;
}

/* ----------------------------------------------------------------------
   set bounce tally flags for current timestep
   tally flag = # of computes needing bounce info on this step
   clear accumulators in computes that will be invoked this step
------------------------------------------------------------------------- */

void Update::bounce_set(bigint ntimestep)
{
  int i,j;

  nsurf_tally = 0;
  for (i = 0; i < nslist_compute; i++)
    if (slist_compute[i]->matchstep(ntimestep)) {
      slist_active[nsurf_tally++] = slist_compute[i];
      slist_compute[i]->clear();
    }

  nboundary_tally = 0;
  for (i = 0; i < nblist_compute; i++)
    if (blist_compute[i]->matchstep(ntimestep)) {
      blist_active[nboundary_tally++] = blist_compute[i];
      blist_compute[i]->clear();
    }
}

/* ----------------------------------------------------------------------
   set global properites via global command in input script
------------------------------------------------------------------------- */

void Update::global(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal global command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"fnum") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal global command");
      fnum = atof(arg[iarg+1]);
      if (fnum <= 0.0) error->all(FLERR,"Illegal global command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"nrho") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal global command");
      nrho = atof(arg[iarg+1]);
      if (nrho <= 0.0) error->all(FLERR,"Illegal global command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"vstream") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal global command");
      vstream[0] = atof(arg[iarg+1]);
      vstream[1] = atof(arg[iarg+2]);
      vstream[2] = atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"temp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal global command");
      temp_thermal = atof(arg[iarg+1]);
      if (temp_thermal <= 0.0) error->all(FLERR,"Illegal global command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"gravity") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal global command");
      double gmag = atof(arg[iarg+1]);
      gravity[0] = atof(arg[iarg+2]);
      gravity[1] = atof(arg[iarg+3]);
      gravity[2] = atof(arg[iarg+4]);
      if (gmag < 0.0) error->all(FLERR,"Illegal global command");
      if (gmag > 0.0 && 
          gravity[0] == 0.0 && gravity[1] == 0.0 && gravity[2] == 0.0)
        error->all(FLERR,"Illegal global command");
      if (gmag > 0.0) MathExtra::snorm3(gmag,gravity);
      else gravity[0] = gravity[1] = gravity[2] = 0.0;
      iarg += 5;
    } else if (strcmp(arg[iarg],"surfmax") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal global command");
      if (surf->exist) 
        error->all(FLERR,
                   "Cannot set global surfmax when surfaces already exist");
      grid->maxsurfpercell = atoi(arg[iarg+1]);
      if (grid->maxsurfpercell <= 0) error->all(FLERR,"Illegal global command");
      grid->allocate_surf_arrays();
      iarg += 2;
    } else if (strcmp(arg[iarg],"gridcut") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal global command");
      if (grid->exist) 
        error->all(FLERR,
                   "Cannot set global gridcut when grid already exists");
      grid->cutoff = atof(arg[iarg+1]);
      if (grid->cutoff < 0.0 && grid->cutoff != -1.0)
        error->all(FLERR,"Illegal global command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"comm/sort") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal global command");
      if (strcmp(arg[iarg+1],"yes") == 0) comm->commsortflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) comm->commsortflag = 0;
      else error->all(FLERR,"Illegal global command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"comm/style") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal global command");
      if (strcmp(arg[iarg+1],"neigh") == 0) comm->commpartstyle = 1;
      else if (strcmp(arg[iarg+1],"all") == 0) comm->commpartstyle = 0;
      else error->all(FLERR,"Illegal global command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"weight") == 0) {
      // for now assume just one arg after "cell"
      // may need to generalize later
      if (iarg+3 > narg) error->all(FLERR,"Illegal global command");
      if (strcmp(arg[iarg+1],"cell") == 0) grid->weight(1,&arg[iarg+2]);
      else error->all(FLERR,"Illegal weight command");
      iarg += 3;
    } else error->all(FLERR,"Illegal global command");
  }
}
