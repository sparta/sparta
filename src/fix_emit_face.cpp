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

#include "stdlib.h"
#include "string.h"
#include "fix_emit_face.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "grid.h"
#include "surf.h"
#include "particle.h"
#include "mixture.h"
#include "comm.h"
#include "modify.h"
#include "geometry.h"
#include "random_park.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};         // same as Domain
enum{PERIODIC,OUTFLOW,REFLECT,SURFACE,AXISYM};  // same as Domain
enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};           // same as Grid
enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};   // several files
enum{NCHILD,NPARENT,NUNKNOWN,NPBCHILD,NPBPARENT,NPBUNKNOWN,NBOUND};  // Grid

#define DELTATASK 256

/* ---------------------------------------------------------------------- */

FixEmitFace::FixEmitFace(SPARTA *sparta, int narg, char **arg) :
  FixEmit(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix emit/face command");

  imix = particle->find_mixture(arg[2]);
  if (imix < 0) error->all(FLERR,"Fix emit/face mixture ID does not exist");

  // flag specified faces

  faces[XLO] = faces[XHI] = faces[YLO] = faces[YHI] =
    faces[ZLO] = faces[ZHI] = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"all") == 0) {
      if (domain->dimension == 3)
	faces[XLO] = faces[XHI] = faces[YLO] = faces[YHI] =
	  faces[ZLO] = faces[ZHI] = 1;
      else faces[XLO] = faces[XHI] = faces[YLO] = faces[YHI] = 1;
    } else if (strcmp(arg[iarg],"xlo") == 0) faces[XLO] = 1;
    else if (strcmp(arg[iarg],"xhi") == 0) faces[XHI] = 1;
    else if (strcmp(arg[iarg],"ylo") == 0) faces[YLO] = 1;
    else if (strcmp(arg[iarg],"yhi") == 0) faces[YHI] = 1;
    else if (strcmp(arg[iarg],"zlo") == 0) faces[ZLO] = 1;
    else if (strcmp(arg[iarg],"zhi") == 0) faces[ZHI] = 1;
    else break;
    iarg++;
  }

  // optional args

  np = 0;

  options(narg-iarg,&arg[iarg]);

  // error checks

  if (domain->dimension == 2 && (faces[ZLO] || faces[ZHI])) 
    error->all(FLERR,"Cannot use fix emit/face in z dimension "
               "for 2d simulation");
  if (domain->axisymmetric && faces[YLO]) 
    error->all(FLERR,"Cannot use fix emit/face on ylo face for "
               "axisymmetric model");
  if (np > 0 && perspecies) 
    error->all(FLERR,"Cannot use fix emit/face n > 0 with perspecies yes");

  // task list

  tasks = NULL;
  ntask = ntaskmax = 0;
}

/* ---------------------------------------------------------------------- */

FixEmitFace::~FixEmitFace()
{
  for (int i = 0; i < ntaskmax; i++) delete [] tasks[i].ntargetsp;
  memory->sfree(tasks);
}

/* ---------------------------------------------------------------------- */

void FixEmitFace::init()
{
  // copies of class data before invoking parent init() and count_task()

  dimension = domain->dimension;
  fnum = update->fnum;
  dt = update->dt;

  nspecies = particle->mixture[imix]->nspecies;
  nrho = particle->mixture[imix]->nrho;
  temp_thermal = particle->mixture[imix]->temp_thermal;
  vstream = particle->mixture[imix]->vstream;
  vscale = particle->mixture[imix]->vscale;
  fraction = particle->mixture[imix]->fraction;
  cummulative = particle->mixture[imix]->cummulative;

  pts = surf->pts;
  lines = surf->lines;
  tris = surf->tris;

  // cannot inflow thru periodic boundary

  for (int i = 0; i < 6; i++)
    if (faces[i] && domain->bflag[i] == PERIODIC)
      error->all(FLERR,"Cannot use fix emit/face on periodic boundary");

  // cannot have inflow on yhi if axisymmetric

  if (domain->axisymmetric && faces[YHI] && vstream[1] != 0.0)
    error->all(FLERR,"Cannot use fix emit on axisymmetric yhi "
               "if streaming velocity has a y-component");

  // warn if any inflow face does not have an inward normal
  //   in direction of streaming velocity

  double normal[3];
  int flag = 0;

  for (int i = 0; i < 6; i++) {
    if (!faces[i]) continue;
    normal[0] = normal[1] = normal[2] = 0.0;
    if (i % 2 == 0) normal[i/2] = 1.0;
    else normal[i/2] = -1.0;
    double indot = vstream[0]*normal[0] + vstream[1]*normal[1] + 
      vstream[2]*normal[2];
    if (indot < 0.0) flag = 1;
  }

  if (flag && comm->me == 0)
    error->warning(FLERR,
                   "One or more fix inflow faces oppose streaming velocity");

  // reallocate ntargetsp for each task 
  // b/c nspecies count of mixture may have changed

  for (int i = 0; i < ntask; i++) {
    delete [] tasks[i].ntargetsp;
    tasks[i].ntargetsp = new double[nspecies];
  }

  // invoke FixEmit::init() to populate task list
  // it calls create_task() for each grid cell

  ntask = 0;
  FixEmit::init();

  // if Np > 0, nper = # of insertions per task
  // set nthresh so as to achieve exactly Np insertions
  // tasks > tasks_with_no_extra need to insert 1 extra particle
  // NOTE: setting same # of insertions per task
  //       should weight by cell face area

  if (np > 0) {
    int all,nupto,tasks_with_no_extra;
    MPI_Allreduce(&ntask,&all,1,MPI_INT,MPI_SUM,world);
    if (all) {
      npertask = np / all;
      tasks_with_no_extra = all - (np % all);
    } else npertask = tasks_with_no_extra = 0;

    MPI_Scan(&ntask,&nupto,1,MPI_INT,MPI_SUM,world);
    if (tasks_with_no_extra < nupto-ntask) nthresh = 0;
    else if (tasks_with_no_extra >= nupto) nthresh = ntask;
    else nthresh = tasks_with_no_extra - (nupto-ntask);
  }
}

/* ---------------------------------------------------------------------- */

int FixEmitFace::create_task(int icell)
{
  int i,j,n,iface,flag,isp,extflag;
  int *cflags;
  double indot,area;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  // corners[i][j] = J corner points of face I of a grid cell
  // works for 2d quads and 3d hexes
  
  int corners[6][4] = {{0,2,4,6}, {1,3,5,7}, {0,1,4,5}, {2,3,6,7}, 
		       {0,1,2,3}, {4,5,6,7}};
  int nface_pts = 4;
  if (domain->dimension == 2) nface_pts = 2;

  // loop over 6 faces of icell

  int ntaskorig = ntask;
  int nmask = cells[icell].nmask;

  for (i = 0; i < 6; i++) {
    if (i == 0) iface = XLO;
    else if (i == 1) iface = XHI;
    else if (i == 2) iface = YLO;
    else if (i == 3) iface = YHI;
    else if (i == 4) iface = ZLO;
    else if (i == 5) iface = ZHI;

    // flag = 1 if insertion happens on iface of cell
    // only if face adjoins global boundary with inflow defined
    // if cell is OVERLAP:
    //   allow if any face corner point is OUTSIDE and none is INSIDE
    //   disallow if any pt of any line/tri in cell touches face

    flag = 0;
    if (faces[iface] && grid->neigh_decode(nmask,iface) == NBOUND) {
      if (cinfo[icell].type == OUTSIDE) flag = 1;
      else if (cinfo[icell].type == OVERLAP) {
        flag = 1;
        cflags = cinfo[icell].corner;
        
        extflag = 0;
        for (j = 0; j < nface_pts; j++) {
          if (cflags[corners[iface][j]] == OUTSIDE) extflag = 1;
          else if (cflags[corners[iface][j]] == INSIDE) flag = 0;
        }
        if (!extflag) flag = 0;
      
        if (flag && dimension == 2) {
          for (j = 0; j < cells[icell].nsurf; j++) {
            n = cells[icell].csurfs[j];
            if (Geometry::
                line_quad_face_touch(pts[lines[n].p1].x,
                                     pts[lines[n].p2].x,
                                     iface,cells[icell].lo,cells[icell].hi)) {
              flag = 0;
              break;
            }
          }
        } else if (flag && dimension == 3) {
          for (j = 0; j < cells[icell].nsurf; j++) {
            n = cells[icell].csurfs[j];
            if (Geometry::
                tri_hex_face_touch(pts[tris[n].p1].x,
                                   pts[tris[n].p2].x,
                                   pts[tris[n].p3].x,
                                   iface,cells[icell].lo,cells[icell].hi)) {
              flag = 0;
              break;
            }
          }
        }
      }
    }
    
    // no insertions on this face

    if (!flag) continue;

    // set cell parameters of task
    // pcell = sub cell for particles if a split cell

    if (ntask == ntaskmax) grow_task();

    tasks[ntask].icell = icell;
    tasks[ntask].iface = iface;
    if (cells[icell].nsplit > 1) tasks[ntask].pcell = split(icell,iface);
    else tasks[ntask].pcell = icell;

    // set face-dependent params of task

    tasks[ntask].lo[0] = cells[icell].lo[0];
    tasks[ntask].hi[0] = cells[icell].hi[0];
    tasks[ntask].lo[1] = cells[icell].lo[1];
    tasks[ntask].hi[1] = cells[icell].hi[1];
    tasks[ntask].lo[2] = cells[icell].lo[2];
    tasks[ntask].hi[2] = cells[icell].hi[2];
    if (dimension == 2) tasks[ntask].lo[2] = tasks[ntask].hi[2] = 0.0;
    tasks[ntask].normal[0] = 0.0;
    tasks[ntask].normal[1] = 0.0;
    tasks[ntask].normal[2] = 0.0;

    if (iface == XLO || iface == XHI) {
      tasks[ntask].ndim = 0;
      tasks[ntask].pdim = 1;
      tasks[ntask].qdim = 2;
      if (iface == XLO) tasks[ntask].hi[0] = cells[icell].lo[0];
      else tasks[ntask].lo[0] = cells[icell].hi[0];
      if (iface == XLO) tasks[ntask].normal[0] = 1.0;
      else tasks[ntask].normal[0] = -1.0;
    } else if (iface == YLO || iface == YHI) {
      tasks[ntask].ndim = 1;
      tasks[ntask].pdim = 0;
      tasks[ntask].qdim = 2;
      if (iface == YLO) tasks[ntask].hi[1] = cells[icell].lo[1];
      else tasks[ntask].lo[1] = cells[icell].hi[1];
      if (iface == YLO) tasks[ntask].normal[1] = 1.0;
      else tasks[ntask].normal[1] = -1.0;
    } else if (iface == ZLO || iface == ZHI) {
      tasks[ntask].ndim = 2;
      tasks[ntask].pdim = 0;
      tasks[ntask].qdim = 1;
      if (iface == ZLO) tasks[ntask].hi[2] = cells[icell].lo[2];
      else tasks[ntask].lo[2] = cells[icell].hi[2];
      if (iface == ZLO) tasks[ntask].normal[2] = 1.0;
      else tasks[ntask].normal[2] = -1.0;
    }

    // indot = dot product of vstream with inward face normal
    // skip task if indot < -3.0, to not allow any particles to be inserted
    // 0.0 up to -3.0 allows backflow influx of particles opposite to
    //   streaming velocity up to a reasonable limit

    indot = vstream[0]*tasks[ntask].normal[0] +
      vstream[1]*tasks[ntask].normal[1] +
      vstream[2]*tasks[ntask].normal[2];
    if (indot < -3.0) continue;

    // area = area for insertion
    // depends on dimension and axisymmetry

    if (iface == XLO || iface == XHI) {
      if (dimension == 3) 
        area = (cells[icell].hi[1]-cells[icell].lo[1]) * 
          (cells[icell].hi[2]-cells[icell].lo[2]);
      else if (domain->axisymmetric)
        area = (cells[icell].hi[1]*cells[icell].hi[1] -
                cells[icell].lo[1]*cells[icell].lo[1])*MY_PI;
      else area = cells[icell].hi[1]-cells[icell].lo[1];
    } else if (iface == YLO || iface == YHI) {
      if (dimension == 3)
        area = (cells[icell].hi[0]-cells[icell].lo[0]) * 
          (cells[icell].hi[2]-cells[icell].lo[2]);
      else if (domain->axisymmetric)
        area = 2.0*MY_PI*cells[icell].hi[1] *
          (cells[icell].hi[0]-cells[icell].lo[0]);
      else area = cells[icell].hi[0]-cells[icell].lo[0];
    } else if (iface == ZLO || iface == ZHI) {
      area = (cells[icell].hi[0]-cells[icell].lo[0]) * 
	(cells[icell].hi[1]-cells[icell].lo[1]);
    }      

    // set ntarget and ntargetsp via mol_inflow()

    tasks[ntask].ntarget = 0.0;
    for (isp = 0; isp < nspecies; isp++) {
      tasks[ntask].ntargetsp[isp] = 
        mol_inflow(indot,vscale[isp],fraction[isp]);
      tasks[ntask].ntargetsp[isp] *= nrho*area*dt / fnum;
      tasks[ntask].ntargetsp[isp] /= cinfo[icell].weight;
      tasks[ntask].ntarget += tasks[ntask].ntargetsp[isp];
    }

    // increment task counter

    ntask++;
  }

  // return # of tasks for this cell

  return ntask-ntaskorig;
}

/* ----------------------------------------------------------------------
   insert particles in grid cells with faces touching inflow boundaries
   NOTE:
     currently not allowing particle insertion on backflow boundaries
     enforced by indot >= 0.0 check in init()
     could allow particle insertion on backflow boundaries
       when streaming velocity is small enough
     need to insure double do-while loops below do not spin endlessly
------------------------------------------------------------------------- */

void FixEmitFace::perform_task()
{
  int pcell,ninsert,nactual,isp,ispecies,ndim,pdim,qdim,id;
  double *lo,*hi,*normal;
  double x[3],v[3];
  double indot,scosine,rn,ntarget,vr;
  double beta_un,normalized_distbn_fn,theta,erot,evib;
  Particle::OnePart *p;

  dt = update->dt;
  Grid::ChildCell *cells = grid->cells;
  int *species = particle->mixture[imix]->species;

  // insert particles for each task = cell/face pair
  // ntarget/ninsert is either perspecies or for all species
  // for one particle:
  //   x = random position on face
  //   v = randomized thermal velocity + vstream
  //       first stage: normal dimension (ndim)
  //       second stage: parallel dimensions (pdim,qdim)

  // double while loop until randomized particle velocity meets 2 criteria
  // inner do-while loop:
  //   v = vstream-component + vthermal is into simulation box
  //   see Bird 1994, p 425
  // outer do-while loop:
  //   shift Maxwellian distribution by stream velocity component
  //   see Bird 1994, p 259, eq 12.5

  int nfix_add_particle = modify->n_add_particle;

  for (int i = 0; i < ntask; i++) {
    pcell = tasks[i].pcell;
    ndim = tasks[i].ndim;
    pdim = tasks[i].pdim;
    qdim = tasks[i].qdim;
    lo = tasks[i].lo;
    hi = tasks[i].hi;
    normal = tasks[i].normal;

    indot = vstream[0]*normal[0] + vstream[1]*normal[1] + vstream[2]*normal[2];

    if (perspecies) {
      for (isp = 0; isp < nspecies; isp++) {
        ispecies = species[isp];
	ntarget = tasks[i].ntargetsp[isp]+random->uniform();
	ninsert = static_cast<int> (ntarget);
        scosine = indot / vscale[isp];

        nactual = 0;
	for (int m = 0; m < ninsert; m++) {
	  x[0] = lo[0] + random->uniform() * (hi[0]-lo[0]);
	  x[1] = lo[1] + random->uniform() * (hi[1]-lo[1]);
	  if (dimension == 3) x[2] = lo[2] + random->uniform() * (hi[2]-lo[2]);
          else x[2] = 0.0;

          if (region && !region->match(x)) continue;

	  do {
	    do beta_un = (6.0*random->gaussian() - 3.0);
	    while (beta_un + scosine < 0.0);
	    normalized_distbn_fn = 2.0 * (beta_un + scosine) / 
	      (scosine + sqrt(scosine*scosine + 2.0)) *
	      exp(0.5 + (0.5*scosine)*(scosine-sqrt(scosine*scosine + 2.0)) - 
		  beta_un*beta_un);
	  } while (normalized_distbn_fn < random->uniform());
	  
          v[ndim] = beta_un*vscale[isp]*normal[ndim] + vstream[ndim];

          theta = MY_2PI * random->gaussian();
          vr = vscale[isp] * sqrt(-log(random->uniform()));
          v[pdim] = vr * sin(theta) + vstream[pdim];
          v[qdim] = vr * cos(theta) + vstream[qdim];
          erot = particle->erot(ispecies,temp_thermal,random);
          evib = particle->evib(ispecies,temp_thermal,random);
          id = MAXSMALLINT*random->uniform();

	  particle->add_particle(id,ispecies,pcell,x,v,erot,evib);
          nactual++;

          p = &particle->particles[particle->nlocal-1];
          p->flag = PINSERT;
          p->dtremain = dt * random->uniform();

          if (nfix_add_particle) 
            modify->add_particle(particle->nlocal-1,temp_thermal,vstream);
	}

	nsingle += nactual;
      }

    } else {
      if (np == 0) {
	ntarget = tasks[i].ntarget+random->uniform();
	ninsert = static_cast<int> (ntarget);
      } else {
	ninsert = npertask;
	if (i >= nthresh) ninsert++;
      }

      nactual = 0;
      for (int m = 0; m < ninsert; m++) {
	rn = random->uniform();
	isp = 0;
	while (cummulative[isp] < rn) isp++;
        ispecies = species[isp];
        scosine = indot / vscale[isp];

	x[0] = lo[0] + random->uniform() * (hi[0]-lo[0]);
	x[1] = lo[1] + random->uniform() * (hi[1]-lo[1]);
        if (dimension == 3) x[2] = lo[2] + random->uniform() * (hi[2]-lo[2]);
        else x[2] = 0.0;

        if (region && !region->match(x)) continue;

	do {
	  do {
	    beta_un = (6.0*random->gaussian() - 3.0);
	  } while (beta_un + scosine < 0.0);
	  normalized_distbn_fn = 2.0 * (beta_un + scosine) / 
	    (scosine + sqrt(scosine*scosine + 2.0)) *
	    exp(0.5 + (0.5*scosine)*(scosine-sqrt(scosine*scosine + 2.0)) - 
		beta_un*beta_un);
	} while (normalized_distbn_fn < random->uniform());
	
        v[ndim] = beta_un*vscale[isp]*normal[ndim] + vstream[ndim];

        theta = MY_2PI * random->gaussian();
        vr = vscale[isp] * sqrt(-log(random->uniform()));
        v[pdim] = vr * sin(theta) + vstream[pdim];
        v[qdim] = vr * cos(theta) + vstream[qdim];
        erot = particle->erot(ispecies,temp_thermal,random);
        evib = particle->evib(ispecies,temp_thermal,random);
        id = MAXSMALLINT*random->uniform();

	particle->add_particle(id,ispecies,pcell,x,v,erot,evib);
        nactual++;

        p = &particle->particles[particle->nlocal-1];
        p->flag = PINSERT;
        p->dtremain = dt * random->uniform();

        if (nfix_add_particle) 
          modify->add_particle(particle->nlocal-1,temp_thermal,vstream);
      }

      nsingle += nactual;
    }
  }
}

/* ----------------------------------------------------------------------
   inserting into split cell icell on face iface
   determine which sub cell the face is part of
   face cannot be touched by surfs, so entire face is part of one sub cell
   compute which via update->split() and return it
------------------------------------------------------------------------- */

int FixEmitFace::split(int icell, int iface)
{
  double x[3];

  Grid::ChildCell *cells = grid->cells;

  // x = center point on face

  x[0] = 0.5 * (cells[icell].lo[0] + cells[icell].hi[0]);
  x[1] = 0.5 * (cells[icell].lo[1] + cells[icell].hi[1]);
  x[2] = 0.5 * (cells[icell].lo[2] + cells[icell].hi[2]);
  if (domain->dimension == 2) x[2] = 0.0;

  if (iface == XLO) x[0] = cells[icell].lo[0];
  else if (iface == XHI) x[0] = cells[icell].hi[0];
  else if (iface == YLO) x[1] = cells[icell].lo[1];
  else if (iface == YHI) x[1] = cells[icell].hi[1];
  else if (iface == ZLO) x[2] = cells[icell].lo[2];
  else if (iface == ZHI) x[2] = cells[icell].hi[2];

  int splitcell;
  if (dimension == 2) splitcell = update->split2d(icell,x);
  else splitcell = update->split3d(icell,x);
  return splitcell;
}

/* ----------------------------------------------------------------------
   pack one task into buf
   return # of bytes packed
   if not memflag, only return count, do not fill buf
------------------------------------------------------------------------- */

int FixEmitFace::pack_task(int itask, char *buf, int memflag)
{
  char *ptr = buf;
  if (memflag) memcpy(ptr,&tasks[itask],sizeof(Task));
  ptr += sizeof(Task);
  ptr = ROUNDUP(ptr);

  // pack task vector

  if (memflag) memcpy(ptr,tasks[itask].ntargetsp,nspecies*sizeof(double));
  ptr += nspecies*sizeof(double);

  return ptr-buf;
}

/* ----------------------------------------------------------------------
   unpack one task from buf
------------------------------------------------------------------------- */

int FixEmitFace::unpack_task(char *buf, int icell)
{
  char *ptr = buf;

  if (ntask == ntaskmax) grow_task();
  double *ntargetsp = tasks[ntask].ntargetsp;

  memcpy(&tasks[ntask],ptr,sizeof(Task));
  ptr += sizeof(Task);
  ptr = ROUNDUP(ptr);

  // unpack task vector

  memcpy(ntargetsp,ptr,nspecies*sizeof(double));
  ptr += nspecies*sizeof(double);

  tasks[ntask].ntargetsp = ntargetsp;

  // reset task icell and pcell
  // if a split cell, set pcell via split() which calls update->split()
  //   which will use current sub cells of icell

  tasks[ntask].icell = icell;
  if (grid->cells[icell].nsplit == 1) tasks[ntask].pcell = icell;
  else tasks[ntask].pcell = split(icell,tasks[ntask].iface);

  ntask++;
  return ptr-buf;
}

/* ----------------------------------------------------------------------
   copy N tasks starting at index oldfirst to index first
------------------------------------------------------------------------- */

void FixEmitFace::copy_task(int icell, int n, int first, int oldfirst)
{
  // reset icell in each copied task
  // copy task vector

  if (first == oldfirst) {
    for (int i = 0; i < n; i++) {
      tasks[first].icell = icell;
      first++;
    }

  } else {
    for (int i = 0; i < n; i++) {
      double *ntargetsp = tasks[first].ntargetsp;

      memcpy(&tasks[first],&tasks[oldfirst],sizeof(Task));
      memcpy(ntargetsp,tasks[oldfirst].ntargetsp,nspecies*sizeof(double));

      tasks[first].ntargetsp = ntargetsp;

      tasks[first].icell = icell;
      first++;
      oldfirst++;
    }
  }

  ntask += n;
}

/* ----------------------------------------------------------------------
   grow task list
------------------------------------------------------------------------- */

void FixEmitFace::grow_task()
{
  int oldmax = ntaskmax;
  ntaskmax += DELTATASK;
  tasks = (Task *) memory->srealloc(tasks,ntaskmax*sizeof(Task),
				    "emit/face:tasks");

  // set all new task bytes to 0 so valgrind won't complain
  // if bytes between fields are uninitialized

  memset(&tasks[oldmax],0,(ntaskmax-oldmax)*sizeof(Task));

  // allocate vector in each new task

  for (int i = oldmax; i < ntaskmax; i++)
    tasks[i].ntargetsp = new double[nspecies];
}

/* ----------------------------------------------------------------------
   reset pcell for all compress task entries
   called from Grid::compress() after grid cells have been compressed
   wait to do this until now b/c split cells accessed by split()
     are setup in Grid::compress() between compress_grid() 
     and post_compress_grid()
------------------------------------------------------------------------- */

void FixEmitFace::post_compress_grid()
{
  Grid::ChildCell *cells = grid->cells;

  for (int i = 0; i < ntask; i++) {
    int icell = tasks[i].icell;
    if (cells[icell].nsplit == 1) tasks[i].pcell = icell;
    else tasks[i].pcell = split(icell,tasks[i].iface);
  }
}

/* ----------------------------------------------------------------------
   process keywords specific to this class
------------------------------------------------------------------------- */

int FixEmitFace::option(int narg, char **arg)
{
  if (strcmp(arg[0],"n") == 0) {
    if (2 > narg) error->all(FLERR,"Illegal fix emit/face command");
    np = atoi(arg[1]);
    if (np <= 0) error->all(FLERR,"Illegal fix emit/face command");
    return 2;
  }

  error->all(FLERR,"Illegal fix emit/face command");
  return 0;
}
