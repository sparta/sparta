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
#include "input.h"
#include "random_knuth.h"
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
enum{NOSUBSONIC,PTBOTH,PONLY};

#define DELTATASK 256
#define TEMPLIMIT 1.0e5

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
  subsonic = 0;
  subsonic_style = NOSUBSONIC;
  subsonic_warning = 0;
  twopass = 0;

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
  if (np > 0 && subsonic)
    error->all(FLERR,"Cannot use fix emit/face n > 0 with subsonic");

  // task list and subsonic data structs

  tasks = NULL;
  ntask = ntaskmax = 0;

  maxactive = 0;
  activecell = NULL;
}

/* ---------------------------------------------------------------------- */

FixEmitFace::~FixEmitFace()
{
  if (copymode) return;

  if (tasks) {
    for (int i = 0; i < ntaskmax; i++) {
      delete [] tasks[i].ntargetsp;
      delete [] tasks[i].vscale;
    }
    memory->sfree(tasks);
  }
  memory->destroy(activecell);
}

/* ---------------------------------------------------------------------- */

void FixEmitFace::init()
{
  // invoke FixEmit::init() to set flags

  FixEmit::init();

  // copies of class data before invoking parent init() and count_task()

  dimension = domain->dimension;
  fnum = update->fnum;
  dt = update->dt;

  nspecies = particle->mixture[imix]->nspecies;
  fraction = particle->mixture[imix]->fraction;
  cummulative = particle->mixture[imix]->cummulative;

  // subsonic prefactor

  tprefactor = update->mvv2e / (3.0*update->boltz);

  // mixture soundspeed, used by subsonic PONLY as default cell property

  double avegamma = 0.0;
  double avemass = 0.0;

  for (int m = 0; m < nspecies; m++) {
    int ispecies = particle->mixture[imix]->species[m];
    avemass += fraction[m] * particle->species[ispecies].mass;
    avegamma += fraction[m] * (1.0 + 2.0 /
                               (3.0 + particle->species[ispecies].rotdof));
  }

  soundspeed_mixture = sqrt(avegamma * update->boltz *
                            particle->mixture[imix]->temp_thermal / avemass);

  // cannot inflow thru periodic boundary

  for (int i = 0; i < 6; i++)
    if (faces[i] && domain->bflag[i] == PERIODIC)
      error->all(FLERR,"Cannot use fix emit/face on periodic boundary");

  // cannot have inflow on yhi if axisymmetric

  double *vstream = particle->mixture[imix]->vstream;

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

  // if used, reallocate ntargetsp and vscale for each task
  // b/c nspecies count of mixture may have changed

  realloc_nspecies();

  // create tasks for all grid cells

  grid_changed();
}

/* ----------------------------------------------------------------------
   grid changed operation
   invoke create_tasks() to rebuild entire task list
   invoked after per-processor list of grid cells has changed
------------------------------------------------------------------------- */

void FixEmitFace::grid_changed()
{
  create_tasks();

  // if Np > 0, nper = # of insertions per task
  // set nthresh so as to achieve exactly Np insertions
  // tasks > tasks_with_no_extra need to insert 1 extra particle
  // NOTE: currently setting same # of insertions per task
  //       could instead weight by cell face area

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

/* ----------------------------------------------------------------------
   create tasks for one grid cell
   add them to tasks list and increment ntasks
------------------------------------------------------------------------- */

void FixEmitFace::create_task(int icell)
{
  int i,j,n,iface,flag,isp,extflag;
  int *cflags;
  double indot,area,ntargetsp;

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  double nrho = particle->mixture[imix]->nrho;
  double *vstream = particle->mixture[imix]->vstream;
  double *vscale = particle->mixture[imix]->vscale;

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
                line_touch_quad_face(lines[n].p1,lines[n].p2,
                                     iface,cells[icell].lo,cells[icell].hi)) {
              flag = 0;
              break;
            }
          }
        } else if (flag && dimension == 3) {
          for (j = 0; j < cells[icell].nsurf; j++) {
            n = cells[icell].csurfs[j];
            if (Geometry::
                tri_touch_hex_face(tris[n].p1,tris[n].p2,tris[n].p3,
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

    indot = vstream[0]*tasks[ntask].normal[0] +
      vstream[1]*tasks[ntask].normal[1] +
      vstream[2]*tasks[ntask].normal[2];

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
    tasks[ntask].area = area;

    // set ntarget and ntargetsp via mol_inflow()
    // skip task if final ntarget = 0.0, due to large outbound vstream
    // do not skip for subsonic since it resets ntarget every step

    tasks[ntask].ntarget = 0.0;
    for (isp = 0; isp < nspecies; isp++) {
      ntargetsp = mol_inflow(indot,vscale[isp],fraction[isp]);
      ntargetsp *= nrho*area*dt / fnum;
      ntargetsp /= cinfo[icell].weight;
      tasks[ntask].ntarget += ntargetsp;
      if (perspecies) tasks[ntask].ntargetsp[isp] = ntargetsp;
    }

    if (!subsonic) {
      if (tasks[ntask].ntarget == 0.0) continue;
      if (tasks[ntask].ntarget >= MAXSMALLINT)
        error->one(FLERR,
                   "Fix emit/face insertion count exceeds 32-bit int");
    }

    // initialize other task values with mixture properties
    // may be overwritten by subsonic methods

    tasks[ntask].nrho = particle->mixture[imix]->nrho;
    tasks[ntask].temp_thermal = particle->mixture[imix]->temp_thermal;
    tasks[ntask].temp_rot = particle->mixture[imix]->temp_rot;
    tasks[ntask].temp_vib = particle->mixture[imix]->temp_vib;
    tasks[ntask].vstream[0] = particle->mixture[imix]->vstream[0];
    tasks[ntask].vstream[1] = particle->mixture[imix]->vstream[1];
    tasks[ntask].vstream[2] = particle->mixture[imix]->vstream[2];

    // increment task counter

    ntask++;
  }
}

/* ----------------------------------------------------------------------
   insert particles in grid cells with faces touching inflow boundaries
------------------------------------------------------------------------- */

void FixEmitFace::perform_task()
{
  if (!twopass) perform_task_onepass();
  else perform_task_twopass();
}

/* ----------------------------------------------------------------------
   perform insertion in one pass thru tasks
   this is simpler, somewhat faster code
   but uses random #s differently than Kokkos, so insertions are different
------------------------------------------------------------------------- */

void FixEmitFace::perform_task_onepass()
{
  int pcell,ninsert,nactual,isp,ispecies,ndim,pdim,qdim,id;
  double indot,scosine,rn,ntarget,vr;
  double beta_un,normalized_distbn_fn,theta,erot,evib;
  double temp_thermal,temp_rot,temp_vib;
  double x[3],v[3];
  double *lo,*hi,*normal,*vstream,*vscale;
  Particle::OnePart *p;

  dt = update->dt;
  int *species = particle->mixture[imix]->species;

  // if subsonic, re-compute particle inflow counts for each task
  // also computes current per-task temp_thermal and vstream

  if (subsonic) subsonic_inflow();

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

  int nfix_update_custom = modify->n_update_custom;

  for (int i = 0; i < ntask; i++) {
    pcell = tasks[i].pcell;
    ndim = tasks[i].ndim;
    pdim = tasks[i].pdim;
    qdim = tasks[i].qdim;
    lo = tasks[i].lo;
    hi = tasks[i].hi;
    normal = tasks[i].normal;

    temp_thermal = tasks[i].temp_thermal;
    temp_rot = tasks[i].temp_rot;
    temp_vib = tasks[i].temp_vib;
    vstream = tasks[i].vstream;

    if (subsonic_style == PONLY) vscale = tasks[i].vscale;
    else vscale = particle->mixture[imix]->vscale;

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
            do beta_un = (6.0*random->uniform() - 3.0);
            while (beta_un + scosine < 0.0);
            normalized_distbn_fn = 2.0 * (beta_un + scosine) /
              (scosine + sqrt(scosine*scosine + 2.0)) *
              exp(0.5 + (0.5*scosine)*(scosine-sqrt(scosine*scosine + 2.0)) -
                  beta_un*beta_un);
          } while (normalized_distbn_fn < random->uniform());

          v[ndim] = beta_un*vscale[isp]*normal[ndim] + vstream[ndim];

          theta = MY_2PI * random->uniform();
          vr = vscale[isp] * sqrt(-log(random->uniform()));
          v[pdim] = vr * sin(theta) + vstream[pdim];
          v[qdim] = vr * cos(theta) + vstream[qdim];
          erot = particle->erot(ispecies,temp_rot,random);
          evib = particle->evib(ispecies,temp_vib,random);
          id = MAXSMALLINT*random->uniform();

          particle->add_particle(id,ispecies,pcell,x,v,erot,evib);
          nactual++;

          p = &particle->particles[particle->nlocal-1];
          p->flag = PINSERT;
          p->dtremain = dt * random->uniform();

          if (nfix_update_custom)
            modify->update_custom(particle->nlocal-1,temp_thermal,
                                 temp_rot,temp_vib,vstream);
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
            beta_un = (6.0*random->uniform() - 3.0);
          } while (beta_un + scosine < 0.0);
          normalized_distbn_fn = 2.0 * (beta_un + scosine) /
            (scosine + sqrt(scosine*scosine + 2.0)) *
            exp(0.5 + (0.5*scosine)*(scosine-sqrt(scosine*scosine + 2.0)) -
                beta_un*beta_un);
        } while (normalized_distbn_fn < random->uniform());

        v[ndim] = beta_un*vscale[isp]*normal[ndim] + vstream[ndim];

        theta = MY_2PI * random->uniform();
        vr = vscale[isp] * sqrt(-log(random->uniform()));
        v[pdim] = vr * sin(theta) + vstream[pdim];
        v[qdim] = vr * cos(theta) + vstream[qdim];
        erot = particle->erot(ispecies,temp_rot,random);
        evib = particle->evib(ispecies,temp_vib,random);
        id = MAXSMALLINT*random->uniform();

        particle->add_particle(id,ispecies,pcell,x,v,erot,evib);
        nactual++;

        p = &particle->particles[particle->nlocal-1];
        p->flag = PINSERT;
        p->dtremain = dt * random->uniform();

        if (nfix_update_custom)
          modify->update_custom(particle->nlocal-1,temp_thermal,
                               temp_rot,temp_vib,vstream);
      }

      nsingle += nactual;
    }
  }
}

/* ----------------------------------------------------------------------
   perform insertion the way Kokkos does in two passes thru tasks
   this uses random #s the same as Kokkos, for easier debugging
------------------------------------------------------------------------- */

void FixEmitFace::perform_task_twopass()
{
  int pcell,ninsert,nactual,isp,ispecies,ndim,pdim,qdim,id;
  double indot,scosine,rn,ntarget,vr;
  double beta_un,normalized_distbn_fn,theta,erot,evib;
  double temp_thermal,temp_rot,temp_vib;
  double x[3],v[3];
  double *lo,*hi,*normal,*vstream,*vscale;
  Particle::OnePart *p;

  dt = update->dt;
  int *species = particle->mixture[imix]->species;

  // if subsonic, re-compute particle inflow counts for each task
  // also computes current per-task temp_thermal and vstream

  if (subsonic) subsonic_inflow();

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

  int nfix_update_custom = modify->n_update_custom;

  int ninsert_dim1 = perspecies ? nspecies : 1;
  int** ninsert_values;
  memory->create(ninsert_values, ntask, ninsert_dim1, "fix_emit_face:ninsert");

  for (int i = 0; i < ntask; i++) {
    if (perspecies) {
      for (isp = 0; isp < nspecies; isp++) {
        ntarget = tasks[i].ntargetsp[isp]+random->uniform();
        ninsert = static_cast<int> (ntarget);
        ninsert_values[i][isp] = ninsert;
      }
    } else {
      if (np == 0) {
        ntarget = tasks[i].ntarget+random->uniform();
        ninsert = static_cast<int> (ntarget);
      } else {
        ninsert = npertask;
        if (i >= nthresh) ninsert++;
      }
      ninsert_values[i][0] = ninsert;
    }
  }

  for (int i = 0; i < ntask; i++) {
    pcell = tasks[i].pcell;
    ndim = tasks[i].ndim;
    pdim = tasks[i].pdim;
    qdim = tasks[i].qdim;
    lo = tasks[i].lo;
    hi = tasks[i].hi;
    normal = tasks[i].normal;

    temp_thermal = tasks[i].temp_thermal;
    temp_rot = tasks[i].temp_rot;
    temp_vib = tasks[i].temp_vib;
    vstream = tasks[i].vstream;

    if (subsonic_style == PONLY) vscale = tasks[i].vscale;
    else vscale = particle->mixture[imix]->vscale;

    indot = vstream[0]*normal[0] + vstream[1]*normal[1] + vstream[2]*normal[2];

    if (perspecies) {
      for (isp = 0; isp < nspecies; isp++) {
        ispecies = species[isp];
        ninsert = ninsert_values[i][isp];
        scosine = indot / vscale[isp];

        nactual = 0;
        for (int m = 0; m < ninsert; m++) {
          x[0] = lo[0] + random->uniform() * (hi[0]-lo[0]);
          x[1] = lo[1] + random->uniform() * (hi[1]-lo[1]);
          if (dimension == 3) x[2] = lo[2] + random->uniform() * (hi[2]-lo[2]);
          else x[2] = 0.0;

          if (region && !region->match(x)) continue;

          do {
            do beta_un = (6.0*random->uniform() - 3.0);
            while (beta_un + scosine < 0.0);
            normalized_distbn_fn = 2.0 * (beta_un + scosine) /
              (scosine + sqrt(scosine*scosine + 2.0)) *
              exp(0.5 + (0.5*scosine)*(scosine-sqrt(scosine*scosine + 2.0)) -
                  beta_un*beta_un);
          } while (normalized_distbn_fn < random->uniform());

          v[ndim] = beta_un*vscale[isp]*normal[ndim] + vstream[ndim];

          theta = MY_2PI * random->uniform();
          vr = vscale[isp] * sqrt(-log(random->uniform()));
          v[pdim] = vr * sin(theta) + vstream[pdim];
          v[qdim] = vr * cos(theta) + vstream[qdim];
          erot = particle->erot(ispecies,temp_rot,random);
          evib = particle->evib(ispecies,temp_vib,random);
          id = MAXSMALLINT*random->uniform();

          particle->add_particle(id,ispecies,pcell,x,v,erot,evib);
          nactual++;

          p = &particle->particles[particle->nlocal-1];
          p->flag = PINSERT;
          p->dtremain = dt * random->uniform();

          if (nfix_update_custom)
            modify->update_custom(particle->nlocal-1,temp_thermal,
                temp_rot,temp_vib,vstream);
        }

        nsingle += nactual;
      }

    } else {
      ninsert = ninsert_values[i][0];

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
            beta_un = (6.0*random->uniform() - 3.0);
          } while (beta_un + scosine < 0.0);
          normalized_distbn_fn = 2.0 * (beta_un + scosine) /
            (scosine + sqrt(scosine*scosine + 2.0)) *
            exp(0.5 + (0.5*scosine)*(scosine-sqrt(scosine*scosine + 2.0)) -
                beta_un*beta_un);
        } while (normalized_distbn_fn < random->uniform());

        v[ndim] = beta_un*vscale[isp]*normal[ndim] + vstream[ndim];

        theta = MY_2PI * random->uniform();
        vr = vscale[isp] * sqrt(-log(random->uniform()));
        v[pdim] = vr * sin(theta) + vstream[pdim];
        v[qdim] = vr * cos(theta) + vstream[qdim];
        erot = particle->erot(ispecies,temp_rot,random);
        evib = particle->evib(ispecies,temp_vib,random);
        id = MAXSMALLINT*random->uniform();

        particle->add_particle(id,ispecies,pcell,x,v,erot,evib);
        nactual++;

        p = &particle->particles[particle->nlocal-1];
        p->flag = PINSERT;
        p->dtremain = dt * random->uniform();

        if (nfix_update_custom)
          modify->update_custom(particle->nlocal-1,temp_thermal,
              temp_rot,temp_vib,vstream);
      }

      nsingle += nactual;
    }
  }

  memory->destroy(ninsert_values);
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
   recalculate task properties based on subsonic BC
------------------------------------------------------------------------- */

void FixEmitFace::subsonic_inflow()
{
  // for grid cells that are part of tasks:
  // calculate local nrho, vstream, and thermal temperature
  // if needed sort particles for grid cells with tasks

  if (!particle->sorted) subsonic_sort();
  subsonic_grid();

  // recalculate particle insertion counts for each task
  // recompute mixture vscale, since depends on temp_thermal

  int isp,icell;
  double mass,indot,area,nrho,temp_thermal,vscale,ntargetsp;
  double *vstream,*normal;

  Particle::Species *species = particle->species;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int *mspecies = particle->mixture[imix]->species;
  double fnum = update->fnum;
  double boltz = update->boltz;

  for (int i = 0; i < ntask; i++) {
    vstream = tasks[i].vstream;
    normal = tasks[i].normal;
    indot = vstream[0]*normal[0] + vstream[1]*normal[1] + vstream[2]*normal[2];

    area = tasks[i].area;
    nrho = tasks[i].nrho;
    temp_thermal = tasks[i].temp_thermal;
    icell = tasks[i].icell;

    tasks[i].ntarget = 0.0;
    for (isp = 0; isp < nspecies; isp++) {
      mass = species[mspecies[isp]].mass;
      vscale = sqrt(2.0 * boltz * temp_thermal / mass);
      ntargetsp = mol_inflow(indot,vscale,fraction[isp]);
      ntargetsp *= nrho*area*dt / fnum;
      ntargetsp /= cinfo[icell].weight;
      tasks[i].ntarget += ntargetsp;
      if (perspecies) tasks[i].ntargetsp[isp] = ntargetsp;
    }
    if (tasks[i].ntarget >= MAXSMALLINT)
      error->one(FLERR,
                 "Fix emit/face subsonic insertion count exceeds 32-bit int");
  }
}

/* ----------------------------------------------------------------------
   identify particles in grid cells associated with a task
   store count and linked list, same as for particle sorting
------------------------------------------------------------------------- */

void FixEmitFace::subsonic_sort()
{
  int i,icell;

  // initialize particle sort lists for grid cells assigned to tasks
  // use task pcell, not icell

  Grid::ChildInfo *cinfo = grid->cinfo;

  for (i = 0; i < ntask; i++) {
    icell = tasks[i].pcell;
    cinfo[icell].first = -1;
    cinfo[icell].count = 0;
  }

  // reallocate particle next list if necessary

  particle->sort_allocate();

  // update list of active grid cells if necessary
  // active cells = those assigned to tasks
  // active_current flag set by parent class

  if (!active_current) {
    if (grid->nlocal > maxactive) {
      memory->destroy(activecell);
      maxactive = grid->nlocal;
      memory->create(activecell,maxactive,"emit/face:active");
    }
    memset(activecell,0,maxactive*sizeof(int));
    for (i = 0; i < ntask; i++) activecell[tasks[i].pcell] = 1;
    active_current = 1;
  }

  // loop over particles to store linked lists for active cells
  // not using reverse loop like Particle::sort(),
  //   since this should only be created/used occasionally

  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;
  int nlocal = particle->nlocal;

  for (i = 0; i < nlocal; i++) {
    icell = particles[i].icell;
    if (!activecell[icell]) continue;
    next[i] = cinfo[icell].first;
    cinfo[icell].first = i;
    cinfo[icell].count++;
  }
}

/* ----------------------------------------------------------------------
   compute number density, thermal temperature, stream velocity
   only for grid cells associated with a task
   first compute for grid cells, then adjust due to boundary conditions
------------------------------------------------------------------------- */

void FixEmitFace::subsonic_grid()
{
  int m,ip,np,icell,ispecies,ndim;
  double mass,masstot,gamma,ke,sign;
  double nrho_cell,massrho_cell,temp_thermal_cell,press_cell;
  double mass_cell,gamma_cell,soundspeed_cell;
  double mv[4];
  double *v,*vstream,*vscale;

  Grid::ChildInfo *cinfo = grid->cinfo;
  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;
  Particle::Species *species = particle->species;
  double boltz = update->boltz;

  int temp_exceed_flag = 0;
  double tempmax = 0.0;

  for (int i = 0; i < ntask; i++) {
    icell = tasks[i].pcell;
    np = cinfo[icell].count;

    // accumulate needed per-particle quantities
    // mv = mass*velocity terms, masstot = total mass
    // gamma = rotational/tranlational DOFs

    mv[0] = mv[1] = mv[2] = mv[3] = 0.0;
    masstot = gamma = 0.0;

    ip = cinfo[icell].first;
    while (ip >= 0) {
      ispecies = particles[ip].ispecies;
      mass = species[ispecies].mass;
      v = particles[ip].v;
      mv[0] += mass*v[0];
      mv[1] += mass*v[1];
      mv[2] += mass*v[2];
      mv[3] += mass * (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
      masstot += mass;
      gamma += 1.0 + 2.0 / (3.0 + species[ispecies].rotdof);
      ip = next[ip];
    }

    // compute/store nrho, 3 temps, vstream for task
    // also vscale for PONLY
    // if sound speed = 0.0 due to <= 1 particle in cell or
    //   all particles having COM velocity, set via mixture properties

    vstream = tasks[i].vstream;
    if (np) {
      vstream[0] = mv[0] / masstot;
      vstream[1] = mv[1] / masstot;
      vstream[2] = mv[2] / masstot;
    } else vstream[0] = vstream[1] = vstream[2] = 0.0;

    if (subsonic_style == PTBOTH) {
      tasks[i].nrho = nsubsonic;
      temp_thermal_cell = tsubsonic;

    } else {
      nrho_cell = np * fnum / cinfo[icell].volume;
      massrho_cell = masstot * fnum / cinfo[icell].volume;
      if (np > 1) {
        ke = mv[3]/np - (mv[0]*mv[0] + mv[1]*mv[1] + mv[2]*mv[2])/np/masstot;
        temp_thermal_cell = tprefactor * ke;
      } else temp_thermal_cell = particle->mixture[imix]->temp_thermal;

      press_cell = nrho_cell * boltz * temp_thermal_cell;
      if (np) {
        mass_cell = masstot / np;
        gamma_cell = gamma / np;
        soundspeed_cell = sqrt(gamma_cell*boltz*temp_thermal_cell / mass_cell);
      } else soundspeed_cell = soundspeed_mixture;

      tasks[i].nrho = nrho_cell +
        (psubsonic - press_cell) / (soundspeed_cell*soundspeed_cell);
      temp_thermal_cell = psubsonic / (boltz * tasks[i].nrho);
      if (temp_thermal_cell > TEMPLIMIT) {
        temp_exceed_flag = 1;
        tempmax = MAX(tempmax,temp_thermal_cell);
      }

      if (np)  {
        ndim = tasks[i].ndim;
        sign = tasks[i].normal[ndim];
        vstream[ndim] += sign *
          (psubsonic - press_cell) / (massrho_cell*soundspeed_cell);
      }

      vscale = tasks[i].vscale;
      for (m = 0; m < nspecies; m++) {
        ispecies = particle->mixture[imix]->species[m];
        vscale[m] = sqrt(2.0 * update->boltz * temp_thermal_cell /
                         species[ispecies].mass);
      }
    }

    tasks[i].temp_thermal = temp_thermal_cell;
    tasks[i].temp_rot = tasks[i].temp_vib = temp_thermal_cell;
  }

  // test if any task has invalid thermal temperature for first time

  if (!subsonic_warning)
    subsonic_warning = subsonic_temperature_check(temp_exceed_flag,tempmax);
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

  // allocate vectors in each new task or set to NULL

  if (perspecies) {
    for (int i = oldmax; i < ntaskmax; i++)
      tasks[i].ntargetsp = new double[nspecies];
  } else {
    for (int i = oldmax; i < ntaskmax; i++)
      tasks[i].ntargetsp = NULL;
  }

  if (subsonic_style == PONLY) {
    for (int i = oldmax; i < ntaskmax; i++)
      tasks[i].vscale = new double[nspecies];
  } else {
    for (int i = oldmax; i < ntaskmax; i++)
      tasks[i].vscale = NULL;
  }
}

/* ----------------------------------------------------------------------
   reallocate nspecies arrays
------------------------------------------------------------------------- */

void FixEmitFace::realloc_nspecies()
{
  if (perspecies) {
    for (int i = 0; i < ntask; i++) {
      delete [] tasks[i].ntargetsp;
      tasks[i].ntargetsp = new double[nspecies];
    }
  }
  if (subsonic_style == PONLY) {
    for (int i = 0; i < ntask; i++) {
      delete [] tasks[i].vscale;
      tasks[i].vscale = new double[nspecies];
    }
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

  if (strcmp(arg[0],"subsonic") == 0) {
    if (3 > narg) error->all(FLERR,"Illegal fix emit/face command");
    subsonic = 1;
    subsonic_style = PTBOTH;
    psubsonic = input->numeric(FLERR,arg[1]);
    if (psubsonic < 0.0) error->all(FLERR,"Illegal fix emit/face command");
    if (strcmp(arg[2],"NULL") == 0) subsonic_style = PONLY;
    else {
      tsubsonic = input->numeric(FLERR,arg[2]);
      if (tsubsonic <= 0.0)
        error->all(FLERR,"Subsonic temperature cannot be <= 0.0");
      nsubsonic = psubsonic / (update->boltz * tsubsonic);
    }
    return 3;
  }

  if (strcmp(arg[0],"twopass") == 0) {
    twopass = 1;
    return 1;
  }

  error->all(FLERR,"Illegal fix emit/face command");
  return 0;
}
