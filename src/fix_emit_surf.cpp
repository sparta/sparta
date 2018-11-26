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
#include "fix_emit_surf.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "particle.h"
#include "mixture.h"
#include "surf.h"
#include "modify.h"
#include "cut2d.h"
#include "cut3d.h"
#include "input.h"
#include "comm.h"
#include "random_park.h"
#include "math_extra.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};   // several files
enum{NOSUBSONIC,PTBOTH,PONLY};

#define DELTATASK 256
#define TEMPLIMIT 1.0e5

/* ---------------------------------------------------------------------- */

FixEmitSurf::FixEmitSurf(SPARTA *sparta, int narg, char **arg) :
  FixEmit(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix emit/surf command");
  
  imix = particle->find_mixture(arg[2]);
  if (imix < 0) 
    error->all(FLERR,"Fix emit/surf mixture ID does not exist");
  
  int igroup = surf->find_group(arg[3]);
  if (igroup < 0) 
    error->all(FLERR,"Fix emit/surf group ID does not exist");
  groupbit = surf->bitmask[igroup];
  
  // optional args
  
  np = 0;
  normalflag = 0;
  subsonic = 0;
  subsonic_style = NOSUBSONIC;
  subsonic_warning = 0;

  int iarg = 4;
  options(narg-iarg,&arg[iarg]);
  
  // error checks
  
  if (!surf->exist) 
    error->all(FLERR,"Fix emit/surf requires surface elements");
  if (np > 0 && perspecies) 
    error->all(FLERR,"Cannot use fix emit/face n > 0 with perspecies yes");
  
  // task list and subsonic data structs
  
  tasks = NULL;
  ntask = ntaskmax = 0;
  
  maxactive = 0;
  activecell = NULL;
}

/* ---------------------------------------------------------------------- */

FixEmitSurf::~FixEmitSurf()
{
  for (int i = 0; i < ntaskmax; i++) {
    delete [] tasks[i].ntargetsp;
    delete [] tasks[i].vscale;
    delete [] tasks[i].path;
    delete [] tasks[i].fracarea;
  }
  memory->sfree(tasks);
  memory->destroy(activecell);
}

/* ---------------------------------------------------------------------- */

void FixEmitSurf::init()
{
  // copies of class data before invoking parent init() and count_task()
  
  dimension = domain->dimension;
  fnum = update->fnum;
  dt = update->dt;
  
  nspecies = particle->mixture[imix]->nspecies;
  fraction = particle->mixture[imix]->fraction;
  cummulative = particle->mixture[imix]->cummulative;
  
  lines = surf->lines;
  tris = surf->tris;
  
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
  
  // create instance of Cut2d,Cut3d for geometry calculations
  
  if (dimension == 3) cut3d = new Cut3d(sparta);
  else cut2d = new Cut2d(sparta,domain->axisymmetric);
  
  // magvstream = magnitude of mxiture vstream vector
  // norm_vstream = unit vector in stream direction
  
  double *vstream = particle->mixture[imix]->vstream;
  magvstream = MathExtra::len3(vstream);
  
  norm_vstream[0] = vstream[0];
  norm_vstream[1] = vstream[1];
  norm_vstream[2] = vstream[2];
  if (norm_vstream[0] != 0.0 || norm_vstream[1] != 0.0 || 
      norm_vstream[2] != 0.0)
    MathExtra::norm3(norm_vstream);
  
  // if used, reallocate ntargetsp and vscale for each task
  // b/c nspecies count of mixture may have changed
  
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
  
  // invoke FixEmit::init() to populate task list
  // it calls create_task() for each grid cell
  
  ntask = 0;
  FixEmit::init();
  
  // if Np > 0, nper = # of insertions per task
  // set nthresh so as to achieve exactly Np insertions
  // tasks > tasks_with_no_extra need to insert 1 extra particle
  // NOTE: setting same # of insertions per task
  //       should weight by overlap area of cell/surf
  
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

  // deallocate Cut2d,Cut3d
  
  if (dimension == 3) delete cut3d;
  else delete cut2d;
}

/* ---------------------------------------------------------------------- */

void FixEmitSurf::setup()
{
  // needed for Kokkos because pointers are changed in UpdateKokkos::setup()

  lines = surf->lines;
  tris = surf->tris;
}

/* ---------------------------------------------------------------------- */

int FixEmitSurf::create_task(int icell)
{
  int i,m,isurf,isp,npoint,isplit,subcell;
  double indot,area,areaone,ntargetsp;
  double *normal,*p1,*p2,*p3,*path;
  double cpath[36],delta[3],e1[3],e2[3];
  
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Grid::SplitInfo *sinfo = grid->sinfo;
  
  double nrho = particle->mixture[imix]->nrho;
  double *vstream = particle->mixture[imix]->vstream;
  double *vscale = particle->mixture[imix]->vscale;
  
  // no tasks if no surfs in cell
  
  if (cells[icell].nsurf == 0) return 0;
  
  // loop over surfs in cell
  // use Cut2d/Cut3d to find overlap area and geoemtry of overlap
  
  int ntaskorig = ntask;
  
  double *lo = cells[icell].lo;
  double *hi = cells[icell].hi;
  int *csurfs = cells[icell].csurfs;
  int nsurf = cells[icell].nsurf;
  
  for (i = 0; i < nsurf; i++) {
    isurf = csurfs[i];
    
    if (dimension == 2) {
      if (!(lines[isurf].mask & groupbit)) continue;
    } else {
      if (!(tris[isurf].mask & groupbit)) continue;
    }
    
    // set cell parameters of task
    // pcell = sub cell for particles if a split cell
    
    if (ntask == ntaskmax) grow_task();
    
    tasks[ntask].icell = icell;
    tasks[ntask].isurf = isurf;
    if (cells[icell].nsplit == 1) tasks[ntask].pcell = icell;
    else {
      isplit = cells[icell].isplit;
      subcell = sinfo[isplit].csplits[i];
      tasks[ntask].pcell = sinfo[isplit].csubs[subcell];
    }
    
    // set geometry-dependent params of task
    // indot = vstream magnitude for normalflag = 1
    // indot = vstream dotted with surface normal for normalflag = 0
    // area = area for insertion = extent of line/triangle inside grid cell
    
    if (dimension == 2) {
      normal = lines[isurf].norm;
      if (normalflag) indot = magvstream;
      else indot = vstream[0]*normal[0] + vstream[1]*normal[1];
      
      p1 = lines[isurf].p1;
      p2 = lines[isurf].p2;
      npoint = cut2d->clip_external(p1,p2,lo,hi,cpath);
      if (npoint < 2) continue;
      
      tasks[ntask].npoint = 2;
      delete [] tasks[ntask].path;
      tasks[ntask].path = new double[6];
      path = tasks[ntask].path;
      path[0] = cpath[0];
      path[1] = cpath[1];
      path[2] = 0.0;
      path[3] = cpath[2];
      path[4] = cpath[3];
      path[5] = 0.0;
      
      // axisymmetric "area" of line segment = surf area of truncated cone
      // PI (y1+y2) sqrt( (y1-y2)^2 + (x1-x2)^2) )
      
      if (domain->axisymmetric) {
        double sqrtarg = (path[1]-path[4])*(path[1]-path[4]) + 
          (path[0]-path[3])*(path[0]-path[3]);
        area = MY_PI * (path[1]+path[4]) * sqrt(sqrtarg);
      } else {
        MathExtra::sub3(&path[0],&path[3],delta);
        area = MathExtra::len3(delta);
      }
      tasks[ntask].area = area;
      
      // set 2 tangent vectors to surf normal
      // tan1 is in xy plane, 90 degrees from normal
      // tan2 is unit +z vector
      
      tasks[ntask].tan1[0] = normal[1];
      tasks[ntask].tan1[1] = -normal[0];
      tasks[ntask].tan1[2] = 0.0;
      tasks[ntask].tan2[0] = 0.0;
      tasks[ntask].tan2[1] = 0.0;
      tasks[ntask].tan2[2] = 1.0;
      
    } else {
      normal = tris[isurf].norm;
      if (normalflag) indot = magvstream;
      else indot = vstream[0]*normal[0] + vstream[1]*normal[1] + 
             vstream[2]*normal[2];
      
      p1 = tris[isurf].p1;
      p2 = tris[isurf].p2;
      p3 = tris[isurf].p3;
      npoint = cut3d->clip_external(p1,p2,p3,lo,hi,cpath);
      if (npoint < 3) continue;
      
      tasks[ntask].npoint = npoint;
      delete [] tasks[ntask].path;
      tasks[ntask].path = new double[npoint*3];
      path = tasks[ntask].path;
      memcpy(path,cpath,npoint*3*sizeof(double));
      delete [] tasks[ntask].fracarea;
      tasks[ntask].fracarea = new double[npoint-2];
      
      area = 0.0;
      p1 = &path[0];
      for (m = 0; m < npoint-2; m++) {
        p2 = &path[3*(m+1)];
        p3 = &path[3*(m+2)];
        MathExtra::sub3(p2,p1,e1);
        MathExtra::sub3(p3,p1,e2);
        MathExtra::cross3(e1,e2,delta);
        areaone = fabs(0.5*MathExtra::len3(delta));
        area += areaone;
        tasks[ntask].fracarea[m] = area;
      }
      tasks[ntask].area = area;
      for (m = 0; m < npoint-2; m++)
        tasks[ntask].fracarea[m] /= area;
      
      // set 2 random tangent vectors to surf normal
      // tangent vectors are also normal to each other
      
      delta[0] = random->uniform();
      delta[1] = random->uniform();
      delta[2] = random->uniform();
      MathExtra::cross3(tris[isurf].norm,delta,tasks[ntask].tan1);
      MathExtra::norm3(tasks[ntask].tan1);
      MathExtra::cross3(tris[isurf].norm,tasks[ntask].tan1,tasks[ntask].tan2);
      MathExtra::norm3(tasks[ntask].tan2);
    }
    
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
                   "Fix emit/surf insertion count exceeds 32-bit int");
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

  // return # of tasks for this cell
  
  return ntask-ntaskorig;
}

/* ----------------------------------------------------------------------
   insert particles in grid cells with emitting surface elements
------------------------------------------------------------------------- */

void FixEmitSurf::perform_task()
{
  int i,m,n,pcell,isurf,ninsert,nactual,isp,ispecies,ntri,id;
  double indot,scosine,rn,ntarget,vr,alpha,beta;
  double beta_un,normalized_distbn_fn,theta,erot,evib;
  double vnmag,vamag,vbmag;
  double *normal,*p1,*p2,*p3,*atan,*btan,*vstream,*vscale;
  double x[3],v[3],e1[3],e2[3];
  Particle::OnePart *p;
  
  double dt = update->dt;
  int *species = particle->mixture[imix]->species;

  // if subsonic, re-compute particle inflow counts for each task
  // also computes current per-task temp_thermal and vstream
  
  if (subsonic) subsonic_inflow();
  
  // insert particles for each task = cell/surf pair
  // ntarget/ninsert is either perspecies or for all species
  // for one particle:
  //   x = random position with overlap of surf with cell
  //   v = randomized thermal velocity + vstream
  //       if normalflag, mag of vstream is applied to surf normal dir
  //       first stage: normal dimension (normal)
  //       second stage: parallel dimensions (tan1,tan2)
  
  // double while loop until randomized particle velocity meets 2 criteria
  // inner do-while loop:
  //   v = vstream-component + vthermal is into simulation box
  //   see Bird 1994, p 425
  // outer do-while loop:
  //   shift Maxwellian distribution by stream velocity component
  //   see Bird 1994, p 259, eq 12.5
  
  int nfix_add_particle = modify->n_add_particle;
  indot = magvstream;
  
  for (i = 0; i < ntask; i++) {
    pcell = tasks[i].pcell;
    isurf = tasks[i].isurf;
    if (dimension == 2) normal = lines[isurf].norm;
    else normal = tris[isurf].norm;
    atan = tasks[i].tan1;
    btan = tasks[i].tan2;
    
    temp_thermal = tasks[i].temp_thermal;
    temp_rot = tasks[i].temp_rot;
    temp_vib = tasks[i].temp_vib;
    vstream = tasks[i].vstream;
    
    if (subsonic_style == PONLY) vscale = tasks[i].vscale;
    else vscale = particle->mixture[imix]->vscale;
    
    if (!normalflag) indot = vstream[0]*normal[0] + vstream[1]*normal[1] + 
                       vstream[2]*normal[2];
    
    if (perspecies) {
      for (isp = 0; isp < nspecies; isp++) {
        ispecies = species[isp];
        ntarget = tasks[i].ntargetsp[isp]+random->uniform();
        ninsert = static_cast<int> (ntarget);
        scosine = indot / vscale[isp];
        
        nactual = 0;
        for (m = 0; m < ninsert; m++) {
          if (dimension == 2) {
            rn = random->uniform();
            p1 = &tasks[i].path[0];
            p2 = &tasks[i].path[3];
            x[0] = p1[0] + rn * (p2[0]-p1[0]);
            x[1] = p1[1] + rn * (p2[1]-p1[1]);
            x[2] = 0.0;
          } else {
            rn = random->uniform();
            ntri = tasks[i].npoint - 2;
            for (n = 0; n < ntri; n++)
              if (rn < tasks[i].fracarea[n]) break;
            p1 = &tasks[i].path[0];
            p2 = &tasks[i].path[3*(n+1)];
            p3 = &tasks[i].path[3*(n+2)];
            MathExtra::sub3(p2,p1,e1);
            MathExtra::sub3(p3,p1,e2);
            alpha = random->uniform();
            beta = random->uniform();
            if (alpha+beta > 1.0) {
              alpha = 1.0 - alpha;
              beta = 1.0 - beta;
            }
            x[0] = p1[0] + alpha*e1[0] + beta*e2[0];
            x[1] = p1[1] + alpha*e1[1] + beta*e2[1];
            x[2] = p1[2] + alpha*e1[2] + beta*e2[2];
          }
          
          if (region && !region->match(x)) continue;
          
          do {
            do beta_un = (6.0*random->uniform() - 3.0);
            while (beta_un + scosine < 0.0);
            normalized_distbn_fn = 2.0 * (beta_un + scosine) / 
              (scosine + sqrt(scosine*scosine + 2.0)) *
              exp(0.5 + (0.5*scosine)*(scosine-sqrt(scosine*scosine + 2.0)) - 
                  beta_un*beta_un);
          } while (normalized_distbn_fn < random->uniform());
          
          if (normalflag) vnmag = beta_un*vscale[isp] + magvstream;
          else vnmag = beta_un*vscale[isp] + indot;
          
          theta = MY_2PI * random->uniform();
          vr = vscale[isp] * sqrt(-log(random->uniform()));
          if (normalflag) {
            vamag = vr * sin(theta);
            vbmag = vr * cos(theta);
          } else {
            vamag = vr * sin(theta) + MathExtra::dot3(vstream,atan);
            vbmag = vr * cos(theta) + MathExtra::dot3(vstream,btan);
          }
          
          v[0] = vnmag*normal[0] + vamag*atan[0] + vbmag*btan[0];
          v[1] = vnmag*normal[1] + vamag*atan[1] + vbmag*btan[1];
          v[2] = vnmag*normal[2] + vamag*atan[2] + vbmag*btan[2];
          
          erot = particle->erot(ispecies,temp_rot,random);
          evib = particle->evib(ispecies,temp_vib,random);
          id = MAXSMALLINT*random->uniform();
          
          particle->add_particle(id,ispecies,pcell,x,v,erot,evib);
          nactual++;
          
          p = &particle->particles[particle->nlocal-1];
          p->flag = PSURF + 1 + isurf;
          p->dtremain = dt * random->uniform();
          
          if (nfix_add_particle) 
            modify->add_particle(particle->nlocal-1,temp_thermal,
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
        
        if (dimension == 2) {
          rn = random->uniform();
          p1 = &tasks[i].path[0];
          p2 = &tasks[i].path[3];
          x[0] = p1[0] + rn * (p2[0]-p1[0]);
          x[1] = p1[1] + rn * (p2[1]-p1[1]);
          x[2] = 0.0;
        } else {
          rn = random->uniform();
          ntri = tasks[i].npoint - 2;
          for (n = 0; n < ntri; n++)
            if (rn < tasks[i].fracarea[n]) break;
          p1 = &tasks[i].path[0];
          p2 = &tasks[i].path[3*(n+1)];
          p3 = &tasks[i].path[3*(n+2)];
          MathExtra::sub3(p2,p1,e1);
          MathExtra::sub3(p3,p1,e2);
          alpha = random->uniform();
          beta = random->uniform();
          if (alpha+beta > 1.0) {
            alpha = 1.0 - alpha;
            beta = 1.0 - beta;
          }
          x[0] = p1[0] + alpha*e1[0] + beta*e2[0];
          x[1] = p1[1] + alpha*e1[1] + beta*e2[1];
          x[2] = p1[2] + alpha*e1[2] + beta*e2[2];
        }
        
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
        
        if (normalflag) vnmag = beta_un*vscale[isp] + magvstream;
        else vnmag = beta_un*vscale[isp] + indot;
        
        theta = MY_2PI * random->uniform();
        vr = vscale[isp] * sqrt(-log(random->uniform()));
        if (normalflag) {
          vamag = vr * sin(theta);
          vbmag = vr * cos(theta);
        } else {
          vamag = vr * sin(theta) + MathExtra::dot3(vstream,atan);
          vbmag = vr * cos(theta) + MathExtra::dot3(vstream,btan);
        }
        
        v[0] = vnmag*normal[0] + vamag*atan[0] + vbmag*btan[0];
        v[1] = vnmag*normal[1] + vamag*atan[1] + vbmag*btan[1];
        v[2] = vnmag*normal[2] + vamag*atan[2] + vbmag*btan[2];

        erot = particle->erot(ispecies,temp_rot,random);
        evib = particle->evib(ispecies,temp_vib,random);
        id = MAXSMALLINT*random->uniform();
        
        particle->add_particle(id,ispecies,pcell,x,v,erot,evib);
        nactual++;
        
        p = &particle->particles[particle->nlocal-1];
        p->flag = PSURF + 1 + isurf;
        p->dtremain = dt * random->uniform();
        
        if (nfix_add_particle) 
          modify->add_particle(particle->nlocal-1,temp_thermal,
                               temp_rot,temp_vib,vstream);
      }
      
      nsingle += nactual;
    }
  }
}

/* ----------------------------------------------------------------------
   recalculate task properties based on subsonic BC
------------------------------------------------------------------------- */

void FixEmitSurf::subsonic_inflow()
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
    
    // indot depends on normalflag
    
    if (dimension == 2) {
      if (normalflag) indot = magvstream;
      else {
        normal = lines[tasks[i].isurf].norm;
        indot = vstream[0]*normal[0] + vstream[1]*normal[1];
      }
    } else {
      if (normalflag) indot = magvstream;
      else {
        normal = tris[tasks[i].isurf].norm;
        indot = vstream[0]*normal[0] + vstream[1]*normal[1] + 
          vstream[2]*normal[2];
      }
    }
    
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
                 "Fix emit/surf subsonic insertion count exceeds 32-bit int");
  }
}

/* ----------------------------------------------------------------------
   identify particles in grid cells associated with a task
   store count and linked list, same as for particle sorting
------------------------------------------------------------------------- */

void FixEmitSurf::subsonic_sort()
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

void FixEmitSurf::subsonic_grid()
{
  int m,ip,np,icell,ispecies;
  double mass,masstot,gamma,ke;
  double nrho_cell,massrho_cell,temp_thermal_cell,press_cell;
  double mass_cell,gamma_cell,soundspeed_cell,vsmag;
  double mv[4];
  double *v,*vstream,*vscale,*normal;
  
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

      // adjust COM vstream by difference bewteen
      //   cell pressure and subsonic target pressure
      // normal = direction of difference which depends on normalflag
      
      if (dimension == 2) {
        if (normalflag) normal = lines[tasks[i].isurf].norm;
        else normal = norm_vstream;
      } else {
        if (normalflag) normal = tris[tasks[i].isurf].norm;
        else normal = norm_vstream;
      }
      
      if (np) {
        vsmag = (psubsonic - press_cell) / (massrho_cell*soundspeed_cell);
        vstream[0] += vsmag*normal[0];
        vstream[1] += vsmag*normal[1];
        vstream[2] += vsmag*normal[2];
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
   pack one task into buf
   return # of bytes packed
   if not memflag, only return count, do not fill buf
------------------------------------------------------------------------- */

int FixEmitSurf::pack_task(int itask, char *buf, int memflag)
{
  char *ptr = buf;
  if (memflag) memcpy(ptr,&tasks[itask],sizeof(Task));
  ptr += sizeof(Task);
  ptr = ROUNDUP(ptr);
  
  // pack task vectors
  // vscale is allocated, but not communicated, since updated every step
  
  if (perspecies) {
    if (memflag) memcpy(ptr,tasks[itask].ntargetsp,nspecies*sizeof(double));
    ptr += nspecies*sizeof(double);
  }
  
  int npoint = tasks[itask].npoint;
  if (memflag) memcpy(ptr,tasks[itask].path,npoint*3*sizeof(double));
  ptr += npoint*3*sizeof(double);
  if (memflag) memcpy(ptr,tasks[itask].fracarea,(npoint-2)*sizeof(double));
  ptr += (npoint-2)*sizeof(double);
  
  return ptr-buf;
}

/* ----------------------------------------------------------------------
   unpack one task from buf
------------------------------------------------------------------------- */

int FixEmitSurf::unpack_task(char *buf, int icell)
{
  char *ptr = buf;
  
  if (ntask == ntaskmax) grow_task();
  double *ntargetsp = tasks[ntask].ntargetsp;
  double *vscale = tasks[ntask].vscale;
  double *path = tasks[ntask].path;
  double *fracarea = tasks[ntask].fracarea;
  
  memcpy(&tasks[ntask],ptr,sizeof(Task));
  ptr += sizeof(Task);
  ptr = ROUNDUP(ptr);
  
  // unpack task vectors
  // vscale is allocated, but not communicated, since updated every step
  
  if (perspecies) {
    memcpy(ntargetsp,ptr,nspecies*sizeof(double));
    ptr += nspecies*sizeof(double);
  }
  
  int npoint = tasks[ntask].npoint;
  delete [] path;
  path = new double[npoint*3];
  memcpy(path,ptr,npoint*3*sizeof(double));
  ptr += npoint*3*sizeof(double);
  
  delete [] fracarea;
  fracarea = new double[npoint-2];
  memcpy(fracarea,ptr,(npoint-2)*sizeof(double));
  ptr += (npoint-2)*sizeof(double);
  
  tasks[ntask].ntargetsp = ntargetsp;
  tasks[ntask].vscale = vscale;
  tasks[ntask].path = path;
  tasks[ntask].fracarea = fracarea;
  
  // reset task icell and pcell
  // if a split cell, set pcell via scan of icell csurfs, extract from sinfo
  
  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;
  
  tasks[ntask].icell = icell;
  if (cells[icell].nsplit == 1) tasks[ntask].pcell = icell;
  else {
    int isurf = tasks[ntask].isurf;
    int nsurf = cells[icell].nsurf;
    int *csurfs = cells[icell].csurfs;
    int i;
    for (i = 0; i < nsurf; i++)
      if (csurfs[i] == isurf) break;
    int isplit = cells[icell].isplit;
    int subcell = sinfo[isplit].csplits[i];
    tasks[ntask].pcell = sinfo[isplit].csubs[subcell];
  }
  
  ntask++;
  return ptr-buf;
}

/* ----------------------------------------------------------------------
   copy N tasks starting at index oldfirst to index first
------------------------------------------------------------------------- */

void FixEmitSurf::copy_task(int icell, int n, int first, int oldfirst)
{
  // reset icell in each copied task
  // copy task vectors
  // vscale is allocated, but not copied, since updated every step
  
  if (first == oldfirst) {
    for (int i = 0; i < n; i++) {
      tasks[first].icell = icell;
      first++;
    }
    
  } else {
    int npoint;
    for (int i = 0; i < n; i++) {
      double *ntargetsp = tasks[first].ntargetsp;
      double *vscale = tasks[first].vscale;
      double *path = tasks[first].path;
      double *fracarea = tasks[first].fracarea;
      
      memcpy(&tasks[first],&tasks[oldfirst],sizeof(Task));
      
      if (perspecies)
        memcpy(ntargetsp,tasks[oldfirst].ntargetsp,nspecies*sizeof(double));
      
      npoint = tasks[first].npoint;
      delete [] path;
      path = new double[npoint*3];
      memcpy(path,tasks[oldfirst].path,npoint*3*sizeof(double));
      
      delete [] fracarea;
      fracarea = new double[npoint-2];
      memcpy(fracarea,tasks[oldfirst].fracarea,(npoint-2)*sizeof(double));
      
      tasks[first].ntargetsp = ntargetsp;
      tasks[first].vscale = vscale;
      tasks[first].path = path;
      tasks[first].fracarea = fracarea;
      
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

void FixEmitSurf::grow_task()
{
  int oldmax = ntaskmax;
  ntaskmax += DELTATASK;
  tasks = (Task *) memory->srealloc(tasks,ntaskmax*sizeof(Task),
				    "emit/face:tasks");
  
  // set all new task bytes to 0 so valgrind won't complain
  // if bytes between fields are uninitialized
  
  memset(&tasks[oldmax],0,(ntaskmax-oldmax)*sizeof(Task));
  
  // allocate vectors in each new task or set to NULL
  // path and fracarea are allocated later to specific sizes
  
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
  
  for (int i = oldmax; i < ntaskmax; i++) {
    tasks[i].path = NULL;
    tasks[i].fracarea = NULL;
  }
}

/* ----------------------------------------------------------------------
   reset pcell for all compress task entries
   called from Grid::compress() after grid cells have been compressed
   wait to do this until now b/c split cells and their sinfo
     are setup in Grid::compress() between compress_grid() 
     and post_compress_grid()
------------------------------------------------------------------------- */

void FixEmitSurf::post_compress_grid()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;
  
  for (int i = 0; i < ntask; i++) {
    int icell = tasks[i].icell;
    if (cells[icell].nsplit == 1) tasks[i].pcell = icell;
    else {
      int isurf = tasks[i].isurf;
      int nsurf = cells[icell].nsurf;
      int *csurfs = cells[icell].csurfs;
      int j;
      for (j = 0; j < nsurf; j++)
        if (csurfs[j] == isurf) break;
      int isplit = cells[icell].isplit;
      int subcell = sinfo[isplit].csplits[j];
      tasks[i].pcell = sinfo[isplit].csubs[subcell];
    }
  }
}

/* ----------------------------------------------------------------------
   process keywords specific to this class
------------------------------------------------------------------------- */

int FixEmitSurf::option(int narg, char **arg)
{
  if (strcmp(arg[0],"n") == 0) {
    if (2 > narg) error->all(FLERR,"Illegal fix emit/surf/normal command");
    np = atoi(arg[1]);
    if (np <= 0) error->all(FLERR,"Illegal fix emit/surf/normal command");
    return 2;
  }
  
  if (strcmp(arg[0],"normal") == 0) {
    if (2 > narg) error->all(FLERR,"Illegal fix emit/surf/normal command");
    if (strcmp(arg[1],"yes") == 0) normalflag = 1;
    else if (strcmp(arg[1],"no") == 0) normalflag = 0;
    else error->all(FLERR,"Illegal fix emit/surf/normal command");
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
  
  error->all(FLERR,"Illegal fix emit/surf/normal command");
  return 0;
}
