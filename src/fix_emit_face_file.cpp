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
#include "stdlib.h"
#include "string.h"
#include "stdlib.h"
#include "fix_emit_face_file.h"
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
enum{NRHO,TEMP_THERMAL,TEMP_ROT,TEMP_VIB,VX,VY,VZ,SPECIES};

#define DELTATASK 256
#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

FixEmitFaceFile::FixEmitFaceFile(SPARTA *sparta, int narg, char **arg) :
  FixEmit(sparta, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix emit/face/file command");

  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;

  imix = particle->find_mixture(arg[2]);
  if (imix < 0) 
    error->all(FLERR,"Fix emit/face/file mixture ID does not exist");

  if (strcmp(arg[3],"xlo") == 0) iface = XLO;
  else if (strcmp(arg[3],"xhi") == 0) iface = XHI;
  else if (strcmp(arg[3],"ylo") == 0) iface = YLO;
  else if (strcmp(arg[3],"yhi") == 0) iface = YHI;
  else if (strcmp(arg[3],"zlo") == 0) iface = ZLO;
  else if (strcmp(arg[3],"zhi") == 0) iface = ZHI;
  else error->all(FLERR,"Illegal fix emit/face/file command");

  if (comm->me == 0) read_file(arg[4],arg[5]);
  bcast_mesh();
  check_mesh_values();

  // optional args

  frac_user = 1.0;

  int iarg = 6;
  options(narg-iarg,&arg[iarg]);

  // error checks
  // NOTE: could allow insertion on iface = YHI for axisymmetric
  //       if vstream[1] = 0.0, that is allowed in FixEmitFace

  if (domain->dimension == 2 && (iface == ZLO || iface == ZHI)) 
    error->all(FLERR,"Cannot use fix emit/face/file in z dimension "
               "for 2d simulation");
  if (domain->axisymmetric && (iface == YLO || iface == YHI)) 
    error->all(FLERR,"Cannot use fix emit/face/file in y dimension "
               "for axisymmetric model");

  // face properties

  normal[0] = normal[1] = normal[2] = 0.0;
  if (iface == XLO || iface == XHI) {
    ndim = 0;
    pdim = 1;
    qdim = 2;
    if (iface == XLO) normal[0] = 1.0;
    else normal[0] = -1.0;
  } else if (iface == YLO || iface == YHI) {
    ndim = 1;
    pdim = 0;
    qdim = 2;
    if (iface == YLO) normal[1] = 1.0;
    else normal[1] = -1.0;
  } else if (iface == ZLO || iface == ZHI) {
    ndim = 2;
    pdim = 0;
    qdim = 1;
    if (iface == ZLO) normal[2] = 1.0;
    else normal[2] = -1.0;
  }

  // task list

  tasks = NULL;
  ntask = ntaskmax = 0;
}

/* ---------------------------------------------------------------------- */

FixEmitFaceFile::~FixEmitFaceFile()
{
  delete [] mesh.which;
  delete [] mesh.imesh;
  delete [] mesh.jmesh;
  memory->destroy(mesh.values);

  for (int i = 0; i < ntaskmax; i++) {
    delete [] tasks[i].ntargetsp;
    delete [] tasks[i].fraction;
    delete [] tasks[i].cummulative;
    delete [] tasks[i].vscale;
  }
  memory->sfree(tasks);
}

/* ---------------------------------------------------------------------- */

void FixEmitFaceFile::init()
{
  // copies of class data before invoking parent init() and count_task()

  dimension = domain->dimension;
  fnum = update->fnum;
  dt = update->dt;

  nspecies = particle->mixture[imix]->nspecies;
  nrho_mix = particle->mixture[imix]->nrho;
  temp_thermal_mix = particle->mixture[imix]->temp_thermal;
  temp_rot_mix = particle->mixture[imix]->temp_rot;
  temp_vib_mix = particle->mixture[imix]->temp_vib;
  vstream_mix = particle->mixture[imix]->vstream;
  vscale_mix = particle->mixture[imix]->vscale;
  fraction_mix = particle->mixture[imix]->fraction;
  fraction_flag_mix = particle->mixture[imix]->fraction_flag;
  fraction_user_mix = particle->mixture[imix]->fraction_user;
  cummulative_mix = particle->mixture[imix]->cummulative;
  species2species_mix = particle->mixture[imix]->species2species;

  pts = surf->pts;
  lines = surf->lines;
  tris = surf->tris;

  // cannot inflow thru periodic boundary

  if (domain->bflag[iface] == PERIODIC)
    error->all(FLERR,"Cannot use fix emit/face/file on periodic boundary");

  // check that any species listed as a mesh value is in mixture

  int isp;

  for (int m = 0; m < mesh.nvalues; m++) {
    if (mesh.which[m] >= 0) continue;
    isp = -mesh.which[m] - 1;
    if (particle->mixture[imix]->species2species[isp] < 0)
      error->all(FLERR,"Fix inflow/file species is not in mixture");
  }

  // reallocate ntargetsp for each task 
  // b/c nspecies count of mixture may have changed

  for (int i = 0; i < ntask; i++) {
    delete [] tasks[i].ntargetsp;
    delete [] tasks[i].fraction;
    delete [] tasks[i].cummulative;
    delete [] tasks[i].vscale;
    tasks[i].ntargetsp = new double[nspecies];
    tasks[i].fraction = new double[nspecies];
    tasks[i].cummulative = new double[nspecies];
    tasks[i].vscale = new double[nspecies];
  }

  // per-species vectors for mesh setting of species fractions
  // initialize to mixture settings

  fflag = new int[nspecies];
  fuser = new double[nspecies];
  for (isp = 0; isp < nspecies; isp++) {
    fflag[isp] = fraction_flag_mix[isp];
    fuser[isp] = fraction_user_mix[isp];
  }

  // invoke FixEmit::init() to populate task list
  // it calls create_task() for each grid cell

  ntask = 0;
  FixEmit::init();

  // clean up

  delete [] fflag;
  delete [] fuser;
}

/* ---------------------------------------------------------------------- */

int FixEmitFaceFile::create_task(int icell)
{
  int j,n,extflag;
  int *cflags;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  // corners[i][j] = J corner points of face I of a grid cell
  // works for 2d quads and 3d hexes
  
  int corners[6][4] = {{0,2,4,6}, {1,3,5,7}, {0,1,4,5}, {2,3,6,7}, 
		       {0,1,2,3}, {4,5,6,7}};
  int nface_pts = 4;
  if (domain->dimension == 2) nface_pts = 2;

  // flag = 1 if insertion happens on iface of cell
  // only if face adjoins global boundary
  // if cell is OVERLAP:
  //   allow if any face corner point is OUTSIDE and none is INSIDE
  //   disallow if any pt of any line/tri in cell touches face
  
  int nmask = cells[icell].nmask;

  int flag = 0;
  if (grid->neigh_decode(nmask,iface) == NBOUND) {
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
  
  if (!flag) return 0;

  // set cell parameters of task
  // pcell = sub cell for particles if a split cell
  
  if (ntask == ntaskmax) grow_task();

  tasks[ntask].icell = icell;
  if (cells[icell].nsplit > 1) tasks[ntask].pcell = split(icell);
  else tasks[ntask].pcell = icell;

  // interpolate remaining task values from mesh to cell face
  // interpolate return 1 if task is valid, else 0

  flag = interpolate(icell);
  if (!flag) return 0;

  // increment task counter

  ntask++;
  return 1;
}

/* ----------------------------------------------------------------------
   insert particles in grid cells with faces touching inflow boundary
   NOTE:
     currently not allowing particle insertion on backflow boundaries
     see comment in FixEmitFace
------------------------------------------------------------------------- */

void FixEmitFaceFile::perform_task()
{
  int pcell,ninsert,nactual,isp,ispecies,ndim,id;
  double temp_thermal,temp_rot,temp_vib;
  double indot,scosine,rn,ntarget,vr;
  double beta_un,normalized_distbn_fn,theta,erot,evib;
  double x[3],v[3];
  double *lo,*hi,*vstream,*cummulative,*vscale;
  Particle::OnePart *p;

  double dt = update->dt;
  int *species = particle->mixture[imix]->species;

  // insert particles for each task = cell
  // ntarget/ninsert is either perspecies or for all species
  // for one particle:
  //   x = random position on subset of face that overlaps with file grid
  //   v = randomized thermal velocity + vstream
  //       first stage: normal dimension (ndim)
  //       second stage: parallel dimensions (pdim1,pdim2)

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
    lo = tasks[i].lo;
    hi = tasks[i].hi;

    vstream = tasks[i].vstream;
    indot = vstream[0]*normal[0] + vstream[1]*normal[1] + vstream[2]*normal[2];
    temp_thermal = tasks[i].temp_thermal;
    temp_rot = tasks[i].temp_rot;
    temp_vib = tasks[i].temp_vib;
    vscale = tasks[i].vscale;

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

          theta = MY_PI * random->gaussian();
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

          if (nfix_add_particle) 
            modify->add_particle(particle->nlocal-1,temp_thermal,
                                 temp_rot,temp_vib,vstream);
	}

	nsingle += nactual;
      }

    } else {
      cummulative = tasks[i].cummulative;
      ntarget = tasks[i].ntarget+random->uniform();
      ninsert = static_cast<int> (ntarget);

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
	  do beta_un = (6.0*random->gaussian() - 3.0);
	  while (beta_un + scosine < 0.0);
	  normalized_distbn_fn = 2.0 * (beta_un + scosine) / 
	    (scosine + sqrt(scosine*scosine + 2.0)) *
	    exp(0.5 + (0.5*scosine)*(scosine-sqrt(scosine*scosine + 2.0)) - 
		beta_un*beta_un);
	} while (normalized_distbn_fn < random->uniform());
	
        v[ndim] = beta_un*vscale[isp]*normal[ndim] + vstream[ndim];

        theta = MY_PI * random->gaussian();
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

        if (nfix_add_particle) 
          modify->add_particle(particle->nlocal-1,temp_thermal,
                               temp_rot,temp_vib,vstream);
      }

      nsingle += nactual;
    }
  }
}

/* ----------------------------------------------------------------------
   scan file for section-ID, read regular grid of values into Mesh data struct
   only called by proc 0
------------------------------------------------------------------------- */

void FixEmitFaceFile::read_file(char *file, char *section)
{
  int i,m,n,ii,jj,offset;
  char line[MAXLINE];
  char *word;

  int dimension = domain->dimension;

  // open file

  FILE *fp = fopen(file,"r");
  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open inflow file %s",file);
    error->one(FLERR,str);
  }

  // loop until section found with matching keyword

  while (1) {
    if (fgets(line,MAXLINE,fp) == NULL)
      error->one(FLERR,"Did not find keyword in inflow file");
    if (strspn(line," \t\n\r") == strlen(line)) continue;  // blank line
    if (line[0] == '#') continue;                          // comment
    word = strtok(line," \t\n\r");
    if (strcmp(word,section) == 0) break;           // matching keyword

    fgets(line,MAXLINE,fp);                         // no match, read NIJ or NI
    word = strtok(line," \t\n\r");                  // skip 2d or 3d section
    int nskip;
    if (strcmp(word,"NIJ") == 0) {
      word = strtok(NULL," \t\n\r");
      nskip = atoi(word);
      word = strtok(NULL," \t\n\r");
      nskip *= atoi(word);
      fgets(line,MAXLINE,fp);                         // values line
      fgets(line,MAXLINE,fp);                         // imesh line
      fgets(line,MAXLINE,fp);                         // jmesh line
    } else if (strcmp(word,"NI") == 0) {
      word = strtok(NULL," \t\n\r");
      nskip = atoi(word);
      fgets(line,MAXLINE,fp);                         // values line
      fgets(line,MAXLINE,fp);                         // imesh line
    } else error->one(FLERR,"Misformatted section in inflow file");

    fgets(line,MAXLINE,fp);                               // blank line
    for (i = 0; i < nskip; i++) fgets(line,MAXLINE,fp);   // value lines
  }

  // read and store the matching section

  fgets(line,MAXLINE,fp);                             // read NIJ or NI
  word = strtok(line," \t\n\r");
  if (strcmp(word,"NIJ") == 0) {
    if (dimension != 3) 
      error->all(FLERR,"Misformatted section in inflow file");
    word = strtok(NULL," \t\n\r");
    mesh.ni = atoi(word);
    word = strtok(NULL," \t\n\r");
    mesh.nj = atoi(word);
  } else if (strcmp(word,"NI") == 0) {
    if (dimension != 2) 
      error->all(FLERR,"Misformatted section in inflow file");
    word = strtok(NULL," \t\n\r");
    mesh.ni = atoi(word);
    mesh.nj = 1;
  }

  if (mesh.ni < 2 || (dimension == 3 && mesh.nj < 2))
    error->one(FLERR,"Inflow file grid is too small");

  // read NV line

  fgets(line,MAXLINE,fp);
  word = strtok(line," \t\n\r");
  if (strcmp(word,"NV") != 0)
    error->one(FLERR,"Misformatted section in inflow file");
  word = strtok(NULL," \t\n\r");
  mesh.nvalues = atoi(word);
  if (mesh.nvalues <= 0)
    error->one(FLERR,"Misformatted section in inflow file");

  // read VALUES line and convert names to which vector

  mesh.which = new int[mesh.nvalues];
  fgets(line,MAXLINE,fp);
  word = strtok(line," \t\n\r");
  for (i = 0; i < mesh.nvalues; i++) {
    word = strtok(NULL," \t\n\r");
    if (strcmp(word,"nrho") == 0) mesh.which[i] = NRHO;
    else if (strcmp(word,"temp") == 0) mesh.which[i] = TEMP_THERMAL;
    else if (strcmp(word,"trot") == 0) mesh.which[i] = TEMP_ROT;
    else if (strcmp(word,"tvib") == 0) mesh.which[i] = TEMP_VIB;
    else if (strcmp(word,"vx") == 0) mesh.which[i] = VX;
    else if (strcmp(word,"vy") == 0) mesh.which[i] = VY;
    else if (strcmp(word,"vz") == 0) mesh.which[i] = VZ;
    else {
      int index = particle->find_species(word);
      if (index < 0) error->one(FLERR,"Unknown species in inflow file");
      mesh.which[i] = -(index+1);
    }
  }

  // read IMESH,JMESH coords

  mesh.imesh = new double[mesh.ni];
  mesh.jmesh = new double[mesh.nj];

  fgets(line,MAXLINE,fp);
  word = strtok(line," \t\n\r");
  if (strcmp(word,"IMESH") != 0) 
    error->one(FLERR,"Misformatted section in inflow file");
  for (i = 0; i < mesh.ni; i++) {
    word = strtok(NULL," \t\n\r");
    mesh.imesh[i] = atof(word);
    if (i && mesh.imesh[i] <= mesh.imesh[i-1])
      error->one(FLERR,"Misformatted section in inflow file");
  }
  mesh.lo[0] = mesh.imesh[0];
  mesh.hi[0] = mesh.imesh[mesh.ni-1];

  if (dimension == 3) {
    fgets(line,MAXLINE,fp);
    word = strtok(line," \t\n\r");
    if (strcmp(word,"JMESH") != 0) 
      error->one(FLERR,"Misformatted section in inflow file");
    for (i = 0; i < mesh.nj; i++) {
      word = strtok(NULL," \t\n\r");
      mesh.jmesh[i] = atof(word);
      if (i && mesh.jmesh[i] <= mesh.jmesh[i-1])
        error->one(FLERR,"Misformatted section in inflow file");
    }
  } else mesh.jmesh[0] = 0.0;
  mesh.lo[1] = mesh.jmesh[0];
  mesh.hi[1] = mesh.jmesh[mesh.nj-1];

  // N = Ni by Nj values lines, store values in mesh
  // II,JJ stored with II varying fastest in 2d

  fgets(line,MAXLINE,fp);    // blank line

  n = mesh.ni * mesh.nj;
  memory->create(mesh.values,n,mesh.nvalues,"inflow/file:values");

  for (i = 0; i < n; i++) {
    fgets(line,MAXLINE,fp);
    word = strtok(line," \t\n\r");
    ii = atoi(word);
    if (dimension == 3) {
      word = strtok(NULL," \t\n\r");
      jj = atoi(word);
    } else jj = 1;
    if (ii < 1 || ii > mesh.ni || jj < 1 || jj > mesh.nj)
      error->one(FLERR,"Misformatted section in inflow file");
    offset = (jj-1)*mesh.ni + (ii-1);   
    for (m = 0; m < mesh.nvalues; m++) {
      word = strtok(NULL," \t\n\r");
      mesh.values[offset][m] = atof(word);
    }
  }

  fclose(fp);
}

/* ----------------------------------------------------------------------
   bcast mesh data struct from proc 0 to all other procs
------------------------------------------------------------------------- */

void FixEmitFaceFile::bcast_mesh()
{
  int n;

  int me = comm->me;
  MPI_Bcast(&mesh.ni,1,MPI_INT,0,world);
  MPI_Bcast(&mesh.nj,1,MPI_INT,0,world);
  MPI_Bcast(&mesh.nvalues,1,MPI_INT,0,world);
  MPI_Bcast(mesh.lo,2,MPI_DOUBLE,0,world);
  MPI_Bcast(mesh.hi,2,MPI_DOUBLE,0,world);

  if (me) {
    mesh.which = new int[mesh.nvalues];
    mesh.imesh = new double[mesh.ni];
    mesh.jmesh = new double[mesh.nj];
    n = mesh.ni * mesh.nj;
    memory->create(mesh.values,n,mesh.nvalues,"inflow/file:values");
  }

  MPI_Bcast(mesh.which,mesh.nvalues,MPI_INT,0,world);
  MPI_Bcast(mesh.imesh,mesh.ni,MPI_DOUBLE,0,world);
  MPI_Bcast(mesh.jmesh,mesh.nj,MPI_DOUBLE,0,world);
  MPI_Bcast(&mesh.values[0][0],mesh.ni*mesh.nj*mesh.nvalues,
            MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
   check that values on mesh points are valid
   all but VX,VY,VZ must be >= 0.0
   species fractions must be <= 1.0
------------------------------------------------------------------------- */

void FixEmitFaceFile::check_mesh_values()
{
  int n = mesh.ni;
  if (domain->dimension == 3) n *= mesh.nj;

  int flag = 0;
  for (int m = 0; m < mesh.nvalues; m++) {
    if (mesh.which[m] == VX || mesh.which[m] == VY || mesh.which[m] == VZ)
      continue;
    for (int i = 0; i < n; i++) {
      if (mesh.values[i][m] < 0.0) flag = 1;
    }
    if (mesh.which[m] < 0) {
      for (int i = 0; i < n; i++)
        if (mesh.values[i][m] > 1.0) flag = 1;
    }
  }

  if (flag) error->all(FLERR,"Fix inflow/file mesh has invalid value");
}

/* ----------------------------------------------------------------------
   interpolate values needed by iface of icell
   use centroid of face to find 4 surrounding mesh points (2 in 2d)
   use bilinear interpolation (linear in 2d) to compute each new value
   return 1 if interpolation successful, 0 if not
------------------------------------------------------------------------- */

int FixEmitFaceFile::interpolate(int icell)
{
  int j,m,isp,plo,phi,qlo,qhi,anyfrac,err;
  double indot,newtemp,area;
  double *lo,*hi;
  double xc[2];

  // task lo/hi = extent of overlap
  // if lo >= hi, then no overlap

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  lo = cells[icell].lo;
  hi = cells[icell].hi;
    
  if (iface % 2 == 0) 
    tasks[ntask].lo[ndim] = tasks[ntask].hi[ndim] = lo[ndim];
  else 
    tasks[ntask].lo[ndim] = tasks[ntask].hi[ndim] = hi[ndim];
  
  tasks[ntask].lo[pdim] = MAX(lo[pdim],mesh.lo[0]);
  tasks[ntask].hi[pdim] = MIN(hi[pdim],mesh.hi[0]);
  if (tasks[ntask].lo[pdim] >= tasks[ntask].hi[pdim]) return 0;
  
  if (dimension == 3) {
    tasks[ntask].lo[qdim] = MAX(lo[qdim],mesh.lo[1]);
    tasks[ntask].hi[qdim] = MIN(hi[qdim],mesh.hi[1]);
    if (tasks[ntask].lo[qdim] >= tasks[ntask].hi[qdim]) return 0;
  }
    
  // default task params from mixture
  
  tasks[ntask].nrho = nrho_mix;
  tasks[ntask].vstream[0] = vstream_mix[0];
  tasks[ntask].vstream[1] = vstream_mix[1];
  tasks[ntask].vstream[2] = vstream_mix[2];
  tasks[ntask].temp_thermal = temp_thermal_mix;
  tasks[ntask].temp_rot = temp_rot_mix;
  tasks[ntask].temp_vib = temp_vib_mix;
  for (j = 0; j < nspecies; j++) {
    tasks[ntask].fraction[j] = fraction_mix[j];
    tasks[ntask].cummulative[j] = cummulative_mix[j];
    tasks[ntask].vscale[j] = vscale_mix[j];
  }

  // xc = centroid of cell face overlap with mesh
  
  xc[0] = 0.5 * (tasks[ntask].lo[pdim]+tasks[ntask].hi[pdim]);
  if (dimension == 3)
    xc[1] = 0.5 * (tasks[ntask].lo[qdim]+tasks[ntask].hi[qdim]);
  
  // find 2/4 points in mesh that surround xc in 2d/3d
  // allow for mesh.lo <= xc <= mesh.hi

  for (m = 1; m < mesh.ni; m++)
    if (xc[0] <= mesh.imesh[m]) break;
  plo = m-1;
  phi = m;

  if (dimension == 3) {
    for (m = 1; m < mesh.nj; m++)
      if (xc[1] <= mesh.jmesh[m]) break;
    qlo = m-1;
    qhi = m;
  }

  // override default params with file values
  // 2d or 3d = linear or bilinear interpolation
  // if temp_thermal is set, vscale must also be reset

  anyfrac = 0;

  if (dimension == 2) {
    for (m = 0; m < mesh.nvalues; m++) {
      if (mesh.which[m] == NRHO)
        tasks[ntask].nrho = linear_interpolation(xc[0],m,plo,phi);
      else if (mesh.which[m] == TEMP_THERMAL) {
        newtemp = tasks[ntask].temp_thermal = 
          linear_interpolation(xc[0],m,plo,phi);
        for (isp = 0; isp < nspecies; isp++)
          tasks[ntask].vscale[isp] = 
            vscale_mix[isp] * sqrt(newtemp/temp_thermal_mix);
      }
      else if (mesh.which[m] == TEMP_ROT)
        tasks[ntask].temp_rot = linear_interpolation(xc[0],m,plo,phi);
      else if (mesh.which[m] == TEMP_VIB)
        tasks[ntask].temp_vib = linear_interpolation(xc[0],m,plo,phi);
      else if (mesh.which[m] == VX)
        tasks[ntask].vstream[0] = linear_interpolation(xc[0],m,plo,phi);
      else if (mesh.which[m] == VY)
        tasks[ntask].vstream[1] = linear_interpolation(xc[0],m,plo,phi);
      else if (mesh.which[m] == VZ)
        tasks[ntask].vstream[2] = linear_interpolation(xc[0],m,plo,phi);
      else {
        anyfrac = 1;
        isp = -mesh.which[m] - 1;
        isp = species2species_mix[isp];
        fuser[isp] = linear_interpolation(xc[0],m,plo,phi);
        fflag[isp] = 1;
      }
    }
  } else {
    for (m = 0; m < mesh.nvalues; m++) {
      if (mesh.which[m] == NRHO)
        tasks[ntask].nrho = 
          bilinear_interpolation(xc[0],xc[1],m,plo,phi,qlo,qhi);
      else if (mesh.which[m] == TEMP_THERMAL) {
        newtemp = tasks[ntask].temp_thermal = 
          bilinear_interpolation(xc[0],xc[1],m,plo,phi,qlo,qhi);
        for (isp = 0; isp < nspecies; isp++)
          tasks[ntask].vscale[isp] = 
            vscale_mix[isp] * sqrt(newtemp/temp_thermal_mix);
      }
      else if (mesh.which[m] == TEMP_ROT)
        tasks[ntask].temp_rot = 
          bilinear_interpolation(xc[0],xc[1],m,plo,phi,qlo,qhi);
      else if (mesh.which[m] == TEMP_VIB)
        tasks[ntask].temp_vib = 
          bilinear_interpolation(xc[0],xc[1],m,plo,phi,qlo,qhi);
      else if (mesh.which[m] == VX)
        tasks[ntask].vstream[0] =
          bilinear_interpolation(xc[0],xc[1],m,plo,phi,qlo,qhi);
      else if (mesh.which[m] == VY)
        tasks[ntask].vstream[1] =
          bilinear_interpolation(xc[0],xc[1],m,plo,phi,qlo,qhi);
      else if (mesh.which[m] == VZ)
        tasks[ntask].vstream[2] =
          bilinear_interpolation(xc[0],xc[1],m,plo,phi,qlo,qhi);
      else {
        anyfrac = 1;
        isp = -mesh.which[m] - 1;
        isp = species2species_mix[isp];
        fuser[isp] = bilinear_interpolation(xc[0],xc[1],m,plo,phi,qlo,qhi);
        fflag[isp] = 1;
      }
    }
  }

  // indot = dot product of vstream with inward face normal

  tasks[ntask].ntarget = 0.0;
  indot = tasks[ntask].vstream[0]*normal[0] +
    tasks[ntask].vstream[1]*normal[1] +
    tasks[ntask].vstream[2]*normal[2];

  // if any species fraction was set,
  // set entire fraction and cummulative vector via Mixture::init_fraction()

  if (anyfrac) {
    err = particle->mixture[imix]->
      init_fraction(fflag,fuser,
                    tasks[ntask].fraction,tasks[ntask].cummulative);
    if (err)
      error->one(FLERR,"Fix inflow/file mixture fractions exceed 1.0");
  }
  
  // area = area for insertion
  // depends on dimension and axisymmetry

  if (dimension == 3) 
    area = (tasks[ntask].hi[pdim]-tasks[ntask].lo[pdim]) * 
      (tasks[ntask].hi[qdim]-tasks[ntask].lo[qdim]);
  else if (domain->axisymmetric)
    area = (tasks[ntask].hi[pdim]*tasks[ntask].hi[pdim] -
            tasks[ntask].lo[pdim]*tasks[ntask].lo[pdim])*MY_PI;
  else area = tasks[ntask].hi[pdim]-tasks[ntask].lo[pdim];

  // set ntarget and ntargetsp via mol_inflow()
  // also scale ntarget by frac_user (0 to 1)
  // skip task if final ntarget = 0.0, due to large outbound vstream

  for (isp = 0; isp < nspecies; isp++) {
    tasks[ntask].ntargetsp[isp] = frac_user *
      mol_inflow(indot,tasks[ntask].vscale[isp],
                 tasks[ntask].fraction[isp]);
    tasks[ntask].ntargetsp[isp] *= tasks[ntask].nrho*area*dt / fnum;
    tasks[ntask].ntargetsp[isp] /= cinfo[icell].weight;
    tasks[ntask].ntarget += tasks[ntask].ntargetsp[isp];
  }
  if (tasks[ntask].ntarget == 0.0) return 0;

  return 1;
}

/* ----------------------------------------------------------------------
   linear interpolation at x between lo and hi bounds, for column M
------------------------------------------------------------------------- */

double FixEmitFaceFile::linear_interpolation(double x, int m, int plo, int phi)
{
  double *imesh = mesh.imesh;
  double **values = mesh.values;
  double value = (values[plo][m]*(imesh[phi]-x) + 
                  values[phi][m]*(x-imesh[plo])) / (imesh[phi]-imesh[plo]);
  return value;
}

/* ----------------------------------------------------------------------
   bilinear interpolation at x,y between (plo,qlo) to (phi,qhi) corners
   for column M
------------------------------------------------------------------------- */

double FixEmitFaceFile::bilinear_interpolation(double x, double y, int m, 
                                             int plo, int phi, int qlo, int qhi)
{
  int ni = mesh.ni;
  double *imesh = mesh.imesh;
  double *jmesh = mesh.jmesh;
  double **values = mesh.values;

  double area = (imesh[phi]-imesh[plo]) * (jmesh[qhi]-jmesh[qlo]);
  double quad1 = values[qlo*ni+plo][m] * (imesh[phi]-x) * (jmesh[qhi]-y);
  double quad2 = values[qlo*ni+phi][m] * (x-imesh[plo]) * (jmesh[qhi]-y);
  double quad3 = values[qhi*ni+phi][m] * (x-imesh[plo]) * (y-jmesh[qlo]);
  double quad4 = values[qhi*ni+plo][m] * (imesh[phi]-x) * (y-jmesh[qlo]);
  double value = (quad1 + quad2 + quad3 + quad4) / area;

  return value;
}

/* ----------------------------------------------------------------------
   inserting into split cell icell on face iface
   determine which sub cell the face is part of
   face cannot be touched by surfs, so entire face is part of one sub cell
   compute which via update->split() and return it
------------------------------------------------------------------------- */

int FixEmitFaceFile::split(int icell)
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

int FixEmitFaceFile::pack_task(int itask, char *buf, int memflag)
{
  char *ptr = buf;
  if (memflag) memcpy(ptr,&tasks[itask],sizeof(Task));
  ptr += sizeof(Task);
  ptr = ROUNDUP(ptr);

  // pack task vectors

  if (memflag) memcpy(ptr,tasks[itask].ntargetsp,nspecies*sizeof(double));
  ptr += nspecies*sizeof(double);
  if (memflag) memcpy(ptr,tasks[itask].fraction,nspecies*sizeof(double));
  ptr += nspecies*sizeof(double);
  if (memflag) memcpy(ptr,tasks[itask].cummulative,nspecies*sizeof(double));
  ptr += nspecies*sizeof(double);
  if (memflag) memcpy(ptr,tasks[itask].vscale,nspecies*sizeof(double));
  ptr += nspecies*sizeof(double);

  return ptr-buf;
}

/* ----------------------------------------------------------------------
   unpack one task from buf
------------------------------------------------------------------------- */

int FixEmitFaceFile::unpack_task(char *buf, int icell)
{
  char *ptr = buf;

  if (ntask == ntaskmax) grow_task();
  double *ntargetsp = tasks[ntask].ntargetsp;
  double *fraction = tasks[ntask].fraction;
  double *cummulative = tasks[ntask].cummulative;
  double *vscale = tasks[ntask].vscale;

  memcpy(&tasks[ntask],ptr,sizeof(Task));
  ptr += sizeof(Task);
  ptr = ROUNDUP(ptr);

  // unpack task vectors

  memcpy(ntargetsp,ptr,nspecies*sizeof(double));
  ptr += nspecies*sizeof(double);
  memcpy(fraction,ptr,nspecies*sizeof(double));
  ptr += nspecies*sizeof(double);
  memcpy(cummulative,ptr,nspecies*sizeof(double));
  ptr += nspecies*sizeof(double);
  memcpy(vscale,ptr,nspecies*sizeof(double));
  ptr += nspecies*sizeof(double);

  tasks[ntask].ntargetsp = ntargetsp;
  tasks[ntask].fraction = fraction;
  tasks[ntask].cummulative = cummulative;
  tasks[ntask].vscale = vscale;

  // reset task icell and pcell
  // if a split cell, set pcell via split() which calls update->split()
  //   which will use current sub cells of icell
 
  tasks[ntask].icell = icell;
  if (grid->cells[icell].nsplit == 1) tasks[ntask].pcell = icell;
  else tasks[ntask].pcell = split(icell);

  ntask++;
  return ptr-buf;
}

/* ----------------------------------------------------------------------
   copy N tasks starting at index oldfirst to index first
------------------------------------------------------------------------- */

void FixEmitFaceFile::copy_task(int icell, int n, int first, int oldfirst)
{
  // reset icell in each copied task
  // copy task vectors

  if (first == oldfirst) {
    for (int i = 0; i < n; i++) {
      tasks[first].icell = icell;
      first++;
    }

  } else {
    for (int i = 0; i < n; i++) {
      double *ntargetsp = tasks[first].ntargetsp;
      double *fraction = tasks[first].fraction;
      double *cummulative = tasks[first].cummulative;
      double *vscale = tasks[first].vscale;

      memcpy(&tasks[first],&tasks[oldfirst],sizeof(Task));
      memcpy(ntargetsp,tasks[oldfirst].ntargetsp,nspecies*sizeof(double));
      memcpy(fraction,tasks[oldfirst].fraction,nspecies*sizeof(double));
      memcpy(cummulative,tasks[oldfirst].cummulative,nspecies*sizeof(double));
      memcpy(vscale,tasks[oldfirst].vscale,nspecies*sizeof(double));

      tasks[first].ntargetsp = ntargetsp;
      tasks[first].fraction = fraction;
      tasks[first].cummulative = cummulative;
      tasks[first].vscale = vscale;

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

void FixEmitFaceFile::grow_task()
{
  int oldmax = ntaskmax;
  ntaskmax += DELTATASK;
  tasks = (Task *) memory->srealloc(tasks,ntaskmax*sizeof(Task),
				    "emit/face/file:tasks");

  // set all new task bytes to 0 so valgrind won't complain
  // if bytes between fields are uninitialized

  memset(&tasks[oldmax],0,(ntaskmax-oldmax)*sizeof(Task));

  // allocate vectors in each new task

  for (int i = oldmax; i < ntaskmax; i++) {
    tasks[i].ntargetsp = new double[nspecies];
    tasks[i].fraction = new double[nspecies];
    tasks[i].cummulative = new double[nspecies];
    tasks[i].vscale = new double[nspecies];
  }
}

/* ----------------------------------------------------------------------
   reset pcell for all compress task entries
   called from Grid::compress() after grid cells have been compressed
   wait to do this until now b/c split cells accessed by split()
     are setup in Grid::compress() between compress_grid() 
     and post_compress_grid()
------------------------------------------------------------------------- */

void FixEmitFaceFile::post_compress_grid()
{
  Grid::ChildCell *cells = grid->cells;

  for (int i = 0; i < ntask; i++) {
    int icell = tasks[i].icell;
    if (cells[icell].nsplit == 1) tasks[i].pcell = icell;
    else tasks[i].pcell = split(icell);
  }
}

/* ----------------------------------------------------------------------
   process keywords specific to this class
------------------------------------------------------------------------- */

int FixEmitFaceFile::option(int narg, char **arg)
{
  if (strcmp(arg[0],"frac") == 0) {
    if (2 > narg) error->all(FLERR,"Illegal fix emit/face/file command");
    frac_user = atof(arg[1]);
    if (frac_user < 0.0 || frac_user > 1.0) 
      error->all(FLERR,"Illegal fix emit/face/file command");
    return 2;
  }

  error->all(FLERR,"Illegal fix emit/face command");
  return 0;
}

/* ----------------------------------------------------------------------
   DEBUG method: print status of cellface I
------------------------------------------------------------------------- */

void FixEmitFaceFile::print_task(int i)
{
  printf("TASK: proc %d index %d pcell %d icell %d\n",
         comm->me,i,tasks[i].pcell,tasks[i].icell);
  printf("TASK: proc %d index %d lo %g %g %g hi %g %g %g\n",
         comm->me,i,
         tasks[i].lo[0],tasks[i].lo[1],tasks[i].lo[2],
         tasks[i].hi[0],tasks[i].hi[1],tasks[i].hi[2]);
  printf("TASK: proc %d index %d ntarget %g %g %g nrho %g vstream %g %g %g\n",
         comm->me,i,tasks[i].ntarget,
         tasks[i].ntargetsp[0],tasks[i].ntargetsp[1],
         tasks[i].nrho,
         tasks[i].vstream[0],tasks[i].vstream[1],tasks[i].vstream[2]);
}
