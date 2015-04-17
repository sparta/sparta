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
#include "fix_inflow_file.h"
#include "update.h"
#include "particle.h"
#include "mixture.h"
#include "domain.h"
#include "modify.h"
#include "grid.h"
#include "surf.h"
#include "comm.h"
#include "geometry.h"
#include "random_mars.h"
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
enum{NO,YES};

enum{NRHO,TEMP_THERMAL,VX,VY,VZ,SPECIES};

#define DELTAFACE 10              // no smaller than 6
#define DELTAGRID 1000            // must be bigger than split cells per cell

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

FixInflowFile::FixInflowFile(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix inflow/file command");

  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  gridmigrate = 1;

  imix = particle->find_mixture(arg[2]);
  if (imix < 0) error->all(FLERR,"Fix inflow/file mixture ID does not exist");

  if (strcmp(arg[3],"xlo") == 0) face = XLO;
  else if (strcmp(arg[3],"xhi") == 0) face = XHI;
  else if (strcmp(arg[3],"ylo") == 0) face = YLO;
  else if (strcmp(arg[3],"yhi") == 0) face = YHI;
  else if (strcmp(arg[3],"zlo") == 0) face = ZLO;
  else if (strcmp(arg[3],"zhi") == 0) face = ZHI;
  else error->all(FLERR,"Illegal fix inflow/file command");

  if (comm->me == 0) read_file(arg[4],arg[5]);
  bcast_mesh();
  check_mesh_values();

  // optional args

  frac_user = 1.0;
  nevery = 1;
  perspecies = YES;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"frac") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix inflow/file command");
      frac_user = atof(arg[iarg+1]);
      if (frac_user < 0.0 || frac_user > 1.0) 
        error->all(FLERR,"Illegal fix inflow/file command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"nevery") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix inflow/file command");
      nevery = atoi(arg[iarg+1]);
      if (nevery <= 0) error->all(FLERR,"Illegal fix inflow/file command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"perspecies") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix inflow/file command");
      if (strcmp(arg[iarg+1],"yes") == 0) perspecies = YES;
      else if (strcmp(arg[iarg+1],"no") == 0) perspecies = NO;
      else error->all(FLERR,"Illegal fix inflow/file command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix inflow/file command");
  }

  // error check

  if (domain->dimension == 2 && (face == ZLO || face == ZHI)) 
    error->all(FLERR,"Cannot use fix inflow/file in z dimension "
               "for 2d simulation");
  if (domain->axisymmetric && (face == YLO || face == YHI)) 
    error->all(FLERR,"Cannot use fix inflow/file in y dimension "
               "for axisymmetric");

  // RNG

  int me = comm->me;
  random = new RanPark(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,me,100);

  // face properties

  normal[0] = normal[1] = normal[2] = 0.0;
  if (face == XLO || face == XHI) {
    ndim = 0;
    pdim = 1;
    qdim = 2;
    if (face == XLO) normal[0] = 1.0;
    else normal[0] = -1.0;
  } else if (face == YLO || face == YHI) {
    ndim = 1;
    pdim = 0;
    qdim = 2;
    if (face == YLO) normal[1] = 1.0;
    else normal[1] = -1.0;
  } else if (face == ZLO || face == ZHI) {
    ndim = 2;
    pdim = 0;
    qdim = 1;
    if (face == ZLO) normal[2] = 1.0;
    else normal[2] = -1.0;
  }

  // local storage

  cellface = NULL;
  ncf = ncfmax = 0;
  c2f = NULL;
  nglocal = nglocalmax = 0;

  // counters

  nsingle = ntotal = 0;
}

/* ---------------------------------------------------------------------- */

FixInflowFile::~FixInflowFile()
{
  delete random;

  delete [] mesh.which;
  delete [] mesh.imesh;
  delete [] mesh.jmesh;
  memory->destroy(mesh.values);

  for (int i = 0; i < ncfmax; i++) {
    delete [] cellface[i].ntargetsp;
    delete [] cellface[i].fraction;
    delete [] cellface[i].cummulative;
    delete [] cellface[i].vscale;
  }
  memory->sfree(cellface);
  memory->destroy(c2f);
}

/* ---------------------------------------------------------------------- */

int FixInflowFile::setmask()
{
  int mask = 0;
  mask |= START_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixInflowFile::init()
{
  int i,j,m,n,isp,icell;
  double xface[3];
  double *vscale = particle->mixture[imix]->vscale;

  particle->exist = 1;

  // corners[i][j] = J corner points of face I of a grid cell
  // works for 2d quads and 3d hexes
  
  int corners[6][4] = {{0,2,4,6}, {1,3,5,7}, {0,1,4,5}, {2,3,6,7}, 
		       {0,1,2,3}, {4,5,6,7}};

  int nface_pts = 4;
  if (domain->dimension == 2) nface_pts = 2;

  // cannot inflow thru periodic boundary

  if (domain->bflag[face] == PERIODIC)
    error->all(FLERR,"Cannot use fix inflow/file on periodic boundary");

  // c2f[I] = 1 if my local cell I allows insertions on its inflow face
  // only allow if face adjoins global boundary for inflow
  // if cell is OUTSIDE, allow insertions
  // if cell is INSIDE, disallow insertions
  // if cell is OVERLAP:
  //   allow if any face corner point is OUTSIDE and none is INSIDE
  //   disallow if any pt of any cell line/tri touches face

  int dimension = domain->dimension;
  Surf::Point *pts = surf->pts;
  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  nglocal = grid->nlocal;

  // realloc c2f if necessary = pointers from cells to cellface

  if (nglocal > nglocalmax) {
    memory->destroy(c2f);
    nglocalmax = nglocal;
    memory->create(c2f,nglocalmax,"inflow/file:c2f");
  }

  int *flags;
  int nmask,dim,extflag;
  double value;

  // set c2f = 1 if face is eligible for insertion, else 0

  for (icell = 0; icell < nglocal; icell++) {
    nmask = cells[icell].nmask;
    if (grid->neigh_decode(nmask,face) == NBOUND && 
        cells[icell].nsplit >= 1) {
      if (cinfo[icell].type == OUTSIDE) c2f[icell] = 1;
      else if (cinfo[icell].type == INSIDE) c2f[icell] = 0;
      else if (cinfo[icell].type == OVERLAP) {
        c2f[icell] = 1;
        flags = cinfo[icell].corner;

        extflag = 0;
        for (j = 0; j < nface_pts; j++) {
          if (flags[corners[face][j]] == OUTSIDE) extflag = 1;
          else if (flags[corners[face][j]] == INSIDE) c2f[icell] = 0;
        }
        if (!extflag) c2f[icell] = 0;

        if (c2f[icell]) {
          if (dimension == 2) {
            for (j = 0; j < cells[icell].nsurf; j++) {
              n = cells[icell].csurfs[j];
              if (Geometry::
                  line_quad_face_touch(pts[lines[n].p1].x,
                                       pts[lines[n].p2].x,
                                       i,cells[icell].lo,cells[icell].hi)) {
                c2f[icell] = 0;
                break;
              }
            }
          } else {
            for (j = 0; j < cells[icell].nsurf; j++) {
              n = cells[icell].csurfs[j];
              if (Geometry::
                  tri_hex_face_touch(pts[tris[n].p1].x,
                                     pts[tris[n].p2].x,
                                     pts[tris[n].p3].x,
                                     i,cells[icell].lo,cells[icell].hi)) {
                c2f[icell] = 0;
                break;
              }
            }
          }
        }
      }
    } else c2f[icell] = 0;
  }

  // ncf = # of my cells to insert onto
  // some may be eliminated in interpolate() if no particles inserted

  ncf = 0;
  for (icell = 0; icell < nglocal; icell++)
    if (c2f[icell]) ncf++;

  // cellface = per-face data struct for all inserts performed on my grid cells
  // reallocate cellface since nspecies count of mixture may have changed

  int nspecies = particle->mixture[imix]->nspecies;

  for (int i = 0; i < ncfmax; i++) delete [] cellface[i].ntargetsp;
  memory->sfree(cellface);
  ncfmax = ncf;
  cellface = (CellFace *) memory->smalloc(ncfmax*sizeof(CellFace),
                                          "inflow/file:cellface");
  memset(cellface,0,ncfmax*sizeof(CellFace));

  for (i = 0; i < ncfmax; i++) {
    cellface[i].ntargetsp = new double[nspecies];
    cellface[i].fraction = new double[nspecies];
    cellface[i].cummulative = new double[nspecies];
    cellface[i].vscale = new double[nspecies];
  }

  // set id,pcell,icell for each cellface
  // convert c2f[I] = index into cellface, -1 if none
  // again, some may be eliminated in interpolate() if no particles inserted

  ncf = 0;
  for (icell = 0; icell < nglocal; icell++) {
    if (!c2f[icell]) {
      c2f[icell] = -1;
      continue;
    }
    cellface[ncf].id = cells[icell].id;
    if (cells[icell].nsplit > 1) cellface[ncf].pcell = split(icell,face);
    else cellface[ncf].pcell = icell;
    cellface[ncf].icell = icell;
    c2f[icell] = ncf++;
  }

  // check that any species listed as a mesh value is in mixture

  for (m = 0; m < mesh.nvalues; m++) {
    if (mesh.which[m] >= 0) continue;
    isp = -mesh.which[m] - 1;
    if (particle->mixture[imix]->species2species[isp] < 0)
      error->all(FLERR,"Fix inflow/file species is not in mixture");
  }

  // interpolate requested values from mesh to cell faces
  // ncf is decremented if cell face is outside mesh extent

  interpolate();

  // cummulative counter

  ntotal = 0;
}

/* ---------------------------------------------------------------------- */

void FixInflowFile::start_of_step()
{
  int pcell,ninsert,isp,ispecies,ndim,pdim1,pdim2,id;
  double *lo,*hi,*vstream,*cummulative,*vscale;
  double x[3],v[3];
  double indot,scosine,rn,ntarget,temp_thermal;
  double beta_un,normalized_distbn_fn,theta,erot,evib;
  Particle::OnePart *p;

  if (update->ntimestep % nevery) return;

  int dimension = domain->dimension;
  Particle::OnePart *particles = particle->particles;
  Grid::ChildCell *cells = grid->cells;
  double dt = update->dt;

  int nspecies = particle->mixture[imix]->nspecies;
  int *species = particle->mixture[imix]->species;

  // insert particles by cell face
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

  // NOTE:
  // currently not allowing particle insertion on backflow boundaries
  // enforced by indot >= 0.0 check in init()
  // could allow particle insertion on backflow boundaries
  //   when streaming velocity is small enough
  // need to insure two do-while loops below do not spin endlessly

  int nfix_add_particle = modify->n_add_particle;
  nsingle = 0;

  for (int i = 0; i < ncf; i++) {
    pcell = cellface[i].pcell;
    lo = cellface[i].lo;
    hi = cellface[i].hi;

    vstream = cellface[i].vstream;
    indot = vstream[0]*normal[0] + vstream[1]*normal[1] + vstream[2]*normal[2];
    temp_thermal = cellface[i].temp_thermal;
    vscale = cellface[i].vscale;

    if (perspecies == YES) {
      for (isp = 0; isp < nspecies; isp++) {
        ispecies = species[isp];
	ntarget = cellface[i].ntargetsp[isp]+random->uniform();
	ninsert = static_cast<int> (ntarget);
        scosine = indot / vscale[isp];

	for (int m = 0; m < ninsert; m++) {
	  x[0] = lo[0] + random->uniform() * (hi[0]-lo[0]);
	  x[1] = lo[1] + random->uniform() * (hi[1]-lo[1]);
	  if (dimension == 3) x[2] = lo[2] + random->uniform() * (hi[2]-lo[2]);
          else x[2] = 0.0;

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
          v[pdim] = vscale[isp]*sin(theta) + vstream[pdim];
          v[qdim] = vscale[isp]*cos(theta) + vstream[qdim];
          erot = particle->erot(ispecies,temp_thermal,random);
          evib = particle->evib(ispecies,temp_thermal,random);
          id = MAXSMALLINT*random->uniform();
	  particle->add_particle(id,ispecies,pcell,x,v,erot,evib);

          p = &particle->particles[particle->nlocal-1];
          p->flag = PINSERT;
          p->dtremain = dt * random->uniform();

          if (nfix_add_particle) 
            modify->add_particle(particle->nlocal-1,temp_thermal,vstream);
	}

	nsingle += ninsert;
      }

    } else {
      cummulative = cellface[i].cummulative;
      ntarget = cellface[i].ntarget+random->uniform();
      ninsert = static_cast<int> (ntarget);

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
        v[pdim] = vscale[isp]*sin(theta) + vstream[pdim];
        v[qdim] = vscale[isp]*cos(theta) + vstream[qdim];

        erot = particle->erot(ispecies,temp_thermal,random);
        evib = particle->evib(ispecies,temp_thermal,random);
        id = MAXSMALLINT*random->uniform();
	particle->add_particle(id,ispecies,pcell,x,v,erot,evib);

        p = &particle->particles[particle->nlocal-1];
        p->flag = PINSERT;
        p->dtremain = dt * random->uniform();

        if (nfix_add_particle) 
          modify->add_particle(particle->nlocal-1,temp_thermal,vstream);
      }

      nsingle += ninsert;
    }
  }

  ntotal += nsingle;
}

/* ----------------------------------------------------------------------
   scan file for section-ID, read regular grid of values into Mesh data struct
   only called by proc 0
------------------------------------------------------------------------- */

void FixInflowFile::read_file(char *file, char *section)
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

void FixInflowFile::bcast_mesh()
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

void FixInflowFile::check_mesh_values()
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
   interpolate values needed be each cell face
   use centroid of face to find 4 surrounding mesh points (2 in 2d)
   use bilinear interpolation (linear in 2d) to compute each new value
   reset ncf to reflect faces that overlap and insert particles
------------------------------------------------------------------------- */

void FixInflowFile::interpolate()
{
  int i,j,m,icell,isp,plo,phi,qlo,qhi,anyfrac,err;
  double indot,newtemp;
  double *lo,*hi;
  double xc[2];

  int dimension = domain->dimension;
  double *meshlo = mesh.lo;
  double *meshhi = mesh.hi;
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  int nspecies = particle->mixture[imix]->nspecies;
  double nrho = particle->mixture[imix]->nrho;
  double *vstream = particle->mixture[imix]->vstream;
  double *vscale = particle->mixture[imix]->vscale;
  double *fraction = particle->mixture[imix]->fraction;
  int *fraction_flag = particle->mixture[imix]->fraction_flag;
  double *fraction_user = particle->mixture[imix]->fraction_user;
  double *cummulative = particle->mixture[imix]->cummulative;
  double temp_thermal = particle->mixture[imix]->temp_thermal;
  int *species2species = particle->mixture[imix]->species2species;
  double fnum = update->fnum;
  double dt = update->dt;

  // per-species vectors for mesh setting of species fractions
  // initialize to mixture settings

  int *fflag = new int[nspecies];
  double *fuser = new double[nspecies];
  for (isp = 0; isp < nspecies; isp++) {
    fflag[isp] = fraction_flag[isp];
    fuser[isp] = fraction_user[isp];
  }

  // loop over all possible faces
  // remove from list if face is outside mesh extent, resetting c2f

  int ninitial = ncf;
  ncf = 0;

  for (i = 0; i < ninitial; i++) {
    cellface[ncf].id = cellface[i].id;
    cellface[ncf].pcell = cellface[i].pcell;
    icell = cellface[ncf].icell = cellface[i].icell;

    c2f[icell] = -1;

    // cellface lo/hi = extent of overlap
    // if lo >= hi, then no overlap

    lo = cells[icell].lo;
    hi = cells[icell].hi;
    
    if (face % 2 == 0) 
      cellface[ncf].lo[ndim] = cellface[ncf].hi[ndim] = lo[ndim];
    else 
      cellface[ncf].lo[ndim] = cellface[ncf].hi[ndim] = hi[ndim];

    cellface[ncf].lo[pdim] = MAX(lo[pdim],meshlo[0]);
    cellface[ncf].hi[pdim] = MIN(hi[pdim],meshhi[0]);
    if (cellface[ncf].lo[pdim] >= cellface[ncf].hi[pdim]) continue;

    if (dimension == 3) {
      cellface[ncf].lo[qdim] = MAX(lo[qdim],meshlo[1]);
      cellface[ncf].hi[qdim] = MIN(hi[qdim],meshhi[1]);
      if (cellface[ncf].lo[qdim] >= cellface[ncf].hi[qdim]) continue;
    }
    
    // cellface defaults from mixture

    cellface[ncf].nrho = nrho;
    cellface[ncf].vstream[0] = vstream[0];
    cellface[ncf].vstream[1] = vstream[1];
    cellface[ncf].vstream[2] = vstream[2];
    cellface[ncf].temp_thermal = temp_thermal;
    for (j = 0; j < nspecies; j++) {
      cellface[ncf].fraction[j] = fraction[j];
      cellface[ncf].cummulative[j] = cummulative[j];
      cellface[ncf].vscale[j] = vscale[j];
    }

    // xc = centroid of cell face overlap with mesh

    xc[0] = 0.5 * (cellface[ncf].lo[pdim]+cellface[ncf].hi[pdim]);
    if (dimension == 3)
      xc[1] = 0.5 * (cellface[ncf].lo[qdim]+cellface[ncf].hi[qdim]);

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

    // override cellface defaults with file values
    // 2d or 3d = linear or bilinear interpolation
    // if temp_thermal is set, vscale must also be reset

    anyfrac = 0;

    if (dimension == 2) {
      for (m = 0; m < mesh.nvalues; m++) {
        if (mesh.which[m] == NRHO)
          cellface[ncf].nrho = linear_interpolation(xc[0],m,plo,phi);
        else if (mesh.which[m] == TEMP_THERMAL) {
          newtemp = cellface[ncf].temp_thermal = 
            linear_interpolation(xc[0],m,plo,phi);
          for (isp = 0; isp < nspecies; isp++)
            cellface[ncf].vscale[isp] = 
              vscale[isp] * sqrt(newtemp/temp_thermal);
        } else if (mesh.which[m] == VX)
          cellface[ncf].vstream[0] = linear_interpolation(xc[0],m,plo,phi);
        else if (mesh.which[m] == VY)
          cellface[ncf].vstream[1] = linear_interpolation(xc[0],m,plo,phi);
        else if (mesh.which[m] == VZ)
          cellface[ncf].vstream[2] = linear_interpolation(xc[0],m,plo,phi);
        else {
          anyfrac = 1;
          isp = -mesh.which[m] - 1;
          isp = species2species[isp];
          fuser[isp] = linear_interpolation(xc[0],m,plo,phi);
          fflag[isp] = 1;
        }
      }
    } else {
      for (m = 0; m < mesh.nvalues; m++) {
        if (mesh.which[m] == NRHO)
          cellface[ncf].nrho = 
            bilinear_interpolation(xc[0],xc[1],m,plo,phi,qlo,qhi);
        else if (mesh.which[m] == TEMP_THERMAL) {
          newtemp = cellface[ncf].temp_thermal = 
            bilinear_interpolation(xc[0],xc[1],m,plo,phi,qlo,qhi);
          for (isp = 0; isp < nspecies; isp++)
            cellface[ncf].vscale[isp] = 
              vscale[isp] * sqrt(newtemp/temp_thermal);
        } else if (mesh.which[m] == VX)
          cellface[ncf].vstream[0] =
            bilinear_interpolation(xc[0],xc[1],m,plo,phi,qlo,qhi);
        else if (mesh.which[m] == VY)
          cellface[ncf].vstream[1] =
            bilinear_interpolation(xc[0],xc[1],m,plo,phi,qlo,qhi);
        else if (mesh.which[m] == VZ)
          cellface[ncf].vstream[2] =
            bilinear_interpolation(xc[0],xc[1],m,plo,phi,qlo,qhi);
        else {
          anyfrac = 1;
          isp = -mesh.which[m] - 1;
          isp = species2species[isp];
          fuser[isp] = bilinear_interpolation(xc[0],xc[1],m,plo,phi,qlo,qhi);
          fflag[isp] = 1;
        }
      }
    }

    // if any species fraction was set,
    // set entire fraction and cummulative vector via Mixture::init_fraction()

    if (anyfrac) {
      err = particle->mixture[imix]->
        init_fraction(fflag,fuser,
                      cellface[ncf].fraction,cellface[ncf].cummulative);
      if (err)
        error->one(FLERR,"Fix inflow/file mixture fractions exceed 1.0");
    }

    // indot = dot product of vstream with inward face normal
    // skip cellface if indot < 0.0, to not allow any particles to be inserted
    // area of cell-face and mesh overlap depends on 2d/3d
    // also scale ntarget by frac_user (0 to 1)

    cellface[ncf].ntarget = 0.0;
    indot = cellface[ncf].vstream[0]*normal[0] +
      cellface[ncf].vstream[1]*normal[1] +
      cellface[ncf].vstream[2]*normal[2];
    if (indot >= 0.0) {
      double area = (cellface[ncf].hi[pdim]-cellface[ncf].lo[pdim]);
      if (dimension == 3) 
        area *= (cellface[ncf].hi[qdim]-cellface[ncf].lo[qdim]);
      for (isp = 0; isp < nspecies; isp++) {
        cellface[ncf].ntargetsp[isp] = frac_user *
          mol_inflow(indot,cellface[ncf].vscale[isp],
                     cellface[ncf].fraction[isp]);
        cellface[ncf].ntargetsp[isp] *= cellface[ncf].nrho*area*dt / fnum;
        cellface[ncf].ntargetsp[isp] /= cinfo[icell].weight;
        cellface[ncf].ntarget += cellface[ncf].ntargetsp[isp];
      }
      c2f[icell] = ncf++;
    }
  }

  // clean up

  delete [] fflag;
  delete [] fuser;
}

/* ----------------------------------------------------------------------
   linear interpolation at x between lo and hi bounds, for column M
------------------------------------------------------------------------- */

double FixInflowFile::linear_interpolation(double x, int m, int plo, int phi)
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

double FixInflowFile::bilinear_interpolation(double x, double y, int m, 
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
   calculate flux of particles of species ISP entering a grid cell
   indot = vstream dotted into face normal, assumed to be >= 0.0
   scosine = s cos(theta) in Bird notation where vscale = 1/beta
   see Bird 1994, eq 4.22
------------------------------------------------------------------------- */

double FixInflowFile::mol_inflow(double indot, double vscale, double fraction)
{
  double scosine = indot / vscale;
  double inward_number_flux = vscale*fraction *
    (exp(-scosine*scosine) + sqrt(MY_PI)*scosine*(1.0 + erf(scosine))) / 
    (2*sqrt(MY_PI));
  return inward_number_flux;
}

/* ----------------------------------------------------------------------
   inserting into split cell parent ICELL on face FLAG
   determine which child split cell the face is part of
   face cannot be touched by surfs, so entire face is part of one split cell
   compute which via update->split() and return it
------------------------------------------------------------------------- */

int FixInflowFile::split(int icell, int flag)
{
  double x[3];

  Grid::ChildCell *cells = grid->cells;

  // x = center point on face

  x[0] = 0.5 * (cells[icell].lo[0] + cells[icell].hi[0]);
  x[1] = 0.5 * (cells[icell].lo[1] + cells[icell].hi[1]);
  x[2] = 0.5 * (cells[icell].lo[2] + cells[icell].hi[2]);
  if (flag == XLO) x[0] = cells[icell].lo[0];
  if (flag == XHI) x[0] = cells[icell].hi[0];
  if (flag == YLO) x[1] = cells[icell].lo[1];
  if (flag == YHI) x[1] = cells[icell].hi[1];
  if (flag == ZLO) x[2] = cells[icell].lo[2];
  if (flag == ZHI) x[2] = cells[icell].hi[2];
  if (domain->dimension == 2) x[2] = 0.0;

  int splitcell;
  if (domain->dimension == 2) splitcell = update->split2d(icell,x);
  else splitcell = update->split3d(icell,x);
  return splitcell;
}

/* ----------------------------------------------------------------------
   pack icell values for per-cell arrays into buf
   also pack cellface data for flagged faces
   return byte count of amount packed
   if memflag, only return count, do not fill buf
------------------------------------------------------------------------- */

int FixInflowFile::pack_grid_one(int icell, char *buf, int memflag)
{
  char *ptr = buf;

  int nspecies = particle->mixture[imix]->nspecies;

  if (memflag) memcpy(ptr,&c2f[icell],sizeof(int));
  ptr += sizeof(int);
  ptr = ROUNDUP(ptr);

  // pack cellface entry if it exists, plus its vectors

  if (c2f[icell] >= 0) {
    int icf = c2f[icell];
    if (memflag) memcpy(ptr,&cellface[icf],sizeof(CellFace));
    ptr += sizeof(CellFace);
    ptr = ROUNDUP(ptr);
    if (memflag) memcpy(ptr,cellface[icf].ntargetsp,nspecies*sizeof(double));
    ptr += nspecies*sizeof(double);
    if (memflag) memcpy(ptr,cellface[icf].fraction,nspecies*sizeof(double));
    ptr += nspecies*sizeof(double);
    if (memflag) memcpy(ptr,cellface[icf].cummulative,nspecies*sizeof(double));
    ptr += nspecies*sizeof(double);
    if (memflag) memcpy(ptr,cellface[icf].vscale,nspecies*sizeof(double));
    ptr += nspecies*sizeof(double);
  }

  return ptr-buf;
}

/* ----------------------------------------------------------------------
   unpack icell values for per-cell arrays from buf
   also unpack cellface data for flagged faces
   return byte count of amount unpacked
------------------------------------------------------------------------- */

int FixInflowFile::unpack_grid_one(int icell, char *buf)
{
  char *ptr = buf;


  int nspecies = particle->mixture[imix]->nspecies;

  grow_percell(1);
  memcpy(&c2f[icell],ptr,sizeof(int));
  ptr += sizeof(int);
  ptr = ROUNDUP(ptr);
  nglocal++;

  // unpack c2f values
  // fill sub cells with -1

  int nsplit = grid->cells[icell].nsplit;
  if (nsplit > 1) {
    grow_percell(nsplit);
    for (int i = 0; i < nsplit; i++) c2f[nglocal++] = -1;
  }
  
  // unpack cellface entry for c2f
  // store 4 vector ptrs to avoid overwriting allocated cellface vectors
  // reset c2f pointer into cellface
  // reset cellface.pcell and cellface.icell
  // pcell setting based on sub cells being immediately after the split cell

  if (c2f[icell] >= 0) {
    c2f[icell] = ncf;
    grow_cellface(1);

    double *ntargetsp = cellface[ncf].ntargetsp;
    double *fraction = cellface[ncf].fraction;
    double *cummulative = cellface[ncf].cummulative;
    double *vscale = cellface[ncf].vscale;

    memcpy(&cellface[ncf],ptr,sizeof(CellFace));
    ptr += sizeof(CellFace);
    ptr = ROUNDUP(ptr);

    cellface[ncf].ntargetsp = ntargetsp;
    cellface[ncf].fraction = fraction;
    cellface[ncf].cummulative = cummulative;
    cellface[ncf].vscale = vscale;
    memcpy(cellface[ncf].ntargetsp,ptr,nspecies*sizeof(double));
    ptr += nspecies*sizeof(double);
    memcpy(cellface[ncf].fraction,ptr,nspecies*sizeof(double));
    ptr += nspecies*sizeof(double);
    memcpy(cellface[ncf].cummulative,ptr,nspecies*sizeof(double));
    ptr += nspecies*sizeof(double);
    memcpy(cellface[ncf].vscale,ptr,nspecies*sizeof(double));
    ptr += nspecies*sizeof(double);

    if (grid->cells[icell].nsplit == 1) cellface[ncf].pcell = icell;
    else cellface[ncf].pcell = split(icell,face);
    cellface[ncf].icell = icell;
    ncf++;
  }

  return ptr-buf;
}

/* ----------------------------------------------------------------------
   compress per-cell arrays due to cells migrating to new procs
   criteria for keeping/discarding a cell is same as in Grid::compress()
   this keeps final ordering of per-cell arrays consistent with Grid class
------------------------------------------------------------------------- */

void FixInflowFile::compress_grid()
{
  int me = comm->me;
  Grid::ChildCell *cells = grid->cells;

  // compress cellface first, preserving/copying 4 cellface vectors
  // reset c2f pointers to new cellface indices

  int nspecies = particle->mixture[imix]->nspecies;
  double *ntargetsp,*fraction,*cummulative,*vscale;

  int ncurrent = ncf;
  ncf = 0;
  for (int icf = 0; icf < ncurrent; icf++) {
    int icell = cellface[icf].icell;
    if (cells[icell].proc != me) continue;
    if (ncf != icf) {
      ntargetsp = cellface[ncf].ntargetsp;
      fraction = cellface[ncf].fraction;
      cummulative = cellface[ncf].cummulative;
      vscale = cellface[ncf].vscale;
      memcpy(&cellface[ncf],&cellface[icf],sizeof(CellFace));
      cellface[ncf].ntargetsp = ntargetsp;
      cellface[ncf].fraction = fraction;
      cellface[ncf].cummulative = cummulative;
      cellface[ncf].vscale = vscale;
      memcpy(ntargetsp,cellface[icf].ntargetsp,nspecies*sizeof(double));
      memcpy(fraction,cellface[icf].fraction,nspecies*sizeof(double));
      memcpy(cummulative,cellface[icf].cummulative,nspecies*sizeof(double));
      memcpy(vscale,cellface[icf].vscale,nspecies*sizeof(double));
    }
    c2f[icell] = ncf;
    ncf++;
  }

  // compress c2f second
  // keep an unsplit or split cell if staying on this proc
  // keep a sub cell if its split cell is staying on this proc
  // reset cellface.icell to new cell index
  // cellface.pcell reset in post_compress(),
  //   since don't know sub cell index at this point

  ncurrent = nglocal;
  nglocal = 0;
  for (int icell = 0; icell < ncurrent; icell++) {
    if (cells[icell].nsplit >= 1) {
      if (cells[icell].proc != me) continue;
    } else {
      int isplit = cells[icell].isplit;
      if (cells[grid->sinfo[isplit].icell].proc != me) continue;
    }

    if (nglocal != icell) c2f[nglocal] = c2f[icell];
    if (c2f[nglocal] >= 0) cellface[c2f[nglocal]].icell = nglocal;
    nglocal++;
  }
}

/* ----------------------------------------------------------------------
   reset cellface.pcell for compressed cellface entries
   called from Grid::compress() after grid cells have been compressed
------------------------------------------------------------------------- */

void FixInflowFile::post_compress_grid()
{
  Grid::ChildCell *cells = grid->cells;

  for (int icf = 0; icf < ncf; icf++) {
    int icell = cellface[icf].icell;
    int nsplit = cells[icell].nsplit;
    if (nsplit == 1) cellface[icf].pcell = icell;
    else cellface[icf].pcell = split(icell,face);
  }
}

/* ----------------------------------------------------------------------
   insure c2f allocated long enough for N new cells
------------------------------------------------------------------------- */

void FixInflowFile::grow_percell(int n)
{
  if (nglocal+n < nglocalmax) return;
  nglocalmax += DELTAGRID;
  memory->grow(c2f,nglocalmax,"inflow/file:c2f");
}

/* ----------------------------------------------------------------------
   insure cellface allocated long enough for N new cell faces
   also allocate new vectors within each cellface
------------------------------------------------------------------------- */

void FixInflowFile::grow_cellface(int n)
{
  if (ncf+n < ncfmax) return;
  int oldmax = ncfmax;
  ncfmax += DELTAFACE;
  cellface = (CellFace *) memory->srealloc(cellface,ncfmax*sizeof(CellFace),
                                           "inflow:cellface");
  memset(&cellface[oldmax],0,(ncfmax-oldmax)*sizeof(CellFace));

  int nspecies = particle->mixture[imix]->nspecies;
  for (int i = oldmax; i < ncfmax; i++) {
    cellface[i].ntargetsp = new double[nspecies];
    cellface[i].fraction = new double[nspecies];
    cellface[i].cummulative = new double[nspecies];
    cellface[i].vscale = new double[nspecies];
  }
}

/* ----------------------------------------------------------------------
   return one-step or total count of particle insertions
------------------------------------------------------------------------- */

double FixInflowFile::compute_vector(int i)
{
  double one,all;
  
  if (i == 0) one = nsingle;
  else one = ntotal;
  MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}

/* ----------------------------------------------------------------------
   DEBUG method: print status of cellface I
------------------------------------------------------------------------- */

void FixInflowFile::print_face(int i)
{
  printf("CFACE: proc %d index %d id %d pcell %d icell %d\n",
         comm->me,i,cellface[i].id,cellface[i].pcell,cellface[i].icell);
  printf("CFACE: proc %d index %d lo %g %g %g hi %g %g %g\n",
         comm->me,i,
         cellface[i].lo[0],cellface[i].lo[1],cellface[i].lo[2],
         cellface[i].hi[0],cellface[i].hi[1],cellface[i].hi[2]);
  printf("CFACE: proc %d index %d ntarget %g %g %g nrho %g vstream %g %g %g\n",
         comm->me,i,cellface[i].
         ntarget,cellface[i].ntargetsp[0],cellface[i].ntargetsp[1],
         cellface[i].nrho,
         cellface[i].vstream[0],cellface[i].vstream[1],cellface[i].vstream[2]);
}
