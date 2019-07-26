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

/* ----------------------------------------------------------------------
   Contributing author: Arnaud Borner (NASA Ames)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "read_isurf.h"
#include "surf.h"
#include "domain.h"
#include "grid.h"
#include "comm.h"
#include "modify.h"
#include "fix_ablate.h"
#include "input.h"
#include "marching_squares.h"
#include "marching_cubes.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{NEITHER,BAD,GOOD};
enum{NONE,CHECK,KEEP};
enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};           // several files
enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};         // same as Domain

#define CHUNK 8192

/* ---------------------------------------------------------------------- */

ReadISurf::ReadISurf(SPARTA *sparta) : Pointers(sparta)
{
  MPI_Comm_rank(world,&me);
  dim = domain->dimension;

  cvalues = NULL;
  svalues = NULL;
}

/* ---------------------------------------------------------------------- */

ReadISurf::~ReadISurf()
{
  memory->destroy(cvalues);
  memory->destroy(svalues);
}

/* ---------------------------------------------------------------------- */

void ReadISurf::command(int narg, char **arg)
{
  // NOTE: at some point could allow another chunk of isurfs to be read
  //       see note below about clear_surf()

  if (!grid->exist)
    error->all(FLERR,"Cannot read_isurf before grid is defined");
  if (!surf->implicit)
    error->all(FLERR,"Cannot read_isurf unless global surfs implicit is set");
  if (surf->exist)
    error->all(FLERR,"Cannot read_isurf when surfs already exist");
  if (particle->exist)
    error->all(FLERR,"Cannot read_isurf when particles exist");
  if (domain->axisymmetric)
    error->all(FLERR,"Cannot read_isurf for axisymmetric domains");

  surf->exist = 1;

  if (narg < 6) error->all(FLERR,"Illegal read_isurf command");

  ggroup = grid->find_group(arg[0]);
  if (ggroup < 0) error->all(FLERR,"Read_isurf grid group ID does not exist");

  nx = input->inumeric(FLERR,arg[1]);
  ny = input->inumeric(FLERR,arg[2]);
  nz = input->inumeric(FLERR,arg[3]);

  if (dim == 2 && nz != 1) error->all(FLERR,"Invalid read_isurf command");

  char *gridfile = arg[4];

  thresh = input->numeric(FLERR,arg[5]);
  if (thresh <= 0 || thresh >= 255) 
    error->all(FLERR,"Invalid read_isurf command");
  int ithresh = static_cast<int> (thresh);
  if (ithresh == thresh) 
    error->all(FLERR,"An integer value for read_isurf thresh is not allowed");

  // process command line args

  process_args(narg-6,&arg[6]);

  // verify that grid group is a set of uniform child cells
  // must comprise a 3d contiguous block

  int nxyz[3];
  int count = grid->check_uniform_group(ggroup,nxyz,corner,xyzsize);
  if (nx != nxyz[0] || ny != nxyz[1] || nz != nxyz[2])
    error->all(FLERR,"Read_isurf grid group does not match nx,ny,nz");

  // read grid corner point values
  // create and destroy dictionary of my grid cells in group
  //   used to assign per-grid values to local grid cells

  if (screen && me == 0) fprintf(screen,"Reading isurf file ...\n");

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  create_hash(count);

  if (dim == 2) {
    memory->create(cvalues,grid->nlocal,4,"readisurf:cvalues");
    for (int i = 0; i < grid->nlocal; i++)
      for (int j = 0; j < 4; j++)
        cvalues[i][j] = 0.0;
  } else {
    memory->create(cvalues,grid->nlocal,8,"readisurf:cvalues");
    for (int i = 0; i < grid->nlocal; i++)
      for (int j = 0; j < 8; j++)
        cvalues[i][j] = 0.0;
  }

  read_corners(gridfile);

  if (typefile) {
    memory->create(svalues,grid->nlocal,"readisurf:svalues");
    read_types(typefile);
  }

  delete hash;

  // create surfs in each grid cell based on corner point values
  // call clear_surf first so cell/corner flags are all set
  // set surf->nsurf and surf->nown
  // if specified, apply group keyword to reset per-surf mask info

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  grid->unset_neighbors();
  grid->remove_ghosts();
  grid->clear_surf();

  MarchingSquares *ms;
  MarchingCubes *mc;

  if (dim == 2) {
    ms = new MarchingSquares(sparta,ggroup,thresh);
    ms->invoke(cvalues,svalues);
  } else {
    mc = new MarchingCubes(sparta,ggroup,thresh);
    mc->invoke(cvalues,svalues);
  }

  surf->nown = surf->nlocal;
  bigint nlocal = surf->nlocal;
  MPI_Allreduce(&nlocal,&surf->nsurf,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  if (sgrouparg) {
    int sgroup = surf->find_group(arg[sgrouparg]);
    if (sgroup < 0) sgroup = surf->add_group(arg[sgrouparg]);
    int sgroupbit = surf->bitmask[sgroup];
 
    int nsurf = surf->nlocal;
    if (dim == 3) {
      Surf::Tri *tris = surf->tris;
      for (int i = 0; i < nsurf; i++) tris[i].mask |= sgroupbit;
    } else {
      Surf::Line *lines = surf->lines;
      for (int i = 0; i < nsurf; i++) lines[i].mask |= sgroupbit;
    }
  }

  // output extent of implicit surfs, some may be tiny

  if (dim == 2) surf->output_extent(0);
  else surf->output_extent(0);

  MPI_Barrier(world);
  double time3 = MPI_Wtime();

  // compute normals of new surfs

  if (dim == 2) surf->compute_line_normal(0);
  else surf->compute_tri_normal(0);

  // cleanup() checks for consistent triangles on grid cell faces
  // needs to come after normals are computed

  if (dim == 2) delete ms;
  else {
    mc->cleanup();
    delete mc;
  }

  // watertight check can be done before surfs are mapped to grid cells

  if (dim == 2) surf->check_watertight_2d();
  else surf->check_watertight_3d();

  MPI_Barrier(world);
  double time4 = MPI_Wtime();

  // -----------------------
  // map surfs to grid cells
  // -----------------------

  // surfs are already assigned to grid cells
  // create split cells due to new surfs

  grid->surf2grid_implicit(1,1);

  MPI_Barrier(world);
  double time5 = MPI_Wtime();

  // re-setup grid ghosts and neighbors

  grid->setup_owned();
  grid->acquire_ghosts();
  grid->reset_neighbors();
  comm->reset_neighbors();

  MPI_Barrier(world);
  double time6 = MPI_Wtime();

  // flag cells and corners as OUTSIDE or INSIDE

  grid->set_inout();
  grid->type_check();

  MPI_Barrier(world);
  double time7 = MPI_Wtime();

  // store corner point cvalues and type values in FixAblate if requested
  // also pass it the geometry of the 3d grid of cells
  // cannot do until now, b/c need split cells to exist via surf2grid()

  if (storeflag) ablate->store_corners(nx,ny,nz,corner,xyzsize,
                                       cvalues,svalues,thresh);

  MPI_Barrier(world);
  double time8 = MPI_Wtime();

  // stats

  double time_total = time8-time1;
  double time_s2g = time5-time4;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  CPU time = %g secs\n",time_total);
      fprintf(screen,"  read/marching/check/surf2grid/ghost/inout/store "
              "percent = %g %g %g %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total,
              100.0*(time6-time5)/time_total,100.0*(time7-time6)/time_total,
              100.0*(time8-time7)/time_total);
      fprintf(screen,"  surf2grid time = %g secs\n",time_s2g);
      fprintf(screen,"  map/rvous/split percent = %g %g %g\n",
              100.0*grid->tmap/time_s2g,100.0*grid->trvous1/time_s2g,
              100.0*grid->tsplit/time_s2g);
    }
    if (logfile) {
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
      fprintf(logfile,"  read/marching/check/surf2grid/ghost/inout/store "
              "percent = %g %g %g %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total,
              100.0*(time6-time5)/time_total,100.0*(time7-time6)/time_total,
              100.0*(time8-time7)/time_total);
      fprintf(logfile,"  surf2grid time = %g secs\n",time_s2g);
      fprintf(logfile,"  map/rvous/split percent = %g %g %g\n",
              100.0*grid->tmap/time_s2g,100.0*grid->trvous1/time_s2g,
              100.0*grid->tsplit/time_s2g);
    }
  }
}

/* ----------------------------------------------------------------------
   read/store all grid corner point values
------------------------------------------------------------------------- */

void ReadISurf::read_corners(char *gridfile)
{
  int nchunk;
  int nxyz[3];
  FILE *fp;

  uint8_t *buf;
  memory->create(buf,CHUNK,"readisurf:buf");

  // proc 0 opens and reads binary file
  // error check the file grid matches input script extent

  if (me == 0) {
    fp = fopen(gridfile,"rb");
    if (fp == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open read_isurf grid corner point file %s",
               gridfile);
      error->one(FLERR,str);
    }
    fread(nxyz,sizeof(int),dim,fp);
  }

  MPI_Bcast(nxyz,dim,MPI_INT,0,world);

  int flag = 0;
  if (nxyz[0] != nx+1) flag = 1;
  if (nxyz[1] != ny+1) flag = 1;
  if (dim == 3 && nxyz[2] != nz+1) flag = 1;
  if (flag) 
    error->all(FLERR,"Grid size in read_isurf grid corner point file "
               "does not match request");

  // read and broadcast one CHUNK of values at a time
  // each proc stores grid corner point values it needs in assign_corners()

  bigint ncorners;
  if (dim == 3) ncorners = (bigint) (nx+1) * (ny+1)*(nz+1);
  else ncorners = (bigint) (nx+1) * (ny+1)*nz;

  bigint nread = 0;

  while (nread < ncorners) {
    if (ncorners-nread > CHUNK) nchunk = CHUNK;
    else nchunk = ncorners-nread;

    if (me == 0) fread(buf,sizeof(uint8_t),nchunk,fp);
    MPI_Bcast(buf,nchunk,MPI_CHAR,0,world);

    assign_corners(nchunk,nread,buf);
    nread += nchunk;
  }

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " corner points\n",ncorners);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " corner points\n",ncorners);
  }

  memory->destroy(buf);

  // close file

  if (me == 0) fclose(fp);
}

/* ----------------------------------------------------------------------
   read/store all grid surface type values
------------------------------------------------------------------------- */

void ReadISurf::read_types(char *typefile)
{
  int nchunk;
  int nxyz[3];
  FILE *fp;

  int *buf;
  memory->create(buf,CHUNK,"readisurf:buf");

  // proc 0 opens and reads binary file
  // error check the file grid matches input script extent

  if (me == 0) {
    fp = fopen(typefile,"rb");
    if (fp == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open read_isurf type file %s",typefile);
      error->one(FLERR,str);
    }
    fread(nxyz,sizeof(int),dim,fp);
  }

  MPI_Bcast(nxyz,dim,MPI_INT,0,world);

  int flag = 0;
  if (nxyz[0] != nx) flag = 1;
  if (nxyz[1] != ny) flag = 1;
  if (dim == 3 && nxyz[2] != nz) flag = 1;
  if (flag) 
    error->all(FLERR,"Grid size in read_isurf type file does not match request");

  // read and broadcast one CHUNK of values at a time
  // each proc stores grid corner point values it needs in assign_corners()

  bigint ntypes = (bigint) nx * ny*nz;
  bigint nread = 0;

  while (nread < ntypes) {
    if (ntypes-nread > CHUNK) nchunk = CHUNK;
    else nchunk = ntypes-nread;

    if (me == 0) fread(buf,sizeof(int),nchunk,fp);
    MPI_Bcast(buf,nchunk,MPI_INT,0,world);

    assign_types(nchunk,nread,buf);
    nread += nchunk;
  }

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " surface types\n",ntypes);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " surface types\n",ntypes);
  }

  memory->destroy(buf);

  // close file

  if (me == 0) fclose(fp);
}

/* ----------------------------------------------------------------------
   create hash for my grid cells in group
   key = index (0 to N-1) of grid cell in Nx by Ny by Nz contiguous block
   value = my local icell
   NOTE: could use count to prealloc the hash size
------------------------------------------------------------------------- */

void ReadISurf::create_hash(int count)
{
  hash = new MyHash;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;
  int groupbit = grid->bitmask[ggroup];

  int ix,iy,iz;
  bigint index;

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    ix = static_cast<int> ((cells[icell].lo[0]-corner[0]) / xyzsize[0] + 0.5);
    iy = static_cast<int> ((cells[icell].lo[1]-corner[1]) / xyzsize[1] + 0.5);
    iz = static_cast<int> ((cells[icell].lo[2]-corner[2]) / xyzsize[2] + 0.5);
    index = (bigint) nx * ny*iz + nx*iy + ix;
    (*hash)[index] = icell;
  }
}

/* ----------------------------------------------------------------------
   store all grid corner point values
   use hash to see if I own any grid cells that contain a corner point
   each corner point can be stored by as many as 4 or 8 grid cells
   check that corner point values = 0 on boundary of grid block
------------------------------------------------------------------------- */

void ReadISurf::assign_corners(int n, bigint offset, uint8_t *buf)
{
  int icell,ncorner,zeroflag;
  int pix,piy,piz;
  bigint pointindex,cellindex;

  for (int i = 0; i < n; i++) {
    pointindex = offset + i;
    pix = pointindex % (nx+1);
    piy = (pointindex / (nx+1)) % (ny+1);
    piz = pointindex / ((nx+1)*(ny+1));

    zeroflag = 0;
    if (buf[i]) {
      if (pix == 0 || piy == 0) zeroflag = 1;
      if (pix == nx || piy == ny) zeroflag = 1;
      if (dim == 3 && (piz == 0 || piz == nz)) zeroflag = 1;
      if (zeroflag) error->all(FLERR,"Grid boundary value != 0");
    }

    // ncorner = 0,1,2,3,4,5,6,7 when corner point is 
    //   bottom-lower-left, bottom-lower-right, 
    //   bottom-upper-left, bottom-upper-right,
    //   top-lower-left, top-lower-right, top-upper-left, top-upper-right
    //   of cell
    // if test on cix,ciy,ciz excludes cells that are outside of grid block

    if (dim == 3) {
      ncorner = 8;
      for (int ciz = piz-1; ciz <= piz; ciz++) {
        for (int ciy = piy-1; ciy <= piy; ciy++) {
          for (int cix = pix-1; cix <= pix; cix++) {
            ncorner--;
            if (cix < 0 || cix >= nx || ciy < 0 || ciy >=ny || 
                ciz < 0 || ciz >= nz) continue;
            cellindex = (bigint) nx * ny*ciz + nx*ciy + cix;
            if (hash->find(cellindex) == hash->end()) continue;
            icell = (*hash)[cellindex];
            cvalues[icell][ncorner] = buf[i];
          }
        }
      }

    // ncorner = 0,1,2,3 when corner point is 
    //   lower-left, lower-right, upper-left, upper-right of cell
    // if test on cix,ciy excludes cells that are outside of grid block

    } else {
      ncorner = 4;
      for (int ciy = piy-1; ciy <= piy; ciy++) {
        for (int cix = pix-1; cix <= pix; cix++) {
          ncorner--;
          if (cix < 0 || cix >= nx || ciy < 0 || ciy >=ny) continue;
          cellindex = (bigint) nx * ciy + cix;
          if (hash->find(cellindex) == hash->end()) continue;
          icell = (*hash)[cellindex];
          cvalues[icell][ncorner] = buf[i];
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   store all grid surf type values
   use hash to see if I own grid cell corresponding to index (0 to N-1)
------------------------------------------------------------------------- */

void ReadISurf::assign_types(int n, bigint offset, int *buf)
{
  int icell;
  bigint cellindex;

  for (int i = 0; i < n; i++) {
    cellindex = offset + i;
    if (hash->find(cellindex) == hash->end()) continue;
    icell = (*hash)[cellindex];
    svalues[icell] = buf[i];
  }
}

/* ----------------------------------------------------------------------
   process command line args
------------------------------------------------------------------------- */

void ReadISurf::process_args(int narg, char **arg)
{
  sgrouparg = 0;
  typefile = NULL;
  storeflag = 0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"group") == 0)  {
      if (iarg+2 > narg) error->all(FLERR,"Invalid read_isurf command");
      sgrouparg = iarg+1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"type") == 0)  {
      if (iarg+2 > narg) error->all(FLERR,"Invalid read_isurf command");
      typefile = arg[iarg+1];
      iarg += 2;
    } else if (strcmp(arg[iarg],"store") == 0)  {
      if (iarg+2 > narg) error->all(FLERR,"Invalid read_isurf command");
      storeflag = 1;
      storeID = arg[iarg+1];
      int ifix = modify->find_fix(storeID);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for read_isurf store does not exist");
      if (strcmp(modify->fix[ifix]->style,"ablate") != 0)
        error->all(FLERR,"Fix for read_isurf store is not a fix ablate");
      ablate = (FixAblate *) modify->fix[ifix];
      if (ggroup != ablate->igroup)
        error->all(FLERR,"Read_isurf group does not match fix ablate group");
      iarg += 2;
    } else error->all(FLERR,"Invalid read_isurf command");
  }
}

