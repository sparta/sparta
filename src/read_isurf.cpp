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
#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "read_isurf.h"
#include "surf.h"
#include "domain.h"
#include "grid.h"
#include "comm.h"
#include "input.h"
#include "my_page.h"
#include "lookup_table.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{NEITHER,BAD,GOOD};
enum{NONE,CHECK,KEEP};
enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};           // several files

#define CHUNK 8192
#define BIG 1.0e20
#define EPSILON 1.0e-16

/* ---------------------------------------------------------------------- */

ReadISurf::ReadISurf(SPARTA *sparta) : Pointers(sparta)
{
  MPI_Comm_rank(world,&me);

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
  if (!grid->exist)
    error->all(FLERR,"Cannot read_isurf before grid is defined");
  if (!surf->implicit)
    error->all(FLERR,"Cannot read_isurf unless global surfs implicit is set");
  if (particle->exist)
    error->all(FLERR,"Cannot read_isurf when particles exist");
  if (domain->axisymmetric)
    error->all(FLERR,"Cannot read_isurf for axisymmetric domains");

  surf->exist = 1;
  dim = domain->dimension;

  if (narg < 6) error->all(FLERR,"Illegal read_isurf command");

  int ggroup = grid->find_group(arg[0]);
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
  count = grid->check_uniform_group(ggroup,nxyz,corner,xyzsize);
  if (nx != nxyz[0] || ny != nxyz[1] || nz != nxyz[2])
    error->all(FLERR,"Read_isurf grid group does not match nx,ny,nz");

  // read grid corner point values
  // create and destroy dictionary of my grid cells in group
  //   used to assign per-grid values to local grid cells

  if (screen && me == 0) fprintf(screen,"Reading isurf file ...\n");

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  create_hash(count,ggroup);

  if (dim == 3) memory->create(cvalues,grid->nlocal,8,"readisurf:cvalues");
  else memory->create(cvalues,grid->nlocal,4,"readisurf:cvalues");

  read_corners(gridfile);

  if (typefile) {
    memory->create(svalues,grid->nlocal,"readisurf:svalues");
    read_types(typefile);
  }

  delete hash;

  // clear out old surfs

  grid->unset_neighbors();
  grid->remove_ghosts();
  grid->clear_surf();

  // create surfs in each grid cell based on corner point values
  // set surf->nsurf and surf->nown
  // if specified, apply group keyword to reset per-surf mask info

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  if (dim == 3) marching_cubes(ggroup);
  else marching_squares(ggroup);

  surf->nown = surf->nlocal;
  bigint nlocal = surf->nlocal;
  MPI_Allreduce(&nlocal,&surf->nsurf,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  
  if (sgrouparg) {
    int sgroup = surf->find_group(arg[sgrouparg]);
    if (sgroup < 0) sgroup = surf->add_group(arg[sgrouparg]);
    int sgroupbit = surf->bitmask[sgroup];
 
    int nslocal = surf->nlocal;
    if (dim == 3) {
      Surf::Tri *tris = surf->tris;
      for (int i = 0; i < nslocal; i++) tris[i].mask |= sgroupbit;
    } else {
      Surf::Line *lines = surf->lines;
      for (int i = 0; i < nslocal; i++) lines[i].mask |= sgroupbit;
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

  // error checks that can be done before surfs are mapped to grid cells

  if (dim == 2) surf->check_watertight_2d(0);
  else surf->check_watertight_3d(0);

  MPI_Barrier(world);
  double time4 = MPI_Wtime();

  // copy surf info into grid cells that own them

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

  // stats

  double time_total = time6-time1;
  double time_s2g = time5-time4;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  CPU time = %g secs\n",time_total);
      fprintf(screen,"  read/marching/check/surf2grid/ghost/inout percent = "
              "%g %g %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total,
              100.0*(time6-time5)/time_total,100.0*(time7-time6)/time_total);
      fprintf(screen,"  surf2grid time = %g secs\n",time_s2g);
      fprintf(screen,"  map/rvous/split percent = %g %g %g\n",
              100.0*grid->tmap/time_s2g,100.0*grid->trvous1/time_s2g,
              100.0*grid->tsplit/time_s2g);
    }
    if (logfile) {
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
      fprintf(logfile,"  read/marching/check/surf2grid/ghost/inout percent = "
              "%g %g %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total,
              100.0*(time6-time5)/time_total,100.0*(time7-time6)/time_total);
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
  int dimension = domain->dimension;

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
    fread(nxyz,sizeof(int),dimension,fp);
  }

  MPI_Bcast(nxyz,dimension,MPI_INT,0,world);

  int flag = 0;
  if (nxyz[0] != nx+1) flag = 1;
  if (nxyz[1] != ny+1) flag = 1;
  if (dimension == 3 && nxyz[2] != nz+1) flag = 1;
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
  int dimension = domain->dimension;

  // proc 0 opens and reads binary file
  // error check the file grid matches input script extent

  if (me == 0) {
    fp = fopen(typefile,"rb");
    if (fp == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open read_isurf type file %s",typefile);
      error->one(FLERR,str);
    }
    fread(nxyz,sizeof(int),dimension,fp);
  }

  MPI_Bcast(nxyz,dimension,MPI_INT,0,world);

  int flag = 0;
  if (nxyz[0] != nx) flag = 1;
  if (nxyz[1] != ny) flag = 1;
  if (dimension == 3 && nxyz[2] != nz) flag = 1;
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
------------------------------------------------------------------------- */

void ReadISurf::create_hash(int count, int igroup)
{
#ifdef SPARTA_MAP
  hash = new std::map<bigint,int>();
#elif defined SPARTA_UNORDERED_MAP
  hash = new std::unordered_map<bigint,int>();
#else
  hash = new std::tr1::unordered_map<bigint,int>();
#endif

  int dim = domain->dimension;
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;
  int groupbit = grid->bitmask[igroup];

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
    } else error->all(FLERR,"Invalid read_isurf command");
  }
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// 2d marching squares algorithm
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   create 2d implicit surfs from grid point values
   treats each grid cell independently
   4 corner points open/solid -> 2^4 = 16 cases
   cases infer 0,1,2 line segments in each grid cell
   order 2 points in each line segment to give normal into flow volume
   treat two saddle point cases (my 9,6) (Wiki 5,10)
     based on ave value at cell center
   follows https://en.wikipedia.org/wiki/Marching_squares
   see 2 sections: Basic algorithm and Disambiguation of saddle points
   treating open circles as flow volume, solid circles as material
   NOTE: Wiki page numbers points counter-clockwise
         SPARTA numbers them in x, then in y
         so bit2 and bit3 are swapped below
         this gives case #s here consistent with Wiki page
------------------------------------------------------------------------- */

void ReadISurf::marching_squares(int igroup)
{
  int i,ipt,isurf,nsurf,which,splitflag;
  int v00,v01,v10,v11,bit0,bit1,bit2,bit3;
  double ave;
  double *lo,*hi;
  surfint *ptr;

  double pt[4][3];
  pt[0][2] = pt[1][2] = pt[2][2] = pt[3][2] = 0.0;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  MyPage<surfint> *csurfs = grid->csurfs;
  int nglocal = grid->nlocal;
  int groupbit = grid->bitmask[igroup];

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    lo = cells[icell].lo;
    hi = cells[icell].hi;

    // cvalues are ordered lower-left, lower-right, upper-left, upper-right
    // Vyx encodes this as 0/1 in each dim

    v00 = cvalues[icell][0];
    v01 = cvalues[icell][1];
    v10 = cvalues[icell][2];
    v11 = cvalues[icell][3];
    
    // make last 2 bits consistent with Wiki page (see NOTE above)

    bit0 = v00 <= thresh ? 0 : 1;
    bit1 = v01 <= thresh ? 0 : 1;
    bit2 = v11 <= thresh ? 0 : 1;
    bit3 = v10 <= thresh ? 0 : 1;
    
    which = (bit3 << 3) + (bit2 << 2) + (bit1 << 1) + bit0;
    splitflag = 0;

    switch (which) {

    case 0: 
      nsurf = 0;
      break;

    case 1:
      nsurf = 1;
      pt[0][0] = lo[0];
      pt[0][1] = interpolate(v00,v10,lo[1],hi[1]);
      pt[1][0] = interpolate(v00,v01,lo[0],hi[0]);
      pt[1][1] = lo[1];
      break;

    case 2:
      nsurf = 1;
      pt[0][0] = interpolate(v00,v01,lo[0],hi[0]);
      pt[0][1] = lo[1];
      pt[1][0] = hi[0];
      pt[1][1] = interpolate(v01,v11,lo[1],hi[1]);
      break;

    case 3:
      nsurf = 1;
      pt[0][0] = lo[0];
      pt[0][1] = interpolate(v00,v10,lo[1],hi[1]);
      pt[1][0] = hi[0];
      pt[1][1] = interpolate(v01,v11,lo[1],hi[1]);
      break;

    case 4:
      nsurf = 1;
      pt[0][0] = hi[0];
      pt[0][1] = interpolate(v01,v11,lo[1],hi[1]);
      pt[1][0] = interpolate(v10,v11,lo[0],hi[0]);
      pt[1][1] = hi[1];
      break;

    case 5:
      nsurf = 2;
      ave = 0.25 * (v00 + v01 + v10 + v11);
      if (ave > thresh) {
        splitflag = 1;
        pt[0][0] = lo[0];
        pt[0][1] = interpolate(v00,v10,lo[1],hi[1]);
        pt[1][0] = interpolate(v10,v11,lo[0],hi[0]);
        pt[1][1] = hi[1];
        pt[2][0] = hi[0];
        pt[2][1] = interpolate(v01,v11,lo[1],hi[1]);
        pt[3][0] = interpolate(v00,v01,lo[0],hi[0]);
        pt[3][1] = lo[1];
      } else {
        pt[0][0] = lo[0];
        pt[0][1] = interpolate(v00,v10,lo[1],hi[1]);
        pt[1][0] = interpolate(v00,v01,lo[0],hi[0]);
        pt[1][1] = lo[1];
        pt[2][0] = hi[0];
        pt[2][1] = interpolate(v01,v11,lo[1],hi[1]);
        pt[3][0] = interpolate(v10,v11,lo[0],hi[0]);
        pt[3][1] = hi[1];
      }
      break;

    case 6:
      nsurf = 1;
      pt[0][0] = interpolate(v00,v01,lo[0],hi[0]);
      pt[0][1] = lo[1];
      pt[1][0] = interpolate(v10,v11,lo[0],hi[0]);
      pt[1][1] = hi[1];
      break;

    case 7:
      nsurf = 1;
      pt[0][0] = lo[0];
      pt[0][1] = interpolate(v00,v10,lo[1],hi[1]);
      pt[1][0] = interpolate(v10,v11,lo[0],hi[0]);
      pt[1][1] = hi[1];
      break;

    case 8:
      nsurf = 1;
      pt[0][0] = interpolate(v10,v11,lo[0],hi[0]);
      pt[0][1] = hi[1];
      pt[1][0] = lo[0];
      pt[1][1] = interpolate(v00,v10,lo[1],hi[1]);
      break;

    case 9:
      nsurf = 1;
      pt[0][0] = interpolate(v10,v11,lo[0],hi[0]);
      pt[0][1] = hi[1];
      pt[1][0] = interpolate(v00,v01,lo[0],hi[0]);
      pt[1][1] = lo[1];
      break;

    case 10:
      nsurf = 2;
      ave = 0.25 * (v00 + v01 + v10 + v11);
      if (ave > thresh) {
        splitflag = 1;
        pt[0][0] = interpolate(v00,v01,lo[0],hi[0]);
        pt[0][1] = lo[1];
        pt[1][0] = lo[0];
        pt[1][1] = interpolate(v00,v10,lo[1],hi[1]);
        pt[2][0] = interpolate(v10,v11,lo[0],hi[0]);
        pt[2][1] = hi[1];
        pt[3][0] = hi[0];
        pt[3][1] = interpolate(v01,v11,lo[1],hi[1]);
      } else {
        pt[0][0] = interpolate(v10,v11,lo[0],hi[0]);
        pt[0][1] = hi[1];
        pt[1][0] = lo[0];
        pt[1][1] = interpolate(v00,v10,lo[1],hi[1]);
        pt[2][0] = interpolate(v00,v01,lo[0],hi[0]);
        pt[2][1] = lo[1];
        pt[3][0] = hi[0];
        pt[3][1] = interpolate(v01,v11,lo[1],hi[1]);
      }
      break;

    case 11:
      nsurf = 1;
      pt[0][0] = interpolate(v10,v11,lo[0],hi[0]);
      pt[0][1] = hi[1];
      pt[1][0] = hi[0];
      pt[1][1] = interpolate(v01,v11,lo[1],hi[1]);
      break;

    case 12:
      nsurf = 1;
      pt[0][0] = hi[0];
      pt[0][1] = interpolate(v01,v11,lo[1],hi[1]);
      pt[1][0] = lo[0];
      pt[1][1] = interpolate(v00,v10,lo[1],hi[1]);
      break;

    case 13: 
      nsurf = 1;
      pt[0][0] = hi[0];
      pt[0][1] = interpolate(v01,v11,lo[1],hi[1]);
      pt[1][0] = interpolate(v00,v01,lo[0],hi[0]);
      pt[1][1] = lo[1];
      break;
    
    case 14:
      nsurf = 1;
      pt[0][0] = interpolate(v00,v01,lo[0],hi[0]);
      pt[0][1] = lo[1];
      pt[1][0] = lo[0];
      pt[1][1] = interpolate(v00,v10,lo[1],hi[1]);
      break;
    
    case 15:
      nsurf = 0;
      break;
    }

    // populate Grid and Surf data structs
    // points will be duplicated, not unique
    // surf ID = cell ID for all surfs in cell
    
    ptr = csurfs->get(nsurf);

    ipt = 0;
    for (i = 0; i < nsurf; i++) {
      if (svalues) surf->add_line(svalues[icell],pt[ipt],pt[ipt+1]);
      else surf->add_line(1,pt[ipt],pt[ipt+1]);
      ipt += 2;
      isurf = surf->nlocal - 1;
      surf->lines[isurf].id = cells[icell].id;
      ptr[i] = isurf;
    }

    cells[icell].nsurf = nsurf;
    if (nsurf) {
      cells[icell].csurfs = ptr;
      cinfo[icell].type = OVERLAP;
    }
  }
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// 3d marching cubes algorithm
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   create 3d implicit surfs from grid point values
------------------------------------------------------------------------- */

void ReadISurf::marching_cubes(int igroup)
{
  int i,j,ipt,isurf,nsurf,icase,which;
  int *ptr;
    
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  MyPage<int> *csurfs = grid->csurfs;
  int nglocal = grid->nlocal;
  int groupbit = grid->bitmask[igroup];
    
  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    lo = cells[icell].lo;
    hi = cells[icell].hi;
        
    // nsurf = # of tris in cell
    // cvalues[8] = 8 corner point values, each is 0 to 255 inclusive
    // thresh = value between 0 and 255 to threshhold on
    // lo[3] = lower left corner pt of grid cell
    // hi[3] = upper right corner pt of grid cell
    // pt = list of 3*nsurf points that are the corner pts of each tri
    
    // cvalues are ordered
    // bottom-lower-left, bottom-lower-right, 
    // bottom-upper-left, bottom-upper-right
    // top-lower-left, top-lower-right, top-upper-left, top-upper-right
    // Vzyx encodes this as 0/1 in each dim
        
    v000 = cvalues[icell][0];
    v001 = cvalues[icell][1];
    v010 = cvalues[icell][2];
    v011 = cvalues[icell][3];
    v100 = cvalues[icell][4];
    v101 = cvalues[icell][5];
    v110 = cvalues[icell][6];
    v111 = cvalues[icell][7];
        
    // make bits 2, 3, 6 and 7 consistent with Lewiner paper (see NOTE above)
        
    bit0 = v000 <= thresh ? 0 : 1;
    bit1 = v001 <= thresh ? 0 : 1;
    bit2 = v011 <= thresh ? 0 : 1;
    bit3 = v010 <= thresh ? 0 : 1;
    bit4 = v100 <= thresh ? 0 : 1;
    bit5 = v101 <= thresh ? 0 : 1;
    bit6 = v111 <= thresh ? 0 : 1;
    bit7 = v110 <= thresh ? 0 : 1;
    
    which = (bit7 << 7) + (bit6 << 6) + (bit5 << 5) + (bit4 << 4) + 
      (bit3 << 3) + (bit2 << 2) + (bit1 << 1) + bit0;
        
    // icase = case of the active cube in [0..15]

    icase = cases[which][0];
    config = cases[which][1];
    subconfig = 0;
    
    //        printf("case %d and config %d lo %1.0f %1.0f %1.0f hi %1.0f %1.0f %1.0f cvalues %d %d %d %d %d %d %d %d\n",icase,config,lo[0],lo[1],lo[2],hi[0],hi[1],hi[2],cvalues[icell][0],cvalues[icell][1],cvalues[icell][2],cvalues[icell][3],cvalues[icell][4],cvalues[icell][5],cvalues[icell][6],cvalues[icell][7]);
    
    switch (icase) {
    case  0:
      nsurf = 0;
      break;
                
    case  1:
      nsurf = add_triangle(tiling1[config], 1);
      break;
                
    case  2:
      nsurf = add_triangle(tiling2[config], 2);
      break;
                
    case  3:
      if (test_face(test3[config]))
        nsurf = add_triangle(tiling3_2[config], 4); // 3.2
      else
        nsurf = add_triangle(tiling3_1[config], 2); // 3.1
      break;
                
    case  4:
      if (test_interior(test4[config],icase))
        nsurf = add_triangle(tiling4_1[config], 2); // 4.1.1
      else
        nsurf = add_triangle(tiling4_2[config], 6); // 4.1.2
      break;
                
    case  5:
      nsurf = add_triangle(tiling5[config], 3);
      break;
                
    case  6:
      if (test_face(test6[config][0]))
        nsurf = add_triangle(tiling6_2[config], 5); // 6.2
      else {
        if (test_interior(test6[config][1],icase))
          nsurf = add_triangle(tiling6_1_1[config], 3); // 6.1.1
        else {
          nsurf = add_triangle(tiling6_1_2[config], 9); // 6.1.2
        }
      }
      break;
                
    case  7:
      if (test_face(test7[config][0])) subconfig +=  1;
      if (test_face(test7[config][1])) subconfig +=  2;
      if (test_face(test7[config][2])) subconfig +=  4;
      switch (subconfig) {
      case 0:
        nsurf = add_triangle(tiling7_1[config], 3); break;
      case 1:
        nsurf = add_triangle(tiling7_2[config][0], 5); break;
      case 2:
        nsurf = add_triangle(tiling7_2[config][1], 5); break;
      case 3:
        nsurf = add_triangle(tiling7_3[config][0], 9); break;
      case 4:
        nsurf = add_triangle(tiling7_2[config][2], 5); break;
      case 5:
        nsurf = add_triangle(tiling7_3[config][1], 9); break;
      case 6:
        nsurf = add_triangle(tiling7_3[config][2], 9); break;
      case 7:
        if (test_interior(test7[config][3],icase))
          nsurf = add_triangle(tiling7_4_2[config], 9);
        else
          nsurf = add_triangle(tiling7_4_1[config], 5);
        break;
      };
      break;
                
    case  8:
      nsurf = add_triangle(tiling8[config], 2);
      break;
                
    case  9:
      nsurf = add_triangle(tiling9[config], 4);
      break;
                
    case 10:
      if (test_face(test10[config][0])) {
        if (test_face(test10[config][1]))
          nsurf = add_triangle(tiling10_1_1_[config], 4); // 10.1.1
        else {
          nsurf = add_triangle(tiling10_2[config], 8); // 10.2
        }
      } else {
        if (test_face(test10[config][1])) {
          nsurf = add_triangle(tiling10_2_[config], 8); // 10.2
        } else {
          if (test_interior(test10[config][2],icase))
            nsurf = add_triangle(tiling10_1_1[config], 4); // 10.1.1
          else
            nsurf = add_triangle(tiling10_1_2[config], 8); // 10.1.2
        }
      }
      break;
                
    case 11:
      nsurf = add_triangle(tiling11[config], 4);
      break;
                
    case 12:
      if (test_face(test12[config][0])) {
        if (test_face(test12[config][1]))
          nsurf = add_triangle(tiling12_1_1_[config], 4); // 12.1.1
        else {
          nsurf = add_triangle(tiling12_2[config], 8); // 12.2
        }
      } else {
        if (test_face(test12[config][1])) {
          nsurf = add_triangle(tiling12_2_[config], 8); // 12.2
        } else {
          if (test_interior(test12[config][2],icase))
            nsurf = add_triangle(tiling12_1_1[config], 4); // 12.1.1
          else
            nsurf = add_triangle(tiling12_1_2[config], 8); // 12.1.2
        }
      }
      break;
                
    case 13:
      if (test_face(test13[config][0])) subconfig +=  1;
      if (test_face(test13[config][1])) subconfig +=  2;
      if (test_face(test13[config][2])) subconfig +=  4;
      if (test_face(test13[config][3])) subconfig +=  8;
      if (test_face(test13[config][4])) subconfig += 16;
      if (test_face(test13[config][5])) subconfig += 32;
                               
      switch (subconfig13[subconfig]) {
      case 0:/* 13.1 */
        nsurf = add_triangle(tiling13_1[config], 4); break;

      case 1:/* 13.2 */
        nsurf = add_triangle(tiling13_2[config][0], 6); break;
      case 2:/* 13.2 */
        nsurf = add_triangle(tiling13_2[config][1], 6); break;
      case 3:/* 13.2 */
        nsurf = add_triangle(tiling13_2[config][2], 6); break;
      case 4:/* 13.2 */
        nsurf = add_triangle(tiling13_2[config][3], 6); break;
      case 5:/* 13.2 */
        nsurf = add_triangle(tiling13_2[config][4], 6); break;
      case 6:/* 13.2 */
        nsurf = add_triangle(tiling13_2[config][5], 6); break;

      case 7:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][0], 10); break;
      case 8:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][1], 10); break;
      case 9:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][2], 10); break;
      case 10:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][3], 10); break;
      case 11:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][4], 10); break;
      case 12:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][5], 10); break;
      case 13:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][6], 10); break;
      case 14:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][7], 10); break;
      case 15:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][8], 10); break;
      case 16:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][9], 10); break;
      case 17:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][10], 10); break;
      case 18:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][11], 10); break;
        
      case 19:/* 13.4 */
        nsurf = add_triangle(tiling13_4[config][0], 12); break;
      case 20:/* 13.4 */
        nsurf = add_triangle(tiling13_4[config][1], 12); break;
      case 21:/* 13.4 */
        nsurf = add_triangle(tiling13_4[config][2], 12); break;
      case 22:/* 13.4 */
        nsurf = add_triangle(tiling13_4[config][3], 12); break;

      case 23:/* 13.5 */
        subconfig = 0;
        if (test_interior(test13[config][6],icase))
          nsurf = add_triangle(tiling13_5_1[config][0], 6);
        else
          nsurf = add_triangle(tiling13_5_2[config][0], 10);
        break;

      case 24:/* 13.5 */
        subconfig = 1;
        if (test_interior(test13[config][6],icase))
          nsurf = add_triangle(tiling13_5_1[config][1], 6);
        else
          nsurf = add_triangle(tiling13_5_2[config][1], 10);
        break;

      case 25:/* 13.5 */
        subconfig = 2;
        if (test_interior(test13[config][6],icase))
          nsurf = add_triangle(tiling13_5_1[config][2], 6);
        else
          nsurf = add_triangle(tiling13_5_2[config][2], 10);
        break;

      case 26:/* 13.5 */
        subconfig = 3;
        if (test_interior(test13[config][6],icase))
          nsurf = add_triangle(tiling13_5_1[config][3], 6);
        else
          nsurf = add_triangle(tiling13_5_2[config][3], 10);
        break;

      case 27:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][0], 10); break;
      case 28:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][1], 10); break;
      case 29:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][2], 10); break;
      case 30:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][3], 10); break;
      case 31:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][4], 10); break;
      case 32:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][5], 10); break;
      case 33:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][6], 10); break;
      case 34:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][7], 10); break;
      case 35:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][8], 10); break;
      case 36:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][9], 10); break;
      case 37:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][10], 10); break;
      case 38:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][11], 10); break;
        
      case 39:/* 13.2 */
        nsurf = add_triangle(tiling13_2_[config][0], 6); break;
      case 40:/* 13.2 */
        nsurf = add_triangle(tiling13_2_[config][1], 6); break;
      case 41:/* 13.2 */
        nsurf = add_triangle(tiling13_2_[config][2], 6); break;
      case 42:/* 13.2 */
        nsurf = add_triangle(tiling13_2_[config][3], 6); break;
      case 43:/* 13.2 */
        nsurf = add_triangle(tiling13_2_[config][4], 6); break;
      case 44:/* 13.2 */
        nsurf = add_triangle(tiling13_2_[config][5], 6); break;
        
      case 45:/* 13.1 */
        nsurf = add_triangle(tiling13_1_[config], 4); break;
        
      default:
        printf("Marching Cubes: Impossible case 13?\n");  print_cube();
      }
      break;
                
    case 14:
      nsurf = add_triangle(tiling14[config], 4);
      break;
    };
        
    // populate Grid and Surf data structs
    // points will be duplicated, not unique
    // surf ID = cell ID for all surfs in cell
        
    ptr = csurfs->get(nsurf);
        
    ipt = 0;
    for (i = 0; i < nsurf; i++) {
      if (svalues) surf->add_tri(svalues[icell],pt[ipt],pt[ipt+1],pt[ipt+2]);
      else surf->add_tri(1,pt[ipt],pt[ipt+1],pt[ipt+2]);
      ipt += 3;
      isurf = surf->nlocal - 1;
      surf->tris[isurf].id = cells[icell].id;
      ptr[i] = isurf;
    }
        
    cells[icell].nsurf = nsurf;
    if (nsurf) {
      cells[icell].csurfs = ptr;
      cinfo[icell].type = OVERLAP;
    }
  }
}

/* ----------------------------------------------------------------------
   adding triangles
------------------------------------------------------------------------- */

int ReadISurf::add_triangle(int *trig, int n)
{
  for(int t = 0; t < 3*n; t++) {
    switch (trig[t]) {
    case 0:
      pt[t][0] = interpolate(v000,v001,lo[0],hi[0]);
      pt[t][1] = lo[1];
      pt[t][2] = lo[2];
      break;
    case 1:
      pt[t][0] = hi[0];
      pt[t][1] = interpolate(v001,v011,lo[1],hi[1]);
      pt[t][2] = lo[2];
      break;
    case 2:
      pt[t][0] = interpolate(v010,v011,lo[0],hi[0]);
      pt[t][1] = hi[1];
      pt[t][2] = lo[2];
      break;
    case 3:
      pt[t][0] = lo[0];
      pt[t][1] = interpolate(v000,v010,lo[1],hi[1]);
      pt[t][2] = lo[2];
      break;
    case 4:
      pt[t][0] = interpolate(v100,v101,lo[0],hi[0]);
      pt[t][1] = lo[1];
      pt[t][2] = hi[2];
      break;
    case 5:
      pt[t][0] = hi[0];
      pt[t][1] = interpolate(v101,v111,lo[1],hi[1]);
      pt[t][2] = hi[2];
      break;
    case 6:
      pt[t][0] = interpolate(v110,v111,lo[0],hi[0]);
      pt[t][1] = hi[1];
      pt[t][2] = hi[2];
      break;
    case 7:
      pt[t][0] = lo[0];
      pt[t][1] = interpolate(v100,v110,lo[1],hi[1]);
      pt[t][2] = hi[2];
      break;
    case 8:
      pt[t][0] = lo[0];
      pt[t][1] = lo[1];
      pt[t][2] = interpolate(v000,v100,lo[2],hi[2]);
      break;
    case 9:
      pt[t][0] = hi[0];
      pt[t][1] = lo[1];
      pt[t][2] = interpolate(v001,v101,lo[2],hi[2]);
      break;
    case 10:
      pt[t][0] = hi[0];
      pt[t][1] = hi[1];
      pt[t][2] = interpolate(v011,v111,lo[2],hi[2]);
      break;
    case 11:
      pt[t][0] = lo[0];
      pt[t][1] = hi[1];
      pt[t][2] = interpolate(v010,v110,lo[2],hi[2]);
      break;
    case 12: {
      int u = 0;
      pt[t][0] = pt[t][1] = pt[t][2] = 0.0;
      if (bit0 ^ bit1) {
        ++u;
        pt[t][0] += interpolate(v000,v001,lo[0],hi[0]);
        pt[t][1] += lo[1];
        pt[t][2] += lo[2];
      }
      if (bit1 ^ bit2) {
        ++u;
        pt[t][0] += hi[0];
        pt[t][1] += interpolate(v001,v011,lo[1],hi[1]);
        pt[t][2] += lo[2];
      }
      if (bit2 ^ bit3) {
        ++u;
        pt[t][0] += interpolate(v010,v011,lo[0],hi[0]);
        pt[t][1] += hi[1];
        pt[t][2] += lo[2];
      }
      if (bit3 ^ bit0) {
        ++u;
        pt[t][0] += lo[0];
        pt[t][1] += interpolate(v000,v010,lo[1],hi[1]);
        pt[t][2] += lo[2];
      }
      if (bit4 ^ bit5) {
        ++u;
        pt[t][0] += interpolate(v100,v101,lo[0],hi[0]);
        pt[t][1] += lo[1];
        pt[t][2] += hi[2];
      }
      if (bit5 ^ bit6) {
        ++u;
        pt[t][0] += hi[0];
        pt[t][1] += interpolate(v101,v111,lo[1],hi[1]);
        pt[t][2] += hi[2];
      }
      if (bit6 ^ bit7) {
        ++u;
        pt[t][0] += interpolate(v110,v111,lo[0],hi[0]);
        pt[t][1] += hi[1];
        pt[t][2] += hi[2];
      }
      if (bit7 ^ bit4) {
        ++u;
        pt[t][0] += lo[0];
        pt[t][1] += interpolate(v100,v110,lo[1],hi[1]);
        pt[t][2] += hi[2];
      }
      if (bit0 ^ bit4) {
        ++u;
        pt[t][0] += lo[0];
        pt[t][1] += lo[1];
        pt[t][2] += interpolate(v000,v100,lo[2],hi[2]);
      }
      if (bit1 ^ bit5) {
        ++u;
        pt[t][0] += hi[0];
        pt[t][1] += lo[1];
        pt[t][2] += interpolate(v001,v101,lo[2],hi[2]);
      }
      if (bit2 ^ bit6) {
        ++u;
        pt[t][0] += hi[0];
        pt[t][1] += hi[1];
        pt[t][2] += interpolate(v011,v111,lo[2],hi[2]);
      }
      if (bit3 ^ bit7) {
        ++u;
        pt[t][0] += lo[0];
        pt[t][1] += hi[1];
        pt[t][2] += interpolate(v010,v110,lo[2],hi[2]);
      }

      pt[t][0] /= static_cast<double> (u);
      pt[t][1] /= static_cast<double> (u);
      pt[t][2] /= static_cast<double> (u);
      break;
    }

    default:
      break;
    }

    //        printf("pt number %d n %d case %d and config %d lo %1.0f %1.0f %1.0f hi %1.0f %1.0f %1.0f cvalues %d %d %d %d %d %d %d %d coords %f %f %f\n" ,t,n,icase,config,lo[0],lo[1],lo[2],hi[0],hi[1],hi[2],v000,v001,v011,v010,v100,v101,v111,v110,pt[t][0],pt[t][1],pt[t][2]);

  }

  return n;
}

/* ----------------------------------------------------------------------
   test a face
   if face > 0 return true if the face contains a part of the surface
------------------------------------------------------------------------- */

bool ReadISurf::test_face(int face)
{
  double A,B,C,D;
    
  switch (face) {
  case -1: 
  case 1:  
    A = static_cast<double>(v000);  B = static_cast<double>(v100);
    C = static_cast<double>(v101);  D = static_cast<double>(v001);  
    break;

  case -2: 
  case 2:  
    A = static_cast<double>(v001);  B = static_cast<double>(v101);
    C = static_cast<double>(v111);  D = static_cast<double>(v011);
    break;

  case -3:
  case 3:  
    A = static_cast<double>(v011);  B = static_cast<double>(v111);
    C = static_cast<double>(v110);  D = static_cast<double>(v010);
    break;

  case -4:
  case 4:  
    A = static_cast<double>(v010);  B = static_cast<double>(v110);
    C = static_cast<double>(v100);  D = static_cast<double>(v000);
    break;

  case -5:
  case 5:
    A = static_cast<double>(v000);  B = static_cast<double>(v010);
    C = static_cast<double>(v011);  D = static_cast<double>(v001);
    break;

  case -6:
  case 6:
    A = static_cast<double>(v100);  B = static_cast<double>(v110);
    C = static_cast<double>(v111);  D = static_cast<double>(v101);
    break;

  default: 
    printf("Invalid face code %d\n", face);
    print_cube(); A = B = C = D = 0.0;
  };
    
  if (fabs(A*C - B*D) < EPSILON) return face >= 0;
  return face * A * (A*C - B*D) >= 0 ;  // face and A invert signs
}

/* ----------------------------------------------------------------------
   test the interior of a cube
   icase = case of the active cube in [0..15]
   if s ==  7, return true if the interior is empty
   if s == -7, return false if the interior is empty
------------------------------------------------------------------------- */

bool ReadISurf::test_interior(int s, int icase)
{
  double t,a,b,At=0.0,Bt=0.0,Ct=0.0,Dt=0.0;
  int test = 0;
  int edge = -1;   // reference edge of the triangulation
    
  switch (icase) {
  case 4:
  case 10:
    a = static_cast<double>(v100 - v000) * static_cast<double>(v111 - v011)
      - static_cast<double>(v110 - v010) * static_cast<double>(v101 - v001);
    b =  static_cast<double>(v011) * static_cast<double>(v100 - v000)
      + static_cast<double>(v000) * static_cast<double>(v111 - v011)
      - static_cast<double>(v001) * static_cast<double>(v110 - v010)
      - static_cast<double>(v010) * static_cast<double>(v101 - v001);
    t = - b / (2*a);
    if (t<0.0 || t>1.0) return s>0;
            
    At = static_cast<double>(v000) + static_cast<double>(v100 - v000) * t;
    Bt = static_cast<double>(v010) + static_cast<double>(v110 - v010) * t;
    Ct = static_cast<double>(v011) + static_cast<double>(v111 - v011) * t;
    Dt = static_cast<double>(v001) + static_cast<double>(v101 - v001) * t;
    break;
            
    case 6:
    case 7:
    case 12:
    case 13:
      switch (icase) {
      case  6: edge = test6 [config][2]; break;
      case  7: edge = test7 [config][4]; break;
      case 12: edge = test12[config][3]; break;
      case 13: edge = tiling13_5_1[config][subconfig][0]; break;
      }
      switch (edge) {
      case 0:
        t  = static_cast<double>(v000) / static_cast<double>(v000 - v001);
        At = 0.0;
        Bt = static_cast<double>(v010) + static_cast<double>(v011 - v010) * t;
        Ct = static_cast<double>(v110) + static_cast<double>(v111 - v110) * t;
        Dt = static_cast<double>(v100) + static_cast<double>(v101 - v100) * t;
        break;
      case 1:
        t  = static_cast<double>(v001) / static_cast<double>(v001 - v011);
        At = 0.0;
        Bt = static_cast<double>(v000) + static_cast<double>(v010 - v000) * t;
        Ct = static_cast<double>(v100) + static_cast<double>(v110 - v100) * t;
        Dt = static_cast<double>(v101) + static_cast<double>(v111 - v101) * t;
        break;
      case 2:
        t  = static_cast<double>(v011) / static_cast<double>(v011 - v010);
        At = 0.0;
        Bt = static_cast<double>(v001) + static_cast<double>(v000 - v001) * t;
        Ct = static_cast<double>(v101) + static_cast<double>(v100 - v101) * t;
        Dt = static_cast<double>(v111) + static_cast<double>(v110 - v111) * t;
        break;
      case 3:
        t  = static_cast<double>(v010) / static_cast<double>(v010 - v000);
        At = 0.0;
        Bt = static_cast<double>(v011) + static_cast<double>(v001 - v011) * t;
        Ct = static_cast<double>(v111) + static_cast<double>(v101 - v111) * t;
        Dt = static_cast<double>(v110) + static_cast<double>(v100 - v110) * t;
        break;
      case 4:
        t  = static_cast<double>(v100) / static_cast<double>(v100 - v101);
        At = 0.0;
        Bt = static_cast<double>(v110) + static_cast<double>(v111 - v110) * t;
        Ct = static_cast<double>(v010) + static_cast<double>(v011 - v010) * t;
        Dt = static_cast<double>(v000) + static_cast<double>(v001 - v000) * t;
        break;
      case 5:
        t  = static_cast<double>(v101) / static_cast<double>(v101 - v111);
        At = 0.0;
        Bt = static_cast<double>(v100) + static_cast<double>(v110 - v100) * t;
        Ct = static_cast<double>(v000) + static_cast<double>(v010 - v000) * t;
        Dt = static_cast<double>(v001) + static_cast<double>(v011 - v001) * t;
        break;
      case 6:
        t  = static_cast<double>(v111) / static_cast<double>(v111 - v110);
        At = 0.0;
        Bt = static_cast<double>(v101) + static_cast<double>(v100 - v101) * t;
        Ct = static_cast<double>(v001) + static_cast<double>(v000 - v001) * t;
        Dt = static_cast<double>(v011) + static_cast<double>(v010 - v011) * t;
        break;
      case 7:
        t  = static_cast<double>(v110) / static_cast<double>(v110 - v100);
        At = 0.0;
        Bt = static_cast<double>(v111) + static_cast<double>(v101 - v111) * t;
        Ct = static_cast<double>(v011) + static_cast<double>(v001 - v011) * t;
        Dt = static_cast<double>(v010) + static_cast<double>(v000 - v010) * t;
        break;
      case 8:
        t  = static_cast<double>(v000) / static_cast<double>(v000 - v100);
        At = 0.0;
        Bt = static_cast<double>(v010) + static_cast<double>(v110 - v010) * t;
        Ct = static_cast<double>(v011) + static_cast<double>(v111 - v011) * t;
        Dt = static_cast<double>(v001) + static_cast<double>(v101 - v001) * t;
        break;
      case 9:
        t  = static_cast<double>(v001) / static_cast<double>(v001 - v101);
        At = 0.0;
        Bt = static_cast<double>(v000) + static_cast<double>(v100 - v000) * t;
        Ct = static_cast<double>(v010) + static_cast<double>(v110 - v010) * t;
        Dt = static_cast<double>(v011) + static_cast<double>(v111 - v011) * t;
        break;
      case 10:
        t  = static_cast<double>(v011) / static_cast<double>(v011 - v111);
        At = 0.0;
        Bt = static_cast<double>(v001) + static_cast<double>(v101 - v001) * t;
        Ct = static_cast<double>(v000) + static_cast<double>(v100 - v000) * t;
        Dt = static_cast<double>(v010) + static_cast<double>(v110 - v010) * t;
        break;
      case 11:
        t  = static_cast<double>(v010) / static_cast<double>(v010 - v110);
        At = 0.0;
        Bt = static_cast<double>(v011) + static_cast<double>(v111 - v011) * t;
        Ct = static_cast<double>(v001) + static_cast<double>(v101 - v001) * t;
        Dt = static_cast<double>(v000) + static_cast<double>(v100 - v000) * t;
        break;

      default: 
        printf("Invalid edge %d\n", edge);  print_cube();  break;
      }
      break;
      
  default: 
    printf("Invalid ambiguous case %d\n", icase);
    print_cube();
    break;
  }
    
  if (At >= 0.0) test ++;
  if (Bt >= 0.0) test += 2;
  if (Ct >= 0.0) test += 4;
  if (Dt >= 0.0) test += 8;
  switch (test) {
  case  0: return s>0;
  case  1: return s>0;
  case  2: return s>0;
  case  3: return s>0;
  case  4: return s>0;
  case  5: 
    if (At * Ct - Bt * Dt <  EPSILON) return s>0;
    break;
  case  6: return s>0;
  case  7: return s<0;
  case  8: return s>0;
  case  9: return s>0;
  case 10: 
    if (At * Ct - Bt * Dt >= EPSILON) return s>0;
    break;
  case 11: return s<0;
  case 12: return s>0;
  case 13: return s<0;
  case 14: return s<0;
  case 15: return s<0;
  }
  
  return s<0;
}

/* ----------------------------------------------------------------------
   print cube for debugging
------------------------------------------------------------------------- */

void ReadISurf::print_cube() { 
  printf("\t %d %d %d %d %d %d %d %d\n",
          v000, v001, v011, v010, v100, v101, v111, v110);
}

/* ----------------------------------------------------------------------
   interpolate for marching squares or cubes
   lo/hi = coordinates of end points of edge of square
   v0/v1 = values at lo/hi end points
   value = interpolated coordinate for thresh value
------------------------------------------------------------------------- */

double ReadISurf::interpolate(int v0, int v1, double lo, double hi)
{
  double value = lo + (hi-lo)*(thresh-v0)/(v1-v0);
  value = MAX(value,lo);
  value = MIN(value,hi);
  return value;
}
