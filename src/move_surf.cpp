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

#include "stdio.h"
#include "string.h"
#include "move_surf.h"
#include "surf.h"
#include "grid.h"
#include "comm.h"
#include "update.h"
#include "domain.h"
#include "input.h"
#include "math_extra.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathExtra;
using namespace MathConst;

enum{READFILE,TRANSLATE,ROTATE};
enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};           // several files

#define MAXLINE 256

/* ---------------------------------------------------------------------- */

MoveSurf::MoveSurf(SPARTA *sparta) : Pointers(sparta)
{
  me = comm->me;
  nprocs = comm->nprocs;

  // pselect = 1 if point is moved, else 0

  if (domain->dimension == 2)
    memory->create(pselect,2*surf->nsurf,"move_surf:pselect");
  else
    memory->create(pselect,3*surf->nsurf,"move_surf:pselect");

  file = NULL;
  fp = NULL;
}

/* ---------------------------------------------------------------------- */

MoveSurf::~MoveSurf()
{
  memory->destroy(pselect);
  delete [] file;
  if (fp) fclose(fp);
}

/* ---------------------------------------------------------------------- */

void MoveSurf::command(int narg, char **arg)
{
  if (!surf->exist)
    error->all(FLERR,"Cannot move_surf with no surf elements defined");

  if (surf->distributed)
    error->all(FLERR,
               "Cannot yet use move_surf with distributed surf elements");

  if (narg < 2) error->all(FLERR,"Illegal move_surf command");

  // process command-line args

  int igroup = surf->find_group(arg[0]);
  if (igroup < 0) error->all(FLERR,"Move_surf group ID does not exist");
  groupbit = surf->bitmask[igroup];

  process_args(narg-1,&arg[1]);
  mode = 0;

  if (action == READFILE)
    error->all(FLERR,"Move_surf file option is not yet implemented");

  // perform surface move

  if (me == 0) {
    if (screen) fprintf(screen,"Moving surfs ...\n");
    if (logfile) fprintf(logfile,"Moving surfs ...\n");
  }

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  dim = domain->dimension;

  // sort particles

  if (particle->exist) particle->sort();

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  // move line/tri points via chosen action by full amount

  if (dim == 2) move_lines(1.0,surf->lines);
  else move_tris(1.0,surf->tris);

  // remake list of surf elements I own
  // assign split cell particles to parent split cell
  // assign surfs to grid cells

  grid->unset_neighbors();
  grid->remove_ghosts();

  if (particle->exist && grid->nsplitlocal) {
    Grid::ChildCell *cells = grid->cells;
    int nglocal = grid->nlocal;
    for (int icell = 0; icell < nglocal; icell++)
      if (cells[icell].nsplit > 1)
        grid->combine_split_cell_particles(icell,1);
  }

  grid->clear_surf();
  grid->surf2grid(1);

  if (dim == 2) surf->check_point_near_surf_2d();
  else surf->check_point_near_surf_3d();

  if (dim == 2) surf->check_watertight_2d();
  else surf->check_watertight_3d();

  MPI_Barrier(world);
  double time3 = MPI_Wtime();

  // re-setup owned and ghost cell info

  grid->setup_owned();
  grid->acquire_ghosts();
  grid->reset_neighbors();
  comm->reset_neighbors();

  MPI_Barrier(world);
  double time4 = MPI_Wtime();

  // flag cells and corners as OUTSIDE or INSIDE
  // reallocate per grid cell arrays in per grid computes
  //   local grid cell counts could have changed due to split cell changes

  grid->set_inout();
  grid->type_check();

  // DEBUG
  //grid->debug();

  MPI_Barrier(world);
  double time5 = MPI_Wtime();

  // remove particles as needed due to surface move

  bigint ndeleted;
  if (particle->exist) ndeleted = remove_particles();

  MPI_Barrier(world);
  double time6 = MPI_Wtime();

  double time_total = time6-time1;

  if (comm->me == 0) {
    if (screen) {
      if (particle->exist)
        fprintf(screen,"  " BIGINT_FORMAT " deleted particles\n",ndeleted);
      fprintf(screen,"  CPU time = %g secs\n",time_total);
      fprintf(screen,"  sort/surf2grid/ghost/inout/particle percent = "
              "%g %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total,
              100.0*(time6-time5)/time_total);
    }
    if (logfile) {
      if (particle->exist)
        fprintf(logfile,"  " BIGINT_FORMAT " deleted particles\n",ndeleted);
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
      fprintf(logfile,"  sort/surf2grid/ghost/inout/particle percent = "
              "%g %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total,
              100.0*(time6-time5)/time_total);
    }
  }
}

/* ----------------------------------------------------------------------
   process command args for both move_surf and fix move/surf
------------------------------------------------------------------------- */

void MoveSurf::process_args(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal move surf command");

  int iarg = 0;
  if (strcmp(arg[0],"file") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal move surf command");
    action = READFILE;
    int n = strlen(arg[1]) + 1;
    file = new char[n];
    strcpy(file,arg[1]);
    n = strlen(arg[2]) + 1;
    entry = new char[n];
    strcpy(entry,arg[2]);
    iarg = 3;
  } else if (strcmp(arg[0],"trans") == 0) {
    if (narg < 4) error->all(FLERR,"Illegal move surf command");
    action = TRANSLATE;
    delta[0] = input->numeric(FLERR,arg[1]);
    delta[1] = input->numeric(FLERR,arg[2]);
    delta[2] = input->numeric(FLERR,arg[3]);
    if (domain->dimension == 2 && delta[2] != 0.0)
      error->all(FLERR,"Invalid move surf translation for 2d simulation");
    iarg = 4;
  } else if (strcmp(arg[0],"rotate") == 0) {
    if (narg < 8) error->all(FLERR,"Illegal move surf command");
    action = ROTATE;
    theta = input->numeric(FLERR,arg[1]);
    rvec[0] = input->numeric(FLERR,arg[2]);
    rvec[1] = input->numeric(FLERR,arg[3]);
    rvec[2] = input->numeric(FLERR,arg[4]);
    origin[0] = input->numeric(FLERR,arg[5]);
    origin[1] = input->numeric(FLERR,arg[6]);
    origin[2] = input->numeric(FLERR,arg[7]);
    if (domain->dimension == 2 && (rvec[0] != 0.0 || rvec[1] != 0.0))
      error->all(FLERR,"Invalid move surf rotation for 2d simulation");
    if (rvec[0] == 0.0 && rvec[1] == 0.0 && rvec[2] == 0.0)
      error->all(FLERR,"Invalid move surf rotation");
    theta *= MY_PI/180.0;
    MathExtra::norm3(rvec);
    iarg = 8;
  } else error->all(FLERR,"Illegal move surf command");

  // optional args

  connectflag = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"connect") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal move surf command");
      if (strcmp(arg[iarg+1],"yes") == 0) connectflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) connectflag = 0;
      iarg += 2;
    } else error->all(FLERR,"Illegal move surf command");
  }
}

/* ----------------------------------------------------------------------
   move points in lines via specified action
   each method sets pselect = 1 for moved points
   fraction = portion of full distance points should move
------------------------------------------------------------------------- */

void MoveSurf::move_lines(double fraction, Surf::Line *origlines)
{
  if (connectflag && groupbit != 1) connect_2d_pre();

  if (action == READFILE) {
    readfile();
    update_points(fraction);
  }
  else if (action == TRANSLATE) translate_2d(fraction,origlines);
  else if (action == ROTATE) rotate_2d(fraction,origlines);

  if (connectflag && groupbit != 1) connect_2d_post();

  surf->compute_line_normal(0);

  // check that all points are still inside simulation box

  surf->check_point_inside(0);
}

/* ----------------------------------------------------------------------
   move points in triangles via specified action
   each method sets pselect = 1 for moved points
   fraction = portion of full distance points should move
------------------------------------------------------------------------- */

void MoveSurf::move_tris(double fraction, Surf::Tri *origtris)
{
  if (connectflag && groupbit != 1) connect_3d_pre();

  if (action == READFILE) {
    readfile();
    update_points(fraction);
  }
  else if (action == TRANSLATE) translate_3d(fraction,origtris);
  else if (action == ROTATE) rotate_3d(fraction,origtris);

  if (connectflag && groupbit != 1) connect_3d_post();

  surf->compute_tri_normal(0);

  // check that all points are still inside simulation box

  surf->check_point_inside(0);
}

/* ----------------------------------------------------------------------
   read entry of new point coords from file
------------------------------------------------------------------------- */

void MoveSurf::readfile()
{
  /*
  int i;
  char line[MAXLINE];
  char *word,*eof;

  // open point file if necessary
  // if already open, will just continue scanning below
  // NOTE: allow for file name with wildcard char

  if (me == 0 && fp == NULL) {
    fp = fopen(file,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open move surf file %s",file);
      error->one(FLERR,str);
    }
  }

  // loop until section found that matches entry

  if (me == 0) {
    while (1) {
      if (fgets(line,MAXLINE,fp) == NULL)
        error->one(FLERR,"Did not find entry in move surf file");
      if (strspn(line," \t\n\r") == strlen(line)) continue;  // blank line
      if (line[0] == '#') continue;                          // comment
      word = strtok(line," \t\n\r");
      if (strcmp(word,entry) != 0) continue;          // non-matching entry
      if (fgets(line,MAXLINE,fp) == NULL)
        error->one(FLERR,"Incompete entry in move surf file");
      word = strtok(line," \t\n\r");        // npoints value after entry
      nread = input->inumeric(FLERR,word);
    }
  }

  // allocate index and coord arrays for nread points

  MPI_Bcast(&nread,1,MPI_INT,0,world);

  if (oldcoord) {
    memory->destroy(readindex);
    memory->destroy(oldcoord);
    memory->destroy(newcoord);
    memory->create(readindex,nread,"move_surf:readindex");
    memory->create(oldcoord,nread,3,"move_surf:oldcoord");
    memory->create(newcoord,nread,3,"move_surf:newcoord");
  }

  // read nread point coords in entry
  // store old current point and new point coords so can move by fraction
  // rindex = ID (index) of this read-in point in master list
  // skip points that are out-of-range

  Surf::Point *pts = surf->pts;
  int npoint = surf->npoint;

  if (me == 0) {
    int id;
    double x,y,z;

    for (int i = 0; i < nread; i++) {
      eof = fgets(line,MAXLINE,fp);
      if (eof == NULL) error->one(FLERR,"Incomplete entry in move surf file");
      id = input->inumeric(FLERR,strtok(line," \t\n\r"));
      x = input->numeric(FLERR,strtok(NULL," \t\n\r"));
      y = input->numeric(FLERR,strtok(NULL," \t\n\r"));
      if (dim == 3) z = input->numeric(FLERR,strtok(NULL," \t\n\r"));
      else z = 0.0;
      if (id < 1 || id > npoint)
        error->one(FLERR,"Invalid point index in move surf file");
      id--;
      readindex[i] = id;
      oldcoord[i][0] = pts[id].x[0];
      oldcoord[i][1] = pts[id].x[1];
      oldcoord[i][2] = pts[id].x[2];
      newcoord[i][0] = x;
      newcoord[i][1] = y;
      newcoord[i][2] = z;
    }
  }

  // broadcast point info to all procs

  MPI_Bcast(readindex,nread,MPI_INT,0,world);
  MPI_Bcast(&oldcoord[0][0],3*nread,MPI_DOUBLE,0,world);
  MPI_Bcast(&newcoord[0][0],3*nread,MPI_DOUBLE,0,world);

  // pselect[I] = index of Ith surf point in nread points (for now)
  // NOTE: check that same surf point does not appear twice in nread list?

  for (i = 0; i < npoint; i++) pselect[i] = -1;
  for (i = 0; i < nread; i++) pselect[readindex[i]] = i;

  int *rflag;
  memory->create(rflag,nread,"move_surf:rflag");
  for (i =0; i < nread; i++) rflag[i] = 0;

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  int nline = surf->nline;
  int ntri = surf->ntri;
  int p1,p2,p3;

  if (dim == 2) {
    for (i = 0; i < nline; i++) {
      if (!(lines[i].mask & groupbit)) continue;
      p1 = lines[i].p1;
      p2 = lines[i].p2;
      if (pselect[p1] >= 0) rflag[pselect[p1]] = 1;
      if (pselect[p2] >= 0) rflag[pselect[p2]] = 1;
    }
  } else {
    for (i = 0; i < ntri; i++) {
      if (!(tris[i].mask & groupbit)) continue;
      p1 = tris[i].p1;
      p2 = tris[i].p2;
      p3 = tris[i].p3;
      if (pselect[p1] >= 0) rflag[pselect[p1]] = 1;
      if (pselect[p2] >= 0) rflag[pselect[p2]] = 1;
      if (pselect[p3] >= 0) rflag[pselect[p3]] = 1;
    }
  }


  // pselect[I] = 1 if Ith surf point is moved by nread points, else 0

  for (i = 0; i < npoint; i++) pselect[i] = 0;
  for (i = 0; i < nread; i++)
    if (rflag[i]) pselect[readindex[i]] = 1;

  // clean up

  memory->destroy(rflag);
  */
}

/* ----------------------------------------------------------------------
   update points using info from file
------------------------------------------------------------------------- */

void MoveSurf::update_points(double fraction)
{
  /*
  int i;

  // update points by fraction of old to new move

  Surf::Point *pts = surf->pts;
  Surf::Point *p;

  for (i = 0; i < nread; i++) {
    if (pselect[readindex[i]] == 0) continue;
    p = &pts[readindex[i]];
    p->x[0] = oldcoord[i][0] + fraction * (newcoord[i][0]-oldcoord[i][0]);
    p->x[1] = oldcoord[i][1] + fraction * (newcoord[i][1]-oldcoord[i][1]);
    p->x[2] = oldcoord[i][2] + fraction * (newcoord[i][2]-oldcoord[i][2]);
  }
  */
}

/* ----------------------------------------------------------------------
   translate surf points in 2d
------------------------------------------------------------------------- */

void MoveSurf::translate_2d(double fraction, Surf::Line *origlines)
{
  double *p1,*p2,*op1,*op2;

  Surf::Line *lines = surf->lines;
  int nsurf = surf->nsurf;

  for (int i = 0; i < 2*nsurf; i++) pselect[i] = 0;

  double dx = fraction * delta[0];
  double dy = fraction * delta[1];

  for (int i = 0; i < nsurf; i++) {
    if (!(lines[i].mask & groupbit)) continue;
    p1 = lines[i].p1;
    p2 = lines[i].p2;
    op1 = origlines[i].p1;
    op2 = origlines[i].p2;

    p1[0] = op1[0] + dx;
    p1[1] = op1[1] + dy;
    pselect[2*i] = 1;

    p2[0] = op2[0] + dx;
    p2[1] = op2[1] + dy;
    pselect[2*i+1] = 1;
  }
}

/* ----------------------------------------------------------------------
   translate surf points in 3d
------------------------------------------------------------------------- */

void MoveSurf::translate_3d(double fraction, Surf::Tri *origtris)
{
  double *p1,*p2,*p3,*op1,*op2,*op3;

  Surf::Tri *tris = surf->tris;
  int nsurf = surf->nsurf;

  for (int i = 0; i < 3*nsurf; i++) pselect[i] = 0;

  double dx = fraction * delta[0];
  double dy = fraction * delta[1];
  double dz = fraction * delta[2];

  for (int i = 0; i < nsurf; i++) {
    if (!(tris[i].mask & groupbit)) continue;
    p1 = tris[i].p1;
    p2 = tris[i].p2;
    p3 = tris[i].p3;
    op1 = origtris[i].p1;
    op2 = origtris[i].p2;
    op3 = origtris[i].p3;

    p1[0] = op1[0] + dx;
    p1[1] = op1[1] + dy;
    p1[2] = op1[2] + dz;
    pselect[3*i] = 1;

    p2[0] = op2[0] + dx;
    p2[1] = op2[1] + dy;
    p2[2] = op2[2] + dz;
    pselect[3*i+1] = 1;

    p3[0] = op3[0] + dx;
    p3[1] = op3[1] + dy;
    p3[2] = op3[2] + dz;
    pselect[3*i+2] = 1;
  }
}

/* ----------------------------------------------------------------------
   rotate surf points in 2d
------------------------------------------------------------------------- */

void MoveSurf::rotate_2d(double fraction, Surf::Line *origlines)
{
  double *p1,*p2,*op1,*op2;
  double q[4],d[3],dnew[3];
  double rotmat[3][3];

  Surf::Line *lines = surf->lines;
  int nsurf = surf->nsurf;

  for (int i = 0; i < 2*nsurf; i++) pselect[i] = 0;

  double angle = fraction * theta;
  MathExtra::axisangle_to_quat(rvec,angle,q);
  MathExtra::quat_to_mat(q,rotmat);

  for (int i = 0; i < nsurf; i++) {
    if (!(lines[i].mask & groupbit)) continue;
    p1 = lines[i].p1;
    p2 = lines[i].p2;
    op1 = origlines[i].p1;
    op2 = origlines[i].p2;

    d[0] = op1[0] - origin[0];
    d[1] = op1[1] - origin[1];
    d[2] = op1[2] - origin[2];
    MathExtra::matvec(rotmat,d,dnew);
    p1[0] = dnew[0] + origin[0];
    p1[1] = dnew[1] + origin[1];
    pselect[2*i] = 1;

    d[0] = op2[0] - origin[0];
    d[1] = op2[1] - origin[1];
    d[2] = op2[2] - origin[2];
    MathExtra::matvec(rotmat,d,dnew);
    p2[0] = dnew[0] + origin[0];
    p2[1] = dnew[1] + origin[1];
    pselect[2*i+1] = 1;
  }
}

/* ----------------------------------------------------------------------
   rotate surf points in 3d
------------------------------------------------------------------------- */

void MoveSurf::rotate_3d(double fraction, Surf::Tri *origtris)
{
  double *p1,*p2,*p3,*op1,*op2,*op3;
  double q[4],d[3],dnew[3];
  double rotmat[3][3];

  Surf::Tri *tris = surf->tris;
  int nsurf = surf->nsurf;

  for (int i = 0; i < 3*nsurf; i++) pselect[i] = 0;

  double angle = fraction * theta;
  MathExtra::axisangle_to_quat(rvec,angle,q);
  MathExtra::quat_to_mat(q,rotmat);

  for (int i = 0; i < nsurf; i++) {
    if (!(tris[i].mask & groupbit)) continue;
    p1 = tris[i].p1;
    p2 = tris[i].p2;
    p3 = tris[i].p3;
    op1 = origtris[i].p1;
    op2 = origtris[i].p2;
    op3 = origtris[i].p3;

    d[0] = op1[0] - origin[0];
    d[1] = op1[1] - origin[1];
    d[2] = op1[2] - origin[2];
    MathExtra::matvec(rotmat,d,dnew);
    p1[0] = dnew[0] + origin[0];
    p1[1] = dnew[1] + origin[1];
    p1[2] = dnew[2] + origin[2];
    pselect[3*i] = 1;

    d[0] = op2[0] - origin[0];
    d[1] = op2[1] - origin[1];
    d[2] = op2[2] - origin[2];
    MathExtra::matvec(rotmat,d,dnew);
    p2[0] = dnew[0] + origin[0];
    p2[1] = dnew[1] + origin[1];
    p2[2] = dnew[2] + origin[2];
    pselect[3*i+1] = 1;

    d[0] = op3[0] - origin[0];
    d[1] = op3[1] - origin[1];
    d[2] = op3[2] - origin[2];
    MathExtra::matvec(rotmat,d,dnew);
    p3[0] = dnew[0] + origin[0];
    p3[1] = dnew[1] + origin[1];
    p3[2] = dnew[2] + origin[2];
    pselect[3*i+2] = 1;
  }
}

/* ----------------------------------------------------------------------
   add points in moved lines to hash
------------------------------------------------------------------------- */

void MoveSurf::connect_2d_pre()
{
  // hash for end points of moved lines
  // key = end point
  // value = global index (0 to 2*Nline-1) of the point
  // NOTE: could prealloc hash to correct size here

  hash = new MyHash();

  // add moved points to hash

  double *p1,*p2;
  OnePoint3d key;

  Surf::Line *lines = surf->lines;
  int nsurf = surf->nsurf;

  for (int i = 0; i < nsurf; i++) {
    if (!(lines[i].mask & groupbit)) continue;
    p1 = lines[i].p1;
    p2 = lines[i].p2;
    key.pt[0] = p1[0]; key.pt[1] = p1[1]; key.pt[2] = 0.0;
    if (hash->find(key) == hash->end()) (*hash)[key] = 2*i+0;
    key.pt[0] = p2[0]; key.pt[1] = p2[1]; key.pt[2] = 0.0;
    if (hash->find(key) == hash->end()) (*hash)[key] = 2*i+1;
  }
}

/* ----------------------------------------------------------------------
   move points in lines connected to line points that were moved
------------------------------------------------------------------------- */

void MoveSurf::connect_2d_post()
{
  // check if non-moved points are in hash
  // if so, set their coords to matching point
  // set pselect for newly moved points so remove_particles() will work

  int m,value,j,jwhich;
  double *p[2],*q;
  OnePoint3d key;

  Surf::Line *lines = surf->lines;
  int nsurf = surf->nsurf;

  for (int i = 0; i < nsurf; i++) {
    if (lines[i].mask & groupbit) continue;
    p[0] = lines[i].p1;
    p[1] = lines[i].p2;

    for (m = 0; m < 2; m++) {
      key.pt[0] = p[m][0]; key.pt[1] = p[m][1]; key.pt[2] = 0.0;
      if (hash->find(key) != hash->end()) {
        value = (*hash)[key];
        j = value/2;
        jwhich = value % 2;
        if (jwhich == 0) q = lines[j].p1;
        else q = lines[j].p2;
        p[m][0] = q[0];
        p[m][1] = q[1];
        if (m == 0) pselect[2*i] = 1;
        else pselect[2*i+1] = 1;
      }
    }
  }

  // free the hash

  delete hash;
}

/* ----------------------------------------------------------------------
   add points in moved triangles to hash
------------------------------------------------------------------------- */

void MoveSurf::connect_3d_pre()
{
  // hash for corner points of moved triangles
  // key = corner point
  // value = global index (0 to 3*Ntri-1) of the point
  // NOTE: could prealloc hash to correct size here

  hash = new MyHash();

  // add moved points to hash

  double *p1,*p2,*p3;
  OnePoint3d key;

  Surf::Tri *tris = surf->tris;
  int nsurf = surf->nsurf;

  for (int i = 0; i < nsurf; i++) {
    if (!(tris[i].mask & groupbit)) continue;
    p1 = tris[i].p1;
    p2 = tris[i].p2;
    p3 = tris[i].p3;
    key.pt[0] = p1[0]; key.pt[1] = p1[1]; key.pt[2] = p1[2];
    if (hash->find(key) == hash->end()) (*hash)[key] = 3*i+0;
    key.pt[0] = p2[0]; key.pt[1] = p2[1]; key.pt[2] = p2[2];
    if (hash->find(key) == hash->end()) (*hash)[key] = 3*i+1;
    key.pt[0] = p3[0]; key.pt[1] = p3[1]; key.pt[2] = p3[2];
    if (hash->find(key) == hash->end()) (*hash)[key] = 3*i+2;
  }
}

/* ----------------------------------------------------------------------
   move points in tris connected to tri points that were moved
------------------------------------------------------------------------- */

void MoveSurf::connect_3d_post()
{
  // check if non-moved points are in hash
  // if so, set their coords to matching point
  // set pselect for newly moved points so remove_particles() will work

  int m,value,j,jwhich;
  double *p[3],*q;
  OnePoint3d key;

  Surf::Tri *tris = surf->tris;
  int nsurf = surf->nsurf;

  for (int i = 0; i < nsurf; i++) {
    if (tris[i].mask & groupbit) continue;
    p[0] = tris[i].p1;
    p[1] = tris[i].p2;
    p[2] = tris[i].p3;

    for (m = 0; m < 3; m++) {
      key.pt[0] = p[m][0]; key.pt[1] = p[m][1]; key.pt[2] = p[m][2];
      if (hash->find(key) != hash->end()) {
        value = (*hash)[key];
        j = value/3;
        jwhich = value % 3;
        if (jwhich == 0) q = tris[j].p1;
        else if (jwhich == 1) q = tris[j].p2;
        else q = tris[j].p3;
        p[m][0] = q[0];
        p[m][1] = q[1];
        p[m][2] = q[2];
        if (m == 0) pselect[3*i] = 1;
        else if (m == 1) pselect[3*i+1] = 1;
        else pselect[3*i+2] = 1;
      }
    }
  }

  // free the hash

  delete hash;
}

/* ----------------------------------------------------------------------
   remove particles in any cell that is now INSIDE or contains moved surfs
   surfs that moved determined by pselect for any of its points
   reassign particles in split cells to sub cell owner
   compress particles if any flagged for deletion
   NOTE: doc this logic better
------------------------------------------------------------------------- */

bigint MoveSurf::remove_particles()
{
  int isurf,nsurf;
  surfint *csurfs;

  dim = domain->dimension;
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;
  int delflag = 0;

  for (int icell = 0; icell < nglocal; icell++) {

    // cell is inside surfs
    // remove particles in case it wasn't before

    if (cinfo[icell].type == INSIDE) {
      if (cinfo[icell].count) delflag = 1;
      particle->remove_all_from_cell(cinfo[icell].first);
      cinfo[icell].count = 0;
      cinfo[icell].first = -1;
      continue;
    }

    // cell has surfs or is split
    // if m < nsurf, loop over csurfs did not finish
    // which means cell contains a moved surf, so delete all its particles

    if (cells[icell].nsurf && cells[icell].nsplit >= 1) {
      nsurf = cells[icell].nsurf;
      csurfs = cells[icell].csurfs;

      int m;
      if (dim == 2) {
        for (m = 0; m < nsurf; m++) {
          isurf = csurfs[m];
          if (pselect[2*isurf]) break;
          if (pselect[2*isurf+1]) break;
        }
      } else {
        for (m = 0; m < nsurf; m++) {
          isurf = csurfs[m];
          if (pselect[3*isurf]) break;
          if (pselect[3*isurf+1]) break;
          if (pselect[3*isurf+2]) break;
        }
      }

      if (m < nsurf) {
        if (cinfo[icell].count) delflag = 1;
        particle->remove_all_from_cell(cinfo[icell].first);
        cinfo[icell].count = 0;
        cinfo[icell].first = -1;
      }
    }

    if (cells[icell].nsplit > 1)
      grid->assign_split_cell_particles(icell);
  }

  int nlocal_old = particle->nlocal;
  if (delflag) particle->compress_rebalance();
  bigint delta = nlocal_old - particle->nlocal;
  bigint ndeleted;
  MPI_Allreduce(&delta,&ndeleted,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  return ndeleted;
}
