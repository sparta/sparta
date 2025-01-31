/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
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
#include "custom.h"
#include "domain.h"
#include "comm.h"
#include "particle.h"
#include "grid.h"
#include "surf.h"
#include "mixture.h"
#include "region.h"
#include "input.h"
#include "variable.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{SET,FILESTYLE,REMOVE};
enum{EQUAL,PARTICLE,GRID,SURF};
enum{INT,DOUBLE};                       // several files

#define MAXLINE 256
#define CHUNK 4     // make this larger after debugging

/* ---------------------------------------------------------------------- */

Custom::Custom(SPARTA *sparta) : Pointers(sparta) {}

/* ---------------------------------------------------------------------- */

void Custom::command(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR,"Illegal custom command");

  // mode

  if (strcmp(arg[0],"particle") == 0) mode = PARTICLE;
  else if (strcmp(arg[0],"grid") == 0) mode = GRID;
  else if (strcmp(arg[0],"surf") == 0) mode = SURF;
  else error->all(FLERR,"Illegal custom command");

  if (mode == PARTICLE && !particle->exist)
    error->all(FLERR,"Cannot use custom particle before particles are defined");
  if (mode == GRID && !grid->exist)
    error->all(FLERR,"Cannot use custom grid before a grid is defined");
  if (mode == SURF && !surf->exist)
    error->all(FLERR,"Cannot use custom surf before surfaces are defined");

  // attribute name

  int n = strlen(arg[1]) + 1;
  aname = new char[n];
  strcpy(aname,arg[1]);

  char *ptr = strchr(aname,'[');
  if (ptr) {
    if (aname[strlen(aname)-1] != ']')
      error->all(FLERR,"Custom attribute name is invalid");
    ccol = atoi(ptr+1);
    *ptr = '\0';
  } else ccol = 0;

  // action

  if (strcmp(arg[2],"set") == 0) action = SET;
  else if (strcmp(arg[2],"file") == 0) action = FILESTYLE;
  else if (strcmp(arg[2],"remove") == 0) action = REMOVE;
  else error->all(FLERR,"Illegal set command");

  // remove a custom attribute and return

  if (action == REMOVE) {
    if (narg > 3) error->all(FLERR,"Illegal custom command");
    if (mode == PARTICLE) {
      int index = particle->find_custom(aname);
      if (index < 0) error->all(FLERR,"Custom attribute name does not exist");
      particle->remove_custom(index);
    } else if (mode == GRID) {
      int index = grid->find_custom(aname);
      if (index < 0) error->all(FLERR,"Custom attribute name does not exist");
      grid->remove_custom(index);
    } else if (mode == SURF) {
      int index = surf->find_custom(aname);
      if (index < 0) error->all(FLERR,"Custom attribute name does not exist");
      surf->remove_custom(index);
    }

    delete [] aname;

    return;
  }

  // read args to reset a custom attribute via a variable

  int iarg;
  vname = fname = NULL;
  
  if (action == SET) {

    if (narg < 6) error->all(FLERR,"Illegal custom command");
    
    // variable name

    variable = input->variable;

    if (strncmp(arg[3],"v_",2) == 0) {
      int n = strlen(arg[3]);
      vname = new char[n];
      strcpy(vname,&arg[3][2]);

      vindex = variable->find(vname);
      if (vindex < 0) error->all(FLERR,"Custom variable name does not exist");
      if (variable->equal_style(vindex)) vstyle = EQUAL;
      else if (variable->particle_style(vindex)) vstyle = PARTICLE;
      else if (variable->grid_style(vindex)) vstyle = GRID;
      else if (variable->surf_style(vindex)) vstyle = SURF;
      else error->all(FLERR,"Custom variable style is invalid");
      if (vstyle != EQUAL && vstyle != mode)
        error->all(FLERR,"Custom variable style is invalid");
      
    } else error->all(FLERR,"Custom variable name is invalid");

    // mixture or group ID

    if (mode == PARTICLE) {
      int imix = particle->find_mixture(arg[4]);
      if (imix < 0) error->all(FLERR,"Custom mixture ID does not exist");
      mixture = particle->mixture[imix];
      mixture->init();
    } else if (mode == GRID) {
      int igroup = grid->find_group(arg[4]);
      if (igroup < 0) error->all(FLERR,"Custom grid group ID does not exist");
      groupbit = grid->bitmask[igroup];
    } else if (mode == SURF) {
      int igroup = surf->find_group(arg[4]);
      if (igroup < 0) error->all(FLERR,"Custom surf group ID does not exist");
      groupbit = surf->bitmask[igroup];
    }

    // region ID

    if (strcmp(arg[5],"NULL") == 0) region = NULL;
    else {
      int iregion = domain->find_region(arg[5]);
      if (iregion < 0) error->all(FLERR,"Custom region ID does not exist");
      region = domain->regions[iregion];
    }

    iarg = 6;
  }

  // read args to reset a custom attribute via a file

  if (action == FILESTYLE) {

    if (narg < 4) error->all(FLERR,"Illegal custom command");
    if (mode == PARTICLE)
      error->all(FLERR,"Custom command cannot use action file with style particle");

    // file name
    
    int n = strlen(arg[3]);
    fname = new char[n];
    strcpy(fname,arg[3]);

    iarg = 4;
  }  

  // optional keywords

  ctype = DOUBLE;
  csize = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid custom command");
      if (strcmp(arg[iarg+1],"int") == 0) ctype = INT;
      else if (strcmp(arg[iarg+1],"float") == 0) ctype = DOUBLE;
      else error->all(FLERR,"Invalid custom command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"size") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid custom command");
      csize = input->inumeric(FLERR,arg[iarg+1]);
      if (csize < 0) error->all(FLERR,"Invalid custom command");
      iarg += 2;
    } else error->all(FLERR,"Invalid custom command");
  }

  // cindex = index of existing custom attribute
  // otherwise create custom attribute if it does not exist
  // add_custom() initializes all values to zero
  
  if (mode == PARTICLE) {
    cindex = particle->find_custom(aname);
    if (cindex >= 0) {
      ctype = particle->etype[cindex];
      csize = particle->esize[cindex];
    } else {
      cindex = particle->add_custom(aname,ctype,csize);
      
    }

  } else if (mode == GRID) {
    cindex = grid->find_custom(aname);
    if (cindex >= 0) {
      ctype = grid->etype[cindex];
      csize = grid->esize[cindex];
    } else {
      cindex = grid->add_custom(aname,ctype,csize);
    }

  } else if (mode == SURF) {
    cindex = surf->find_custom(aname);
    if (cindex >= 0) {
      ctype = surf->etype[cindex];
      csize = surf->esize[cindex];
    } else {
      cindex = surf->add_custom(aname,ctype,csize);
    }
  }

  // error check that name or name[] matches new or existing attribute

  int eflag = 0;
  if (csize == 0 && ccol != 0) eflag = 1;
  if (csize != 0 && ccol == 0) eflag = 1;
  if (csize != 0 && ccol > csize) eflag = 1;
  if (eflag) error->all(FLERR,"Custom name does not match new or existing custom attribute");

  // for action SET, evaluate variable
  // store result as floating point scalar or vector
  // count = # of changed attributes

  int count;
  double *vector = NULL;

  if (action == SET) {
    double scalar = 0.0;

    if (vstyle == EQUAL) {
      scalar = variable->compute_equal(vindex);
    } else if (vstyle == PARTICLE) {
      memory->create(vector,particle->nlocal,"custom:vector");
      variable->compute_particle(vindex,vector,1,0);
    } else if (vstyle == GRID) {
      memory->create(vector,grid->nlocal,"custom:vector");
      variable->compute_grid(vindex,vector,1,0);
    } else if (vstyle == SURF) {
      memory->create(vector,surf->nown,"custom:vector");
      variable->compute_surf(vindex,vector,1,0);
    }

    // assign value(s) to custom attribute
    // convert to integer if necessary
    // no assignment if particle/grid/surf not in mixture or group or region
    
    if (mode == PARTICLE) count = set_particle(scalar,vector);
    else if (mode == GRID) count = set_grid(scalar,vector);
    else if (mode == SURF) count = set_surf(scalar,vector);
  }

  // for action FILESTYLE, read file and set attributes
  // count = # of changed attributes

  if (action = FILESTYLE) {
    count = read_file(mode,cindex,ctype,csize,fname,0);
  }
    
  // print stats

  bigint bcount = count;
  bigint bcountall;
  MPI_Allreduce(&bcount,&bcountall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  const char *mname;
  if (mode == PARTICLE) mname = "particle";
  else if (mode == GRID) mname = "grid";
  else if (mode == SURF) mname = "surf";

  if (comm->me == 0) {
    if (screen)
      fprintf(screen,"Custom %s %s attributes set = " BIGINT_FORMAT "\n",
              mname,aname,bcountall);
    if (logfile)
      fprintf(logfile,"Custom %s %s attributes set = " BIGINT_FORMAT "\n",
              mname,aname,bcountall);
  }

  // clean up

  delete [] aname;
  delete [] vname;
  delete [] fname;
  memory->destroy(vector);
}

/* ---------------------------------------------------------------------- */

int Custom::set_particle(double scalar, double *vector)
{
  Particle::OnePart *particles = particle->particles;
  int *species2species = mixture->species2species;
  int nlocal = particle->nlocal;

  int *choose;
  memory->create(choose,nlocal,"set:choose");
  memset(choose,0,nlocal*sizeof(int));

  int flag;

  int count = 0;
  for (int i = 0; i < nlocal; i++) {
    flag = 1;
    if (species2species[particles[i].ispecies] < 0) flag = 0;
    if (flag && region) {
      if (!region->inside(particles[i].x)) flag = 0;
    }
    if (!flag) continue;

    choose[i] = 1;
    count++;
  }

  // set custom values via scalar or vector

  if (ctype == INT) {
    int iscalar = static_cast<int> (scalar);

    if (csize == 0) {
      int *cvector = particle->eivec[particle->ewhich[cindex]];
      if (vector) {
	for (int i = 0 ; i < nlocal; i++) {
	  if (choose[i]) cvector[i] = static_cast<int> (vector[i]);
	}
      } else {
	for (int i = 0 ; i < nlocal; i++) {
	  if (choose[i]) cvector[i] = iscalar;
	}
      }

    } else {
      int **carray = particle->eiarray[particle->ewhich[cindex]];
      ccol--;
      if (vector) {
	for (int i = 0 ; i < nlocal; i++) {
	  if (choose[i]) carray[i][ccol] = static_cast<int> (vector[i]);
	}
      } else {
	for (int i = 0 ; i < nlocal; i++) {
	  if (choose[i]) carray[i][ccol] = iscalar;
	}
      }
    }

  } else if (ctype == DOUBLE) {
    if (csize == 0) {
      double *cvector = particle->edvec[particle->ewhich[cindex]];
      if (vector) {
	for (int i = 0 ; i < nlocal; i++) {
	  if (choose[i]) cvector[i] = vector[i];
	}
      } else {
	for (int i = 0 ; i < nlocal; i++) {
	  if (choose[i]) cvector[i] = scalar;
	}
      }

    } else {
      double **carray = particle->edarray[particle->ewhich[cindex]];
      ccol--;
      if (vector) {
	for (int i = 0 ; i < nlocal; i++) {
	  if (choose[i]) carray[i][ccol] = vector[i];
	}
      } else {
	for (int i = 0 ; i < nlocal; i++) {
	  if (choose[i]) carray[i][ccol] = scalar;
	}
      }
    }
  }

  memory->destroy(choose);

  return count;
}

/* ----------------------------------------------------------------------
   read a custom attribute file for fmode = GRID or SURF
   assign values to custom grid or surf vectors/arrays
   return count of attributes assigned by this proc
   external = 1 when called from fix custom
---------------------------------------------------------------------- */

int Custom::read_file(int fmode, int cindex, int ctype, int csize, char *filename, int external)
{
  // setup read buffers
  // NOTE: should these be created ?
  
  char *line = new char[MAXLINE];
  char *buffer = new char[CHUNK*MAXLINE];

  // set ivec,dvec,iarray,darray pointers
  // NOTE: what about when setting one column of an array
  
  int *ivec;
  double *dvec;
  int **iarray;
  double **darray;
  
  // nwords_required = # of words per line

  //int nwords_required = 1 + nvalues_custom;   // NOTE: need to set this
  int nwords_required = 2;
  
  // read file

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  int me = comm->me;
  int nprocs = comm->nprocs;
  FILE *fp;

  // NOTE: support compressed files like read_grid ?
  
  if (me == 0) {
    if (!external && screen)
      fprintf(screen,"Reading custom grid file ... %s\n",filename);
    fp = fopen(filename,"r");
  }

  // read header portion of file
  // nfile = count of attribute lines in file

  // NOTE: allow for nfile to be a bigint ?
  
  int nfile;
  
  if (me == 0) {
    char *eof,*ptr;
    
    while (1) {
      eof = fgets(line,MAXLINE,fp);
      if (eof == NULL) error->one(FLERR,"Unexpected end of custom attribute file");

      // trim anything from '#' onward
      // if line is blank, continue
      // else break and read nfile

      if ((ptr = strchr(line,'#'))) *ptr = '\0';
      if (strspn(line," \t\n\r") == strlen(line)) continue;
      break;
    }

    sscanf(line,"%d",&nfile);
    //sscanf(line,BIGINT_FORMAT,&nfile);
  }

  MPI_Bcast(&nfile,1,MPI_INT,0,world);
  
  // read and broadcast one CHUNK of lines at a time

  bigint nread = 0;
  int i,m,nchunk;
  char *next,*idptr;
  cellint id;
  
  bigint fcount = 0;
  while (nread < nfile) {
    if (nfile-nread > CHUNK) nchunk = CHUNK;
    else nchunk = nfile-nread;
    if (me == 0) {
      char *eof;
      m = 0;
      for (i = 0; i < nchunk; i++) {
        eof = fgets(&buffer[m],MAXLINE,fp);
        if (eof == NULL) error->one(FLERR,"Unexpected end of custom attribute file");
        m += strlen(&buffer[m]);
      }
      if (buffer[m-1] != '\n') strcpy(&buffer[m++],"\n");
      m++;
    }
    MPI_Bcast(&m,1,MPI_INT,0,world);
    MPI_Bcast(buffer,m,MPI_CHAR,0,world);

    // add occasional barrier to prevent issues from having too many
    //  outstanding MPI recv requests (from the broadcast above)

    if (fcount % 1024 == 0) MPI_Barrier(world);

    // process nchunk lines and assign attribute value(s) if I own grid cell ID

    for (i = 0; i < nchunk; i++) {
      next = strchr(buffer,'\n');

      *next = '\0';
      int nwords = input->count_words(buffer);
      *next = '\n';
      
      if (nwords != nwords_required)
	error->all(FLERR,"Incorrect line format in custom attribute file");
      
      idptr = strtok(buffer," \t\n\r\f");
      id = ATOCELLINT(idptr);
      if (id < 0) error->all(FLERR,"Invalid cell ID in custom attribute grid file");

      // NOTE: need to look up grid cell ID on my proc to get local index
      int index;
      
      if (ctype == INT) {
        if (csize == 0)
          ivec[index] = input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
        else
          for (int iv = 0; iv < csize; iv++)
            iarray[index][iv] = input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
      } else if (ctype == DOUBLE) {
        if (csize == 0)
          dvec[index] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
        else
          for (int iv = 0; iv < csize; iv++)
            darray[index][iv] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
      }
    }

    // increment nread and fcount and continue to next chunck
    
    nread += nchunk;
    fcount++;
  }
  
  // close file

  if (me == 0) {
    //if (compressed) pclose(fp);
    //else fclose(fp);
    fclose(fp);
  }

  
  int count = 0;


  return count;

  // free read buffers
  
  delete [] line;
  delete [] buffer;

  // NOTE: print stats if not external ?
}

/* ---------------------------------------------------------------------- */

int Custom::set_grid(double scalar, double *vector)
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  int *choose;
  memory->create(choose,nglocal,"set:choose");
  memset(choose,0,nglocal*sizeof(int));

  int flag;
  double point[3];

  int count = 0;
  for (int i = 0; i < nglocal; i++) {
    flag = 1;
    if (!(cinfo[i].mask & groupbit)) flag = 0;
    if (flag && region) {
      point[0] = 0.5 * (cells[i].lo[0] + cells[i].hi[0]);
      point[1] = 0.5 * (cells[i].lo[1] + cells[i].hi[1]);
      point[2] = 0.5 * (cells[i].lo[2] + cells[i].hi[2]);
      if (!region->inside(point)) flag = 0;
    }
    if (!flag) continue;

    choose[i] = 1;
    count++;
  }

  // set custom values via scalar or vector

  if (ctype == INT) {
    int iscalar = static_cast<int> (scalar);

    if (csize == 0) {
      int *cvector = grid->eivec[grid->ewhich[cindex]];
      if (vector) {
	for (int i = 0 ; i < nglocal; i++) {
	  if (choose[i]) cvector[i] = static_cast<int> (vector[i]);
	}
      } else {
	for (int i = 0 ; i < nglocal; i++) {
	  if (choose[i]) cvector[i] = iscalar;
	}
      }

    } else {
      int **carray = grid->eiarray[grid->ewhich[cindex]];
      ccol--;
      if (vector) {
	for (int i = 0 ; i < nglocal; i++) {
	  if (choose[i]) carray[i][ccol] = static_cast<int> (vector[i]);
	}
      } else {
	for (int i = 0 ; i < nglocal; i++) {
	  if (choose[i]) carray[i][ccol] = iscalar;
	}
      }
    }

  } else if (ctype == DOUBLE) {
    if (csize == 0) {
      double *cvector = grid->edvec[grid->ewhich[cindex]];
      if (vector) {
	for (int i = 0 ; i < nglocal; i++) {
	  if (choose[i]) cvector[i] = vector[i];
	}
      } else {
	for (int i = 0 ; i < nglocal; i++) {
	  if (choose[i]) cvector[i] = scalar;
	}
      }

    } else {
      double **carray = grid->edarray[grid->ewhich[cindex]];
      ccol--;
      if (vector) {
	for (int i = 0 ; i < nglocal; i++) {
	  if (choose[i]) carray[i][ccol] = vector[i];
	}
      } else {
	for (int i = 0 ; i < nglocal; i++) {
	  if (choose[i]) carray[i][ccol] = scalar;
	}
      }
    }
  }

  memory->destroy(choose);

  return count;
}

/* ---------------------------------------------------------------------- */

int Custom::set_surf(double scalar, double *vector)
{
  int dim = domain->dimension;
  int distributed = surf->distributed;

  Surf::Line *lines;
  Surf::Tri *tris;
  int start,stop,skip;

  if (!distributed) {
    lines = surf->lines;
    tris = surf->tris;
    start = comm->me;
    stop = surf->nlocal;
    skip = comm->nprocs;
  } else {
    lines = surf->mylines;
    tris = surf->mytris;
    start = 0;
    stop = surf->nown;
    skip = 1;
  }

  int nsown = surf->nown;
  int *choose;
  memory->create(choose,nsown,"set:choose");
  memset(choose,0,nsown*sizeof(int));

  int flag;
  double point[3];

  int count = 0;
  for (int i = start ; i < stop; i += skip) {
    flag = 1;
    if (dim == 2) {
      if (!(lines[i].mask & groupbit)) flag = 0;
    } else {
      if (!(tris[i].mask & groupbit)) flag = 0;
    }
    if (flag && region) {
      if (dim == 2) {
	point[0] = 0.5 * (lines[i].p1[0] + lines[i].p2[0]);
	point[1] = 0.5 * (lines[i].p1[0] + lines[i].p2[0]);
	point[2] = 0.0;
      } else {
	point[0] = MathConst::THIRD *
          (tris[i].p1[0] + tris[i].p2[0] + tris[i].p3[0]);
	point[1] = MathConst::THIRD *
          (tris[i].p1[1] + tris[i].p2[1] + tris[i].p3[1]);
	point[2] = MathConst::THIRD *
          (tris[i].p1[2] + tris[i].p2[2] + tris[i].p3[2]);
      }
      if (!region->inside(point)) flag = 0;
    }
    if (!flag) continue;

    choose[count++] = 1;
  }

  // set custom values via scalar or vector

  if (ctype == INT) {
    int iscalar = static_cast<int> (scalar);

    if (csize == 0) {
      int *cvector = surf->eivec[surf->ewhich[cindex]];
      if (vector) {
	for (int i = 0 ; i < nsown; i++) {
	  if (choose[i]) cvector[i] = static_cast<int> (vector[i]);
	}
      } else {
	for (int i = 0 ; i < nsown; i++) {
	  if (choose[i]) cvector[i] = iscalar;
	}
      }

    } else {
      int **carray = surf->eiarray[surf->ewhich[cindex]];
      ccol--;
      if (vector) {
	for (int i = 0 ; i < nsown; i++) {
	  if (choose[i]) carray[i][ccol] = static_cast<int> (vector[i]);
	}
      } else {
	for (int i = 0 ; i < nsown; i++) {
	  if (choose[i]) carray[i][ccol] = iscalar;
	}
      }
    }

  } else if (ctype == DOUBLE) {
    if (csize == 0) {
      double *cvector = surf->edvec[surf->ewhich[cindex]];
      if (vector) {
	for (int i = 0 ; i < nsown; i++) {
	  if (choose[i]) cvector[i] = vector[i];
	}
      } else {
	for (int i = 0 ; i < nsown; i++) {
	  if (choose[i]) cvector[i] = scalar;
	}
      }

    } else {
      double **carray = surf->edarray[surf->ewhich[cindex]];
      ccol--;
      if (vector) {
	for (int i = 0 ; i < nsown; i++) {
	  if (choose[i]) carray[i][ccol] = vector[i];
	}
      } else {
	for (int i = 0 ; i < nsown; i++) {
	  if (choose[i]) carray[i][ccol] = scalar;
	}
      }
    }
  }

  memory->destroy(choose);

  return count;
}
