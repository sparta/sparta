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
#include "mixture.h"
#include "region.h"
#include "input.h"
#include "surf.h"
#include "update.h"
#include "variable.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{CREATE,REMOVE,SET,FILESTYLE,FILECOARSE};
enum{EQUAL,PARTICLE,GRID,SURF};
enum{INT,DOUBLE};                       // several files
enum{TEXT,BINARY};

#define MAXLINE 256
#define CHUNK 4     // NOTE: make this larger after debugging
#define BIG 1.0e20
#define MAXTIE 2    // NOTE: make this larger after debugging
#define EPSCUT 1.0e-6

/* ---------------------------------------------------------------------- */

Custom::Custom(SPARTA *sparta) : Pointers(sparta) {}

/* ---------------------------------------------------------------------- */

Custom::~Custom()
{
  // delete data in Action list which consumes memory
  
  for (int i = 0; i < naction; i++) {
    int action = actions[i].action;
    if (action == FILESTYLE) {
      delete [] actions[i].fname;
      delete [] actions[i].cindex_file;
      delete [] actions[i].ctype_file;
      delete [] actions[i].csize_file;
      delete [] actions[i].ccol_file;
    }
  }
  
  memory->sfree(actions);
}

/* ---------------------------------------------------------------------- */

void Custom::command(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR,"Illegal custom command");

  // process mode and sequence of actions
  // perform them one by one
  // final arg = 0 for calling from Custom
  
  bigint count = process_actions(narg,arg,0);

  // print stats

  bigint countall;
  MPI_Allreduce(&count,&countall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  const char *mname;
  if (mode == PARTICLE) mname = "particle";
  else if (mode == GRID) mname = "grid";
  else if (mode == SURF) mname = "surf";

  if (comm->me == 0) {
    if (screen)
      fprintf(screen,"Custom %s attributes set = " BIGINT_FORMAT "\n",
              mname,countall);
    if (logfile)
      fprintf(logfile,"Custom %s attributes set = " BIGINT_FORMAT "\n",
              mname,countall);
  }
}

/* ----------------------------------------------------------------------
   parse mode and list of actions from input script command
   external = 0 if called from Custom, 1 if called from FixCustom
   if external = 0
     invoke each action as soon as parsed
     b/c an action may depend on previous actions
   if external = 1
     CREATE and REMOVE actions are not allowed
     store action info in Action list for later invocation by FixCustom
---------------------------------------------------------------------- */

bigint Custom::process_actions(int narg, char **arg, int external)
{
  // mode
  
  if (strcmp(arg[0],"particle") == 0) mode = PARTICLE;
  else if (strcmp(arg[0],"grid") == 0) mode = GRID;
  else if (strcmp(arg[0],"surf") == 0) mode = SURF;
  else error->all(FLERR,"Illegal custom command");

  if (mode == PARTICLE && !particle->exist)
    error->all(FLERR,"Cannot use custom particle command before particles are defined");
  if (mode == GRID && !grid->exist)
    error->all(FLERR,"Cannot use custom grid command before a grid is defined");
  if (mode == SURF && !surf->exist)
    error->all(FLERR,"Cannot use custom surf command before surfaces are defined");

  // loop over actions

  naction = 0;
  actions = NULL;
  bigint count = 0;
  
  int iarg = 1;
    
  while (iarg < narg) {

    // grow list of Actions for FixCustom

    if (external)
      actions = (Action *) memory->srealloc(actions,(naction+1)*sizeof(Action),
                                            "custom:actions");

    int action;
    if (strcmp(arg[iarg],"create") == 0) action = CREATE;
    else if (strcmp(arg[iarg],"remove") == 0) action = REMOVE;
    else if (strcmp(arg[iarg],"set") == 0) action = SET;
    else if (strcmp(arg[iarg],"file") == 0) action = FILESTYLE;
    else if (strcmp(arg[iarg],"file/coarse") == 0) action = FILECOARSE;
    else error->all(FLERR,"Illegal custom command action");

    // create new custom attribute
    
    if (action == CREATE) {

      if (external) error->all(FLERR,"Fix custom cannot use create action");
      if (iarg+4 > narg) error->all(FLERR,"Illegal custom command");
      
      int n = strlen(arg[iarg+1]) + 1;
      char *aname = new char[n];
      strcpy(aname,arg[iarg+1]);
      int ccol = attribute_bracket(aname);
      if (ccol) error->all(FLERR,"Illegal custom attribute syntax");

      int cindex;
      if (mode == PARTICLE) cindex = particle->find_custom(aname);
      else if (mode == GRID) cindex = grid->find_custom(aname);
      else if (mode == SURF) cindex = surf->find_custom(aname);
      if (cindex >= 0) error->all(FLERR,"Custom attribute name already exists");

      int ctype;
      if (strcmp(arg[iarg+2],"int") == 0) ctype = INT;
      else if (strcmp(arg[iarg+2],"float") == 0) ctype = DOUBLE;
      else error->all(FLERR,"Illegal custom create datatype");

      int csize = input->inumeric(FLERR,arg[iarg+3]);
      if (csize < 0) error->all(FLERR,"Invalid custom create size");

      if (mode == PARTICLE) particle->add_custom(aname,ctype,csize);
      else if (mode == GRID) grid->add_custom(aname,ctype,csize);
      else if (mode == SURF) surf->add_custom(aname,ctype,csize);

      delete [] aname;
      iarg += 4;
      
    // remove a custom attribute

    } else if (action == REMOVE) {
      
      if (external) error->all(FLERR,"Fix custom cannot use remove action");
      if (iarg+2 > narg) error->all(FLERR,"Illegal custom command");

      int n = strlen(arg[iarg+1]) + 1;
      char *aname = new char[n];
      strcpy(aname,arg[iarg+1]);
      int ccol = attribute_bracket(aname);
      if (ccol) error->all(FLERR,"Illegal custom attribute syntax");

      int cindex;
      if (mode == PARTICLE) cindex = particle->find_custom(aname);
      else if (mode == GRID) cindex = grid->find_custom(aname);
      else if (mode == SURF) cindex = surf->find_custom(aname);
      if (cindex < 0) error->all(FLERR,"Custom attribute name does not exist");

      if (mode == PARTICLE) particle->remove_custom(cindex);
      else if (mode == GRID) grid->remove_custom(cindex);
      else if (mode == SURF) surf->remove_custom(cindex);

      delete [] aname;
      iarg += 2;

    // set a custom vector or column of custom array via a variable

    } else if (action == SET) {

      if (iarg+5 > narg) error->all(FLERR,"Illegal custom command");

      int n = strlen(arg[iarg+1]) + 1;
      char *aname = new char[n];
      strcpy(aname,arg[iarg+1]);
      int ccol = attribute_bracket(aname);

      int cindex,ctype,csize;
      if (mode == PARTICLE) {
        cindex = particle->find_custom(aname);
        if (cindex < 0) error->all(FLERR,"Custom attribute name does not exist");
        ctype = particle->etype[cindex];
        csize = particle->esize[cindex];
      } else if (mode == GRID) {
        cindex = grid->find_custom(aname);
        if (cindex < 0) error->all(FLERR,"Custom attribute name does not exist");
        ctype = grid->etype[cindex];
        csize = grid->esize[cindex];
      } else if (mode == SURF) {
        cindex = surf->find_custom(aname);
        if (cindex < 0) error->all(FLERR,"Custom attribute name does not exist");
        ctype = surf->etype[cindex];
        csize = surf->esize[cindex];
      }

      if (csize && ccol == 0)
        error->all(FLERR,"Custom attribute array requires bracketed index");
      if (csize == 0 && ccol)
        error->all(FLERR,"Custom attribute vector cannot use bracketed index");
    
      // variable name

      if (strncmp(arg[iarg+2],"v_",2) != 0)
        error->all(FLERR,"Custom variable name is invalid");
        
      n = strlen(arg[iarg+2]);
      char *vname = new char[n];
      strcpy(vname,&arg[iarg+2][2]);
      
      Variable *variable = input->variable;
      int vindex = variable->find(vname);
      if (vindex < 0) error->all(FLERR,"Custom variable name does not exist");

      int vstyle;
      if (variable->equal_style(vindex)) vstyle = EQUAL;
      else if (variable->particle_style(vindex)) vstyle = PARTICLE;
      else if (variable->grid_style(vindex)) vstyle = GRID;
      else if (variable->surf_style(vindex)) vstyle = SURF;
      else error->all(FLERR,"Custom variable style is invalid");
      if (vstyle != EQUAL && vstyle != mode)
        error->all(FLERR,"Custom variable style is invalid");

      // mixture or group ID

      Mixture *mixture;
      int groupbit;
      
      if (mode == PARTICLE) {
        int imix = particle->find_mixture(arg[iarg+3]);
        if (imix < 0) error->all(FLERR,"Custom mixture ID does not exist");
        mixture = particle->mixture[imix];
        mixture->init();
      } else if (mode == GRID) {
        int igroup = grid->find_group(arg[iarg+3]);
        if (igroup < 0) error->all(FLERR,"Custom grid group ID does not exist");
        groupbit = grid->bitmask[igroup];
      } else if (mode == SURF) {
        int igroup = surf->find_group(arg[iarg+3]);
        if (igroup < 0) error->all(FLERR,"Custom surf group ID does not exist");
        groupbit = surf->bitmask[igroup];
      }

      // region ID

      Region *region;
      if (strcmp(arg[iarg+4],"NULL") == 0) region = NULL;
      else {
        int iregion = domain->find_region(arg[iarg+4]);
        if (iregion < 0) error->all(FLERR,"Custom region ID does not exist");
        region = domain->regions[iregion];
      }

      // if not external: invoke action now
      // else: store info in Action list for FixCustom

      if (!external)
        action_set(vstyle,vindex,cindex,ctype,csize,ccol,groupbit,mixture,region);
      else {
        actions[naction].action = action;
        actions[naction].vstyle = vstyle;
        actions[naction].vindex = vindex;
        actions[naction].cindex = cindex;
        actions[naction].ctype = ctype;
        actions[naction].csize = csize;
        actions[naction].ccol = ccol;
        actions[naction].groupbit = groupbit;
        actions[naction].mixture = mixture;
        actions[naction].region = region;
        naction++;
      }
      
      delete [] aname;
      delete [] vname;
      iarg += 5;
      
    // set multiple custom attributes from a file

    } else if (action == FILESTYLE) {
      
      if (iarg+3 > narg) error->all(FLERR,"Illegal custom command");

      if (mode == PARTICLE)
        error->all(FLERR,"Custom command cannot use action file with style particle");
      
      // file name
      // if external check that fname has wildcard char
      
      int n = strlen(arg[iarg+1]) + 1;
      char *fname = new char[n];
      strcpy(fname,arg[iarg+1]);

      if (external && strchr(fname,'*') == NULL)
        error->all(FLERR,"Fix custom filename must have * wildcard char");

      // colcount = # of attribute values per file line
      
      int colcount = input->inumeric(FLERR,arg[iarg+2]);
      if (colcount < 1) error->all(FLERR,"Custom command file column count is invalid");
      if (iarg+3+colcount > narg) error->all(FLERR,"Illegal custom command");

      // create vectors of attribute name settings

      int *cindex = new int[colcount];
      int *ctype = new int[colcount];
      int *csize = new int[colcount];
      int *ccol = new int[colcount];

      for (int i = 0; i < colcount; i++) {
        int n = strlen(arg[iarg+3+i]) + 1;  
        char *aname = new char[n];
        strcpy(aname,arg[iarg+3+i]);
        ccol[i] = attribute_bracket(aname);
        if (mode == GRID) {
          cindex[i] = grid->find_custom(aname);
          if (cindex[i] < 0)
            error->all(FLERR,"Custom attribute name does not exist");
          ctype[i] = grid->etype[cindex[i]];
          csize[i] = grid->esize[cindex[i]];
        } else if (mode == SURF) {
          cindex[i] = surf->find_custom(aname);
          if (cindex[i] < 0)
            error->all(FLERR,"Custom attribute name does not exist");
          ctype[i] = surf->etype[cindex[i]];
          csize[i] = surf->esize[cindex[i]];
        }
        if (csize[i] && ccol[i] == 0)
          error->all(FLERR,"Custom attribute array requires bracketed index");
        if (csize[i] == 0 && ccol[i])
          error->all(FLERR,"Custom attribute vector cannot use bracketed index");
        delete [] aname;
      }

      // if not external: invoke action now
      // else: store info in Action list for FixCustom
      // for mode = GRID
      //   set estatus of all changed custom vecs/arrays to 1
      //   b/c file read stores values for owned+ghost cells
      // for mode = SURF
      //   set estatus of all changed custom vecs/arrays to 0

      if (!external) {
        if (comm->me == 0 && screen)
          fprintf(screen,"Reading custom file %s ...\n",fname);
        count += read_file(mode,colcount,cindex,ctype,csize,ccol,fname);

        if (mode == GRID)
          for (int i = 0; i < colcount; i++)
            grid->estatus[cindex[i]] = 1;
        else if (mode == SURF)
          for (int i = 0; i < colcount; i++)
            surf->estatus[cindex[i]] = 0;

        delete [] fname;
        delete [] cindex;
        delete [] ctype;
        delete [] csize;
        delete [] ccol;

      } else {
        actions[naction].action = action;
        actions[naction].fname = fname;
        actions[naction].colcount = colcount;
        actions[naction].cindex_file = cindex;
        actions[naction].ctype_file = ctype;
        actions[naction].csize_file = csize;
        actions[naction].ccol_file = ccol;
        naction++;
      }
      
      iarg += 3 + colcount;

    // set multiple custom attributes from a coarse file

    } else if (action == FILECOARSE) {
      
      if (iarg+5 > narg) error->all(FLERR,"Illegal custom command");

      if (mode == PARTICLE || mode == SURF)
        error->all(FLERR,"Custom command cannot use action file/coarse "
                   "with style particle or surf");

      // # of coarse files and filestyle

      int numfile = input->inumeric(FLERR,arg[iarg+1]);
      int filestyle;
      if (strcmp(arg[iarg+2],"text") == 0) filestyle = TEXT;
      else if (strcmp(arg[iarg+2],"binary") == 0) filestyle = BINARY;
      
      // file name
      // if numfile > 1, fname must have "%" wildcard char
      // if external, fname must have "*" wildcard char
      
      int n = strlen(arg[iarg+3]) + 1;
      char *fname = new char[n];
      strcpy(fname,arg[iarg+3]);

      if (numfile > 1 && strchr(fname,'%') == NULL)
        error->all(FLERR,"Fix custom filename must have % wildcard char");
      if (external && strchr(fname,'*') == NULL)
        error->all(FLERR,"Fix custom filename must have * wildcard char");

      // colcount = # of attribute values per file line
      
      int colcount = input->inumeric(FLERR,arg[iarg+4]);
      if (colcount < 1)
        error->all(FLERR,"Custom command file/coarse column count is invalid");
      if (iarg+5+colcount > narg) error->all(FLERR,"Illegal custom command");

      // create vectors of attribute name settings

      int *cindex = new int[colcount];
      int *ctype = new int[colcount];
      int *csize = new int[colcount];
      int *ccol = new int[colcount];

      for (int i = 0; i < colcount; i++) {
        int n = strlen(arg[iarg+5+i]) + 1;  
        char *aname = new char[n];
        strcpy(aname,arg[iarg+5+i]);
        ccol[i] = attribute_bracket(aname);
        cindex[i] = grid->find_custom(aname);
        if (cindex[i] < 0)
          error->all(FLERR,"Custom attribute name does not exist");
        ctype[i] = grid->etype[cindex[i]];
        csize[i] = grid->esize[cindex[i]];
        if (csize[i] && ccol[i] == 0)
          error->all(FLERR,"Custom attribute array requires bracketed index");
        if (csize[i] == 0 && ccol[i])
          error->all(FLERR,"Custom attribute vector cannot use bracketed index");
        delete [] aname;
      }

      // if not external: invoke action now
      // else: store info in Action list for FixCustom
      
      // make a function to expand both wildcards with an integer value
      // use it everywhere

      // NOTE: make remainder of this a function callable from process_actions()

      if (!external) {
        if (comm->me == 0 && screen)
          fprintf(screen,"Reading custom file/coarse %s ...\n",fname);

        // read all coarse files, share results with all procs
        // sets ncoarse, xyz_coarse, values_coarse

        read_coarse_files(fname,external,numfile,colcount);

        if (comm->me == 0 && screen)
          fprintf(screen,"  %d coarse grid points with %d values each\n",ncoarse,colcount);

        // warn if > 1M points ?

        // form k-d tree for 2d or 3d points

        KDTree *kdtree = new KDTree(sparta,domain->dimension,ncoarse,xyz_coarse);

        int *list;
        memory->create(list,ncoarse,"custom:list");
        for (int i = 0; i < ncoarse; i++) list[i] = i;
        
        kdtree->create_tree(-1,ncoarse,list);

        memory->destroy(list);
        
        // print stats on created KD tree

        if (comm->me == 0) kdtree->stats_tree();
        
        // loop over my owned grid cells
        // find closest coarse point
        // find all coarse points within cutoff = closest distance + EPSILON
        // error if any of closest coarse points are duplicates
        // set grid cell attributes to average of closest coarse point values
        
        // NOTE: gather search stats internal to KDTree class, so can print them
        // NOTE: do other uses of grid-style vars set inside = 0.0 ?
        //       maybe should just use a far-away coarse point?
        //       or have user flag to choose zero or far-away

        Grid::ChildCell *cells = grid->cells;
        int nglocal = grid->nlocal;
        double ctr[3];
        double distsq;

        int *closest;
        double *closest_distsq;
        memory->create(closest,nglocal,"KDTree:closest");
        memory->create(closest_distsq,nglocal,"KDTree:closest_distsq");
        
        for (int i = 0; i < nglocal; i++) {
          if (cells[i].nsplit <= 0) continue;
          ctr[0] = 0.5 * (cells[i].lo[0] + cells[i].hi[0]);
          ctr[1] = 0.5 * (cells[i].lo[1] + cells[i].hi[1]);
          ctr[2] = 0.5 * (cells[i].lo[2] + cells[i].hi[2]);
          closest[i] = kdtree->find_nearest(ctr,0,closest_distsq[i]);
          printf("CLOSEST i %d icoarse %d distsq %g\n",i,closest[i],
                 sqrt(closest_distsq[i]));
        }

        int ncount;
        double cut,cutsq;
        int plist[MAXTIE];
        double dlist[MAXTIE];

        int count = 0;
        
        for (int i = 0; i < nglocal; i++) { 
          if (cells[i].nsplit <= 0) continue;
          ctr[0] = 0.5 * (cells[i].lo[0] + cells[i].hi[0]);
          ctr[1] = 0.5 * (cells[i].lo[1] + cells[i].hi[1]);
          ctr[2] = 0.5 * (cells[i].lo[2] + cells[i].hi[2]);
          cut = sqrt(closest_distsq[i]) * (1.0*EPSCUT);
          cutsq = cut*cut;
          ncount = 0;
          kdtree->find_within_cutoff(ctr,0,cutsq,ncount,plist,dlist);
          count += ncount;
        }

        // print stats on tree searches

        int count_all;
        MPI_Allreduce(&count,&count_all,1,MPI_INT,MPI_SUM,world);

        if (comm->me == 0) printf("NEIGHS FOUND %g\n",1.0*count_all/nglocal);
        
        // if (comm->me == 0) kdtree->stats_search();
        
        // use search results to set custom attributes of each owned grid cell

        // for mode = GRID
        //   set estatus of all changed custom vecs/arrays to 1
        //   b/c file read stores values for owned+ghost cells
        //   NOTE: should I still do this for ghost cells ?
        //         either here or in FILESTYLE
        //         because owned+ghost can be ALL grid cells
        
        for (int i = 0; i < colcount; i++)
          grid->estatus[cindex[i]] = 1;

        // clean up
        
        delete [] fname;
        delete [] cindex;
        delete [] ctype;
        delete [] csize;
        delete [] ccol;

        delete kdtree;
        memory->destroy(xyz_coarse);
        memory->destroy(values_coarse);
        memory->destroy(closest);
        memory->destroy(closest_distsq);
        
      } else {
        actions[naction].action = action;
        actions[naction].numfile = numfile;
        actions[naction].filestyle = filestyle;
        actions[naction].fname = fname;
        actions[naction].colcount = colcount;
        actions[naction].cindex_file = cindex;
        actions[naction].ctype_file = ctype;
        actions[naction].csize_file = csize;
        actions[naction].ccol_file = ccol;
        naction++;
      }
      
      iarg += 5 + colcount;
    }
  }

  return count;
}


/* ----------------------------------------------------------------------
   read one or more coarse files
   external = 0/1 for custom versus fix custom
   numfile = # of coarse grid files to read
   colcount = # of values per coarse grid point file must have
   commumicate grid info to all procs
   each proc stores copy of all coarse-grid info:
     ncoarse = # of coarse grid points
     xyz_coarse = coords of each coarse grid point
     values_coarse = values for each coarse grid point
---------------------------------------------------------------------- */

void Custom::read_coarse_files(char *fname, int external, int numfile, int colcount)
{
  // binary files not yet supported

  if (mode == BINARY) error->all(FLERR,"Custom file/coarse binary files not yet supported");

  // setup

  char *line = new char[MAXLINE];

  // if external:
  // filecurrent = replace '*' in fname with current timestep
  
  char *filecurrent;
        
  if (external) {
    filecurrent = new char[strlen(fname) + 16];
    char *ptr = strchr(fname,'*');
    *ptr = '\0';
    sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
            fname,update->ntimestep,ptr+1);
    *ptr = '*';
  } else filecurrent = fname;

  // storage for coarse grid points I read in
  
  int ncoarse_me = 0;
  double **xyz_coarse_me = NULL;
  double **values_coarse_me = NULL;
  
  for (int iproc = comm->me; iproc < numfile; iproc += comm->nprocs) {

    // if multiple files:
    // filewhich = replace '*' in filecurrent with current timestep 

    char *filewhich;
    
    if (numfile > 1) {
      filewhich = new char[strlen(filecurrent) + 16];
      char *ptr = strchr(fname,'%');
      *ptr = '\0';
      sprintf(filewhich,"%s%d%s",filecurrent,iproc+1,ptr+1);
      *ptr = '%';
    } else filewhich = filecurrent;
    
    // read a single coarse file = filewhich
  
    // NOTE: print name of bad file ?
    FILE *fp = fopen(filewhich,"r");
    if (fp == NULL) error->one(FLERR,"Could not open custom coarse file");
  
    // read past header portion of file
    // comments or blank lines are allowed
  
    char *eof,*ptr;
    
    while (1) {
      eof = fgets(line,MAXLINE,fp);
      if (eof == NULL) error->one(FLERR,"Unexpected end of custom coarse file");

      // trim anything from '#' onward
      // if line is blank, continue
      // else break and read nfile

      if ((ptr = strchr(line,'#'))) *ptr = '\0';
      if (strspn(line," \t\n\r") == strlen(line)) continue;
      break;
    }

    // line: Npoints Nvalues

    int npoints,nvalues;
    sscanf(line,"%d %d",&npoints,&nvalues);

    if (nvalues != colcount)
      error->one(FLERR,"Incorrect line format in custom coarse file");

    // reallocate local storage to store more grid point info

    memory->grow(xyz_coarse_me,ncoarse_me+npoints,domain->dimension,"custom:xyz_coarse_me");
    memory->grow(values_coarse_me,ncoarse_me+npoints,colcount,"custom:values_coarse_me");

    for (int i = ncoarse_me; i < ncoarse_me+npoints; i++) {
      eof = fgets(line,MAXLINE,fp);
      if (eof == NULL) error->one(FLERR,"Unexpected end of custom coarse file");

      if (i == 0) {
        int nwords = input->count_words(line);
        if (nwords != colcount + 1 + domain->dimension)
          error->one(FLERR,"Incorrect line format in custom coarse file");
      }

      ptr = strtok(line," \t\n\r\f");   // skip coarse cell ID
      xyz_coarse_me[i][0] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
      xyz_coarse_me[i][1] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
      if (domain->dimension == 3)
        xyz_coarse_me[i][2] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));

      for (int j = 0; j < colcount; j++)
        values_coarse_me[i][j] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
    }

    fclose(fp);
    if (numfile > 1) delete [] filewhich;
    ncoarse_me += npoints;
  }

  if (external) delete [] filecurrent;

  // check that each coarse point is inside or on simulation box

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  int count = 0;
  int flag;
  
  for (int i = 0; i < ncoarse_me; i++) {
    flag = 0;
    if (xyz_coarse_me[i][0] < boxlo[0] || xyz_coarse_me[i][0] > boxhi[0]) flag = 1;
    if (xyz_coarse_me[i][1] < boxlo[1] || xyz_coarse_me[i][1] > boxhi[1]) flag = 1;
    if (domain->dimension == 3)
      if (xyz_coarse_me[i][2] < boxlo[2] || xyz_coarse_me[i][2] > boxhi[2]) flag = 1;
    if (flag) count++;
  }

  int count_all;
  MPI_Allreduce(&count,&count_all,1,MPI_INT,MPI_SUM,world);

  if (count_all) {
    error->all(FLERR,"How many coarse grid points are outside simulation box");
  }

  // perform Allgatherv() so each proc has copy of all coarse points read by all procs
  // done once for xyz_coarse, another for values_coarse

  MPI_Allreduce(&ncoarse_me,&ncoarse,1,MPI_INT,MPI_SUM,world);
  memory->create(xyz_coarse,ncoarse,domain->dimension,"custom:xyz_coarse");
  memory->create(values_coarse,ncoarse,colcount,"custom:values_coarse");
  
  int nprocs = comm->nprocs;
  int *recvcounts,*displs;
  memory->create(recvcounts,nprocs,"custom:recvcounts");
  memory->create(displs,nprocs,"custom:displs");

  int nper = domain->dimension * ncoarse_me;

  MPI_Allgather(&nper,1,MPI_INT,recvcounts,1,MPI_INT,world);
  displs[0] = 0;
  for (int i = 1; i < nprocs; i++) displs[i] = displs[i-1] + recvcounts[i-1];

  if (ncoarse_me)
    MPI_Allgatherv(xyz_coarse_me[0],nper,MPI_DOUBLE,
                   xyz_coarse[0],recvcounts,displs,MPI_DOUBLE,world);
  else
    MPI_Allgatherv(NULL,nper,MPI_DOUBLE,
                   xyz_coarse[0],recvcounts,displs,MPI_DOUBLE,world);

  nper = colcount * ncoarse_me;

  MPI_Allgather(&nper,1,MPI_INT,recvcounts,1,MPI_INT,world);
  displs[0] = 0;
  for (int i = 1; i < nprocs; i++) displs[i] = displs[i-1] + recvcounts[i-1];

  if (ncoarse_me)
    MPI_Allgatherv(values_coarse_me[0],nper,MPI_DOUBLE,
                   values_coarse[0],recvcounts,displs,MPI_DOUBLE,world);
  else
    MPI_Allgatherv(NULL,nper,MPI_DOUBLE,
                   values_coarse[0],recvcounts,displs,MPI_DOUBLE,world);
    
  // clean up

  delete [] line;
  memory->destroy(recvcounts);
  memory->destroy(displs);
  memory->destroy(xyz_coarse_me);
  memory->destroy(values_coarse_me);
}

/* ----------------------------------------------------------------------
   invoke Action list of stored SET and FILESTYLE actions
   invoked by FixCustom
   return count of attribute values changed by this proc
---------------------------------------------------------------------- */

bigint Custom::process_actions()
{
  bigint count = 0;
  
  // process list of stored actions
  // only SET and FILESTYLE actions for

  for (int i = 0; i < naction; i++) {
      
   // set a custom vector or column of custom array via a variable

    if (actions[i].action == SET) {
      int vstyle = actions[i].vstyle;
      int vindex = actions[i].vindex;
      int cindex = actions[i].cindex;
      int ctype = actions[i].ctype;
      int csize = actions[i].csize;
      int ccol = actions[i].ccol;
      int groupbit = actions[i].groupbit;
      Mixture *mixture = actions[i].mixture;
      Region *region = actions[i].region;

      count += action_set(vstyle,vindex,cindex,ctype,csize,ccol,
                          groupbit,mixture,region);
      
    // set multiple custom attributes from a file

    } else if (actions[i].action == FILESTYLE) {

      char *fname = actions[i].fname;
      int colcount = actions[i].colcount;
      int *cindex = actions[i].cindex_file;
      int *ctype = actions[i].ctype_file;
      int *csize = actions[i].csize_file;
      int *ccol = actions[i].ccol_file;

      // replace '*' in fname with current timestep (logic from Dump class)
      // read the file and set attributes via input script column names

      char *filecurrent = new char[strlen(fname) + 16];
      char *ptr = strchr(fname,'*');
      *ptr = '\0';
      sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
              fname,update->ntimestep,ptr+1);
      *ptr = '*';

      count += read_file(mode,colcount,
                         cindex,ctype,csize,ccol,filecurrent);
      
      delete [] filecurrent;

      // for mode = GRID
      //   set estatus of all changed custom vecs/arrays to 1
      //   b/c file read stores values for owned+ghost cells
      // for mode = SURF
      //   set estatus of all changed custom vecs/arrays to 0
      
      if (mode == GRID)
        for (int i = 0; i < colcount; i++)
          grid->estatus[cindex[i]] = 1;
      else if (mode == SURF)
        for (int i = 0; i < colcount; i++)
          surf->estatus[cindex[i]] = 0;

    } else if (actions[i].action == FILECOARSE) {

      char *fname = actions[i].fname;
      int colcount = actions[i].colcount;
      int *cindex = actions[i].cindex_file;
      int *ctype = actions[i].ctype_file;
      int *csize = actions[i].csize_file;
      int *ccol = actions[i].ccol_file;





      
      // for mode = GRID
      //   set estatus of all changed custom vecs/arrays to 1
      //   b/c file read stores values for owned+ghost cells
      
      if (mode == GRID)
        for (int i = 0; i < colcount; i++)
          grid->estatus[cindex[i]] = 1;
    }
  }

  return count;
}

/* ----------------------------------------------------------------------
   perform set action for mode = PARTICLE, GRID, SURF
   set a custom vector or column of custom array via a variable
   return count of attribute values changed by this proc
---------------------------------------------------------------------- */

bigint Custom::action_set(int vstyle, int vindex,
                          int cindex, int ctype, int csize, int ccol,
                          int groupbit, Mixture *mixture, Region *region)
{
  bigint count = 0;
  
  double scalar = 0.0;
  double *vector = NULL;

  Variable *variable = input->variable;

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

  if (mode == PARTICLE)
    count = set_particle(mixture,region,
                         cindex,ctype,csize,ccol,scalar,vector);
  else if (mode == GRID)
    count = set_grid(groupbit,region,
                     cindex,ctype,csize,ccol,scalar,vector);
  else if (mode == SURF)
    count = set_surf(groupbit,region,
                     cindex,ctype,csize,ccol,scalar,vector);
  
  memory->destroy(vector);

  // for mode = GRID
  //   set estatus of custom vec/array to 0
  //   b/c variable evalulation only sets values for owned cells
  // for mode = SURF
  //   set estatus of all changed custom vecs/arrays to 0

  if (mode == GRID) grid->estatus[cindex] = 0;
  else if (mode == SURF) surf->estatus[cindex] = 0;

  return count;
}    

/* ----------------------------------------------------------------------
   set a PARTICLE custom vector or column of custom array via a variable
   scalar/vector = evaulated variable result
   return count of attribute values changed by this proc
---------------------------------------------------------------------- */

bigint Custom::set_particle(Mixture *mixture, Region *region,
                            int cindex, int ctype, int csize, int ccol,
                            double scalar, double *vector)
{
  Particle::OnePart *particles = particle->particles;
  int *species2species = mixture->species2species;
  int nlocal = particle->nlocal;

  int *choose;
  memory->create(choose,nlocal,"set:choose");
  memset(choose,0,nlocal*sizeof(int));

  int flag;

  bigint count = 0;
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
   set a GRID custom vector or column of custom array via a variable
   scalar/vector = evaulated variable result
   return count of attribute values changed by this proc
---------------------------------------------------------------------- */

bigint Custom::set_grid(int groupbit, Region *region,
                        int cindex, int ctype, int csize, int ccol,
                        double scalar, double *vector)
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  // select which cells to set, skip sub cells
  // region criterion is based on cell center point
  
  int *choose;
  memory->create(choose,nglocal,"set:choose");
  memset(choose,0,nglocal*sizeof(int));

  int flag;
  double point[3];

  bigint count = 0;
  for (int i = 0; i < nglocal; i++) {
    flag = 1;
    if (cells[i].nsplit <= 0) flag = 0;
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

/* ----------------------------------------------------------------------
   set a SURF custom vector or column of custom array via a variable
   scalar/vector = evaulated variable result
   return count of attribute values changed by this proc
---------------------------------------------------------------------- */

bigint Custom::set_surf(int groupbit, Region *region,
                        int cindex, int ctype, int csize, int ccol,
                        double scalar, double *vector)
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

  bigint count = 0;
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

/* ----------------------------------------------------------------------
   read a custom attribute file for mode = GRID or SURF
   assign values to custom grid or surf vectors/arrays
   return count of attributes assigned by this proc
---------------------------------------------------------------------- */

bigint Custom::read_file(int mode, int colcount,
                         int *cindex, int *ctype, int *csize, int *ccol,
                         char *filename)
{
  // setup read buffers
  // NOTE: should these be created by Memory class ?
  
  char *line = new char[MAXLINE];
  char *buffer = new char[CHUNK*MAXLINE];

  // set ivec,dvec,iarray,darray pointers
  // only one will be active for each value in input line
  
  int **ivec = new int*[colcount];
  double **dvec = new double*[colcount];
  int ***iarray = new int**[colcount];
  double ***darray = new double**[colcount];

  for (int j = 0; j < colcount; j++) {
    if (ctype[j] == INT) {
      if (csize[j] == 0) {
        if (mode == GRID)
          ivec[j] = grid->eivec[grid->ewhich[cindex[j]]];
        else if (mode == SURF)
          ivec[j] = surf->eivec[surf->ewhich[cindex[j]]];
      } else {
        if (mode == GRID)
          iarray[j] = grid->eiarray[grid->ewhich[cindex[j]]];
        else if (mode == SURF)
          iarray[j] = surf->eiarray[surf->ewhich[cindex[j]]];
      }
    } else if (ctype[j] == DOUBLE) {
      if (csize[j] == 0) {
        if (mode == GRID)
          dvec[j] = grid->edvec[grid->ewhich[cindex[j]]];
        else if (mode == SURF)
          dvec[j] = surf->edvec[surf->ewhich[cindex[j]]];
      } else {
        if (mode == GRID)
          darray[j] = grid->edarray[grid->ewhich[cindex[j]]];
        else if (mode == SURF)
          darray[j] = surf->edarray[surf->ewhich[cindex[j]]];
      }
    }
  }
  
  // ensure grid cell IDs are hashed
  
  Grid::MyHash *hash;

  if (mode == GRID) {
    if (!grid->hashfilled) grid->rehash();
    hash = grid->hash;
  }

  // nsurf = max ID of a surf in file
  
  bigint nsurf;
  if (mode == SURF) nsurf = surf->nsurf;
  
  // read file

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  // NOTE: support compressed files like read_grid ?

  int me = comm->me;
  int nprocs = comm->nprocs;
  FILE *fp;

  if (me == 0) {
    fp = fopen(filename,"r");
    if (fp == NULL) error->one(FLERR,"Could not open custom attribute file");
  }
  
  // read header portion of file
  // comments or blank lines are allowed
  // nfile = count of attribute lines in file
  // NOTE: allow for nfile to be a bigint ?
  
  int nfile,nvalues;
  
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

    // line: Nfile Nvalues

    sscanf(line,"%d %d",&nfile,&nvalues);
    //sscanf(line,BIGINT_FORMAT,&nfile);
  }

  MPI_Bcast(&nfile,1,MPI_INT,0,world);
  MPI_Bcast(&nvalues,1,MPI_INT,0,world);

  if (nvalues != colcount)
    error->all(FLERR,"Incorrect line format in custom attribute file");

  // read and broadcast one CHUNK of lines at a time

  bigint count = 0;
  bigint nread = 0;
  bigint fcount = 0;

  int i,m,nchunk,index,iproc;
  char *next,*buf,*idptr;
  cellint id;
  
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

    // process nchunk lines and assign attribute values if I store grid/surf ID
    // store means as owned or ghost grid cell
    
    buf = buffer;

    for (i = 0; i < nchunk; i++) {
      next = strchr(buf,'\n');
      if (next == NULL) printf("NULL proc %d i %d\n",me,i);
      *next = '\0';
      int nwords = input->count_words(buf);
      *next = '\n';

      if (nwords != colcount + 1)
	error->all(FLERR,"Incorrect line format in custom attribute file");

      // grid ID will match either an owned or ghost grid cell
      
      if (mode == GRID) {
        idptr = strtok(buf," \t\n\r\f");
        id = ATOCELLINT(idptr);
        if (id <= 0) error->all(FLERR,"Invalid cell ID in custom attribute grid file");

        if (hash->find(id) == hash->end()) {
          buf = next + 1;
          continue;
        }
        index = (*hash)[id];

      // surf ID will only match for the owning proc

      } else if (mode == SURF) {
        idptr = strtok(buf," \t\n\r\f");
        id = ATOSURFINT(idptr);
        if (id <= 0 || id > nsurf)
          error->all(FLERR,"Invalid surf ID in custom attribute surf file");

        iproc = (id-1) % nprocs;
        if (iproc != me) {
          buf = next + 1;
          continue;
        }
        index = (id-1) / nprocs;
      }

      // assign all attribute values for this grid cell or surf

      for (int j = 0; j < colcount; j++) {
        if (ctype[j] == INT) {
          if (csize[j] == 0)
            ivec[j][index] = input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
          else
            iarray[j][index][ccol[j]-1] = input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
        } else if (ctype[j] == DOUBLE) {
          if (csize[j] == 0)
            dvec[j][index] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
          else
            darray[j][index][ccol[j]-1] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
        }
      }

      count += colcount;
      buf = next + 1;
    }

    // increment nread and fcount and continue to next chunk
    
    nread += nchunk;
    fcount++;
  }
  
  // close file

  if (me == 0) {
    //if (compressed) pclose(fp);
    //else fclose(fp);
    fclose(fp);
  }
  
  // free read buffers and vec/array ptrs
  
  delete [] line;
  delete [] buffer;
  delete [] ivec;
  delete [] dvec;
  delete [] iarray;
  delete [] darray;

  return count;
}

/* ---------------------------------------------------------------------- 
  process an attribute name with optional bracketed index name[N]
   ccol = 0 if no brackets (vector attribute)
   ccol = N if bracktes (Nth column of array attribute)
   return ccol
---------------------------------------------------------------------- */

int Custom::attribute_bracket(char *aname)
{
  int ccol;
  char *ptr = strchr(aname,'[');
  if (ptr) {
    if (aname[strlen(aname)-1] != ']')
      error->all(FLERR,"Custom command attribute name is invalid");
    ccol = atoi(ptr+1);
    *ptr = '\0';
  } else ccol = 0;

  return ccol;
} 

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// KDTree class
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

enum{STUB,BRANCH,LEAF};
#define DELTANODE 4       // NOTE: make this larger after debugging

/* ---------------------------------------------------------------------- */

KDTree::KDTree(SPARTA *sparta, int dimension, int n, double **coords) : Pointers(sparta)
{
  dim = dimension;
  points = coords;
  npoints = n;
  
  tree = NULL;
  ntree = maxtree = 0;
}

/* ---------------------------------------------------------------------- */

KDTree::~KDTree()
{
  memory->sfree(tree);

  // NOTE: check if all memory is cleaned up
}

/* ----------------------------------------------------------------------
   recursive build of KD tree
   n = # of points
   plist = list of N indices into coords array
---------------------------------------------------------------------- */

void KDTree::create_tree(int iparent, int n, int *plist)
{
  // zero points, create STUB node

  if (n == 0) {
    if (ntree == maxtree) {
      maxtree += DELTANODE;
      tree = (Node *) memory->srealloc(tree,maxtree*sizeof(Node),"KDTree:tree");
    }
    tree[ntree].which = STUB;
    tree[ntree].iparent = iparent;
    ntree++;
    return;
  }

  // single point, create LEAF node
  
  if (n == 1) {
    if (ntree == maxtree) {
      maxtree += DELTANODE;
      tree = (Node *) memory->srealloc(tree,maxtree*sizeof(Node),"KDTree:tree");
    }
    tree[ntree].which = LEAF;
    tree[ntree].iparent = iparent;
    tree[ntree].ipoint = plist[0];
    ntree++;
    return;
  }

  // 2 or more points, create BRANCH node, continue recursing

  if (ntree == maxtree) {
    maxtree += DELTANODE;
    tree = (Node *) memory->srealloc(tree,maxtree*sizeof(Node),"KDTree:tree");
  }
  tree[ntree].which = BRANCH;
  tree[ntree].iparent = iparent;

  // calculate bbox around list of points
  
  double bboxlo[3] = {BIG,BIG,BIG};
  double bboxhi[3] = {-BIG,-BIG,-BIG};

  for (int i = 0; i < n; i++) {
    bboxlo[0] = MIN(bboxlo[0],points[plist[i]][0]);
    bboxhi[0] = MAX(bboxhi[0],points[plist[i]][0]);
    bboxlo[1] = MIN(bboxlo[1],points[plist[i]][1]);
    bboxhi[1] = MAX(bboxhi[1],points[plist[i]][1]);
    if (dim == 3) {
      bboxlo[2] = MIN(bboxlo[2],points[plist[i]][2]);
      bboxhi[2] = MAX(bboxhi[2],points[plist[i]][2]);
    }
  }

  /* NOTE: no need to store bbox in Tree ?
  tree[ntree].bboxlo[0] = bboxlo[0];
  tree[ntree].bboxlo[1] = bboxlo[1];
  tree[ntree].bboxlo[2] = bboxlo[2];
  tree[ntree].bboxhi[0] = bboxhi[0];
  tree[ntree].bboxhi[1] = bboxhi[1];
  tree[ntree].bboxhi[2] = bboxhi[2];
  */
  
  // splitdim = which dim to split points in (minimum bbox edge)
  // split = splitting value in that dim
  
  double xdelta = bboxhi[0] - bboxlo[0];
  double ydelta = bboxhi[1] - bboxlo[1];
  double zdelta = bboxhi[2] - bboxlo[2];

  int splitdim;
  if (xdelta >= ydelta) splitdim = 0;
  else splitdim = 1;
  if (dim == 3) {
    if (splitdim == 0 && zdelta > xdelta) splitdim = 2;
    if (splitdim == 1 && zdelta > ydelta) splitdim = 2;
  }
      
  double split = 0.5 * (bboxlo[splitdim] + bboxhi[splitdim]);
  
  tree[ntree].splitdim = splitdim;
  tree[ntree].split = split;

  // split plist into left and right lists

  int *leftlist,*rightlist;
  memory->create(leftlist,n,"KDTree:leftlist");
  memory->create(rightlist,n,"KDTree:rightlist");

  int nleft = 0;
  int nright = 0;
  int maxleft = 0;
  int maxright = 0;
  
  for (int i = 0; i < n; i++) {
    if (points[plist[i]][splitdim] <= split) leftlist[nleft++] = plist[i];
    else rightlist[nright++] = plist[i];
  }

  // recurse on left and right lists

  ntree++;
  int new_parent = ntree-1;

  /*
  printf("INITIAL TREE %d\n",itree);

  printf("CREATE LEFT: itree %d nleft %d leftlist %p %d\n",itree,nleft,leftlist,leftlist[0]);
  printf("   ");
  for (int i = 0; i < nleft; i++) printf(" %d",leftlist[i]);
  printf("\n");
  */
  
  tree[new_parent].left = ntree;
  create_tree(new_parent,nleft,leftlist);
  //printf("DESTROY LEFT: itree %d treeleft %d leftlist %p\n",itree,tree[itree].left,leftlist);
  memory->destroy(leftlist);

  /*
  printf("CREATE RIGHT: itree %d nright %d rightlist %p %d\n",itree,nright,rightlist,rightlist[0]);
  printf("   ");
  for (int i = 0; i < nright; i++) printf(" %d",rightlist[i]);
  printf("\n");
  */
  
  //printf("ITREE on RIGHT %d\n",itree);

  
  tree[new_parent].right = ntree;
  create_tree(new_parent,nright,rightlist);

  /*
  printf("DESTROY RIGHT: itree %d treeright %d rightlist %p\n",itree,tree[itree].right,rightlist);
  */
  memory->destroy(rightlist);

  //printf("FINAL TREE %d\n",itree);
  
  //return itree;
}

/* ----------------------------------------------------------------------
   recursive search of KD tree for point nearest to X
   find nearest coarse point to point X, starting at Node inode of KD tree
   returns index of nearest leaf Node and its squared distance (distsq) to point X
   algorithm:
   initial walk to leaf from inode to find distance to leaf node's point
     this can miss closer points in other branches
     to find all of them, walk back up to starting node
       each time examine opposite branch
       if it could possibly have a closer point, then walk down it
       do this recursively, so walk down untaken branches of every node visited
     efficiency comes form criteria for "could possibly have"
       if distance of point to split plane > current closest point
       then no need to walk that branch
---------------------------------------------------------------------- */

int KDTree::find_nearest(double *x, int inode, double &distsq)
{
  int startnode = inode;
  inode = walk_to_leaf(inode,x);
  int ipoint = tree[inode].ipoint;

  double dx,dy,dz;
  dx = x[0] - points[ipoint][0];
  dy = x[1] - points[ipoint][1];
  if (dim == 3) dz = x[2] - points[ipoint][2];
  else dz = 0.0;
  distsq = dx*dx + dy*dy + dz*dz;

  int ipoint_new;
  double distsq_new;
  
  while (inode != startnode) {
    int ichild = inode;
    inode = tree[inode].iparent;
    int splitdim = tree[inode].splitdim;
    double split = tree[inode].split;
    double delta = x[splitdim] - split;
    if (delta*delta <= distsq) {
      if (tree[inode].left == ichild)
        ipoint_new = find_nearest(x,tree[inode].right,distsq_new);
      else
        ipoint_new = find_nearest(x,tree[inode].left,distsq_new);
      if (distsq_new < distsq) {
        ipoint = ipoint_new;
        distsq = distsq_new;
      }
    }
  }

  return ipoint;
}

/* ----------------------------------------------------------------------
   recursive search of KD tree for all points within cutoff distance of X
   start search starting at Node inode of KD tree
   returns count of neighbors and plist,dlist
     plist = list of coarse point indices
     dlist = list of distances to coarse points
   algorithm:
   initial walk to leaf from inode to find distance to leaf node's point
     add it to list if wihin cutoff
     this can miss neighbor points in other branches
     to find all of them, walk back up to starting node
       each time examine opposite branch
       if it could possibly have a neighbor point, then walk down it
       do this recursively, so walk down untaken branches of every node visited
     efficiency comes form criteria for "could possibly have"
       if distance of point to split plane > cutoff
       then no need to walk that branch
---------------------------------------------------------------------- */

void KDTree::find_within_cutoff(double *x, int inode, double cutsq,
                                int &ncount, int *plist, double *dlist)
{
  int startnode = inode;
  inode = walk_to_leaf(inode,x);
  int ipoint = tree[inode].ipoint;

  double dx,dy,dz;
  dx = x[0] - points[ipoint][0];
  dy = x[1] - points[ipoint][1];
  if (dim == 3) dz = x[2] - points[ipoint][2];
  else dz = 0.0;
  double distsq = dx*dx + dy*dy + dz*dz;

  if (distsq < cutsq) {
    if (ncount == MAXTIE)
      error->one(FLERR,"KDTree cutoff induces too may ties");
    plist[ncount] = ipoint;
    dlist[ncount] = sqrt(distsq);
    ncount++;
  }

  while (inode != startnode) {
    int ichild = inode;
    inode = tree[inode].iparent;
    int splitdim = tree[inode].splitdim;
    double split = tree[inode].split;
    double delta = x[splitdim] - split;
    if (delta*delta <= cutsq) {
      if (tree[inode].left == ichild)
        find_within_cutoff(x,tree[inode].right,cutsq,ncount,plist,dlist);
      else
        find_within_cutoff(x,tree[inode].left,cutsq,ncount,plist,dlist);
    }
  }
}

/* ---------------------------------------------------------------------- */

int KDTree::walk_to_leaf(int inode, double *x)
{
  if (tree[inode].which == LEAF) return inode;
  
  int splitdim = tree[inode].splitdim;
  double split = tree[inode].split;
  if (x[splitdim] <= split) inode = walk_to_leaf(tree[inode].left,x);
  else inode = walk_to_leaf(tree[inode].right,x);
  return inode;
}

/* ---------------------------------------------------------------------- */

void KDTree::stats_tree()
{
  int nbranch = 0;
  int nleaf = 0;
  for (int i = 0; i < ntree; i++) {
    if (tree[i].which == BRANCH) nbranch++;
    if (tree[i].which == LEAF) nleaf++;
    if (tree[i].which == BRANCH) {
      printf("BRANCH NODE %d iparent %d which %d left %d right %d splitdim %d split %g\n",
             i,tree[i].iparent,tree[i].which,tree[i].left,tree[i].right,tree[i].splitdim,tree[i].split);
    }
    if (tree[i].which == LEAF) {
      printf("LEAF NODE %d iparent %d which %d ipoint %d\n",
             i,tree[i].iparent,tree[i].which,tree[i].ipoint);
    }
  }
  
  //int maxdepth = depthwalk(0,0,0);

  printf("KD tree stats:\n");
  printf("  %d = points stored\n",npoints);
  printf("  %d = # of nodes\n",ntree);
  printf("  %d = # of branches\n",nbranch);
  printf("  %d = # of leaves\n",nleaf);
  //printf("  %d = maxdepth\n",maxdepth);
}

/* ---------------------------------------------------------------------- */

int KDTree::depthwalk(int inode, int depth, int maxdepth)
{
  if (tree[inode].which == LEAF) return depth+1;

  //printf("DWALK inode %d %d %d\n",inode,depth,maxdepth);
  
  //int maxdepth = caller_maxdepth;
  maxdepth = MAX(maxdepth,depthwalk(tree[inode].left,depth+1,maxdepth));
  maxdepth = MAX(maxdepth,depthwalk(tree[inode].right,depth+1,maxdepth));
  return maxdepth;
}

/* ---------------------------------------------------------------------- */

void KDTree::stats_search()
{
}
