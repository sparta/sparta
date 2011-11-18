/* ----------------------------------------------------------------------
   DSMC - Sandia parallel DSMC code
   www.sandia.gov/~sjplimp/dsmc.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2011) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level DSMC directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "string.h"
#include "dsmc.h"
#include "style_command.h"
#include "universe.h"
#include "input.h"
#include "particle.h"
#include "update.h"
#include "comm.h"
#include "domain.h"
#include "grid.h"
#include "surf.h"
#include "collide.h"
#include "output.h"
#include "timer.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;

/* ----------------------------------------------------------------------
   start up DSMC
   allocate fundamental classes (memory, error, universe, input)
   parse input switches
   initialize communicators, screen & logfile output
   input is allocated at end after MPI info is setup
------------------------------------------------------------------------- */

DSMC::DSMC(int narg, char **arg, MPI_Comm communicator)
{
  memory = new Memory(this);
  error = new Error(this);
  universe = new Universe(this,communicator);
  output = NULL;

  screen = NULL;
  logfile = NULL;

  // parse input switches

  int inflag = 0;
  int screenflag = 0;
  int logflag = 0;
  int partscreenflag = 0;
  int partlogflag = 0;
  int cudaflag = -1;
  int helpflag = 0;
  suffix = NULL;
  suffix_enable = 0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"-partition") == 0 || 
	strcmp(arg[iarg],"-p") == 0) {
      universe->existflag = 1;
      if (iarg+2 > narg) 
	error->universe_all(FLERR,"Invalid command-line argument");
      iarg++;
      while (iarg < narg && arg[iarg][0] != '-') {
	universe->add_world(arg[iarg]);
	iarg++;
      }
    } else if (strcmp(arg[iarg],"-in") == 0 || 
	       strcmp(arg[iarg],"-i") == 0) {
      if (iarg+2 > narg) 
	error->universe_all(FLERR,"Invalid command-line argument");
      inflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-screen") == 0 || 
	       strcmp(arg[iarg],"-sc") == 0) {
      if (iarg+2 > narg) 
	error->universe_all(FLERR,"Invalid command-line argument");
      screenflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-log") == 0 || 
	       strcmp(arg[iarg],"-l") == 0) {
      if (iarg+2 > narg) 
	error->universe_all(FLERR,"Invalid command-line argument");
      logflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-var") == 0 || 
	       strcmp(arg[iarg],"-v") == 0) {
      if (iarg+3 > narg)
	error->universe_all(FLERR,"Invalid command-line argument");
      iarg += 2;
      while (iarg < narg && arg[iarg][0] != '-') iarg++;
    } else if (strcmp(arg[iarg],"-echo") == 0 || 
	       strcmp(arg[iarg],"-e") == 0) {
      if (iarg+2 > narg)
	error->universe_all(FLERR,"Invalid command-line argument");
      iarg += 2;
    } else if (strcmp(arg[iarg],"-pscreen") == 0 || 
	       strcmp(arg[iarg],"-ps") == 0) {
      if (iarg+2 > narg) 
       error->universe_all(FLERR,"Invalid command-line argument");
      partscreenflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-plog") == 0 || 
	       strcmp(arg[iarg],"-pl") == 0) {
      if (iarg+2 > narg) 
       error->universe_all(FLERR,"Invalid command-line argument");
      partlogflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-cuda") == 0 || 
	       strcmp(arg[iarg],"-c") == 0) {
      if (iarg+2 > narg)
	error->universe_all(FLERR,"Invalid command-line argument");
      if (strcmp(arg[iarg+1],"on") == 0) cudaflag = 1;
      else if (strcmp(arg[iarg+1],"off") == 0) cudaflag = 0;
      else error->universe_all(FLERR,"Invalid command-line argument");
      iarg += 2;
    } else if (strcmp(arg[iarg],"-suffix") == 0 || 
	       strcmp(arg[iarg],"-sf") == 0) {
      if (iarg+2 > narg)
	error->universe_all(FLERR,"Invalid command-line argument");
      delete [] suffix;
      int n = strlen(arg[iarg+1]) + 1;
      suffix = new char[n];
      strcpy(suffix,arg[iarg+1]);
      suffix_enable = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-help") == 0 || 
	       strcmp(arg[iarg],"-h") == 0) {
      if (iarg+1 > narg)
	error->universe_all(FLERR,"Invalid command-line argument");
      helpflag = 1;
      iarg += 1;
    } else error->universe_all(FLERR,"Invalid command-line argument");
  }

  // if no partition command-line switch, universe is one world w/ all procs

  if (universe->existflag == 0) universe->add_world(NULL);

  // sum of procs in all worlds must equal total # of procs

  if (!universe->consistent())
    error->universe_all(FLERR,"Processor partitions are inconsistent");

  // universe cannot use stdin for input file

  if (universe->existflag && inflag == 0)
    error->universe_all(FLERR,"Must use -in switch with multiple partitions");

  // if no partition command-line switch, cannot use -pscreen option

  if (universe->existflag == 0 && partscreenflag)
    error->universe_all(FLERR,"Can only use -pscreen with multiple partitions");

  // if no partition command-line switch, cannot use -plog option

  if (universe->existflag == 0 && partlogflag)
    error->universe_all(FLERR,"Can only use -plog with multiple partitions");

  // set universe screen and logfile

  if (universe->me == 0) {
    if (screenflag == 0)
      universe->uscreen = stdout;
    else if (strcmp(arg[screenflag],"none") == 0)
      universe->uscreen = NULL;
    else {
      universe->uscreen = fopen(arg[screenflag],"w");
      if (universe->uscreen == NULL) 
	error->universe_one(FLERR,"Cannot open universe screen file");
    }
    if (logflag == 0) {
      universe->ulogfile = fopen("log.dsmc","w");
      if (universe->ulogfile == NULL) 
	error->universe_one(FLERR,"Cannot open log.dsmc");
    } else if (strcmp(arg[logflag],"none") == 0)
      universe->ulogfile = NULL;
    else {
      universe->ulogfile = fopen(arg[logflag],"w");
      if (universe->ulogfile == NULL) 
	error->universe_one(FLERR,"Cannot open universe log file");
    }
  }

  if (universe->me > 0) {
    if (screenflag == 0) universe->uscreen = stdout;
    else universe->uscreen = NULL;
    universe->ulogfile = NULL;
  }

  // universe does not exist on its own, only a single world
  // inherit settings from universe
  // set world screen, logfile, communicator, infile
  // open input script if from file

  if (universe->existflag == 0) {
    screen = universe->uscreen;
    logfile = universe->ulogfile;
    world = universe->uworld;
    infile = NULL;

    if (universe->me == 0) {
      if (inflag == 0) infile = stdin;
      else infile = fopen(arg[inflag],"r");
      if (infile == NULL) {
	char str[128];
	sprintf(str,"Cannot open input script %s",arg[inflag]);
	error->one(FLERR,str);
      }
    }

    if (universe->me == 0) {
      if (screen) fprintf(screen,"DSMC (%s)\n",universe->version);
      if (logfile) fprintf(logfile,"DSMC (%s)\n",universe->version);
    }

  // universe is one or more worlds
  // split into separate communicators
  // set world screen, logfile, communicator, infile
  // open input script

  } else {
    int me;
    MPI_Comm_split(universe->uworld,universe->iworld,0,&world);
    MPI_Comm_rank(world,&me);

    if (me == 0)
      if (partscreenflag == 0)
       if (screenflag == 0) {
         char str[32];
         sprintf(str,"screen.%d",universe->iworld);
         screen = fopen(str,"w");
         if (screen == NULL) error->one(FLERR,"Cannot open screen file");
       } else if (strcmp(arg[screenflag],"none") == 0)
         screen = NULL;
       else {
         char str[128];
         sprintf(str,"%s.%d",arg[screenflag],universe->iworld);
         screen = fopen(str,"w");
         if (screen == NULL) error->one(FLERR,"Cannot open screen file");
       }
      else if (strcmp(arg[partscreenflag],"none") == 0)
	screen = NULL;
      else {
	char str[128];
	sprintf(str,"%s.%d",arg[partscreenflag],universe->iworld);
	screen = fopen(str,"w");
	if (screen == NULL) error->one(FLERR,"Cannot open screen file");
      } else screen = NULL;
    
    if (me == 0)
      if (partlogflag == 0)
       if (logflag == 0) {
         char str[32];
         sprintf(str,"log.dsmc.%d",universe->iworld);
         logfile = fopen(str,"w");
         if (logfile == NULL) error->one(FLERR,"Cannot open logfile");
       } else if (strcmp(arg[logflag],"none") == 0)
         logfile = NULL;
       else {
         char str[128];
         sprintf(str,"%s.%d",arg[logflag],universe->iworld);
         logfile = fopen(str,"w");
         if (logfile == NULL) error->one(FLERR,"Cannot open logfile");
       }
      else if (strcmp(arg[partlogflag],"none") == 0)
	logfile = NULL;
      else {
	char str[128];
	sprintf(str,"%s.%d",arg[partlogflag],universe->iworld);
	logfile = fopen(str,"w");
	if (logfile == NULL) error->one(FLERR,"Cannot open logfile");
      } else logfile = NULL;
    
    if (me == 0) {
      infile = fopen(arg[inflag],"r");
      if (infile == NULL) {
	char str[128];
	sprintf(str,"Cannot open input script %s",arg[inflag]);
	error->one(FLERR,str);
      }
    } else infile = NULL;
    
    // screen and logfile messages for universe and world
    
    if (universe->me == 0) {
      if (universe->uscreen) {
	fprintf(universe->uscreen,"DSMC (%s)\n",universe->version);
	fprintf(universe->uscreen,"Running on %d partitions of processors\n",
		universe->nworlds);
      }
      if (universe->ulogfile) {
	fprintf(universe->ulogfile,"DSMC (%s)\n",universe->version);
	fprintf(universe->ulogfile,"Running on %d partitions of processors\n",
		universe->nworlds);
      }
    }
    
    if (me == 0) {
      if (screen) {
	fprintf(screen,"DSMC (%s)\n",universe->version);
	fprintf(screen,"Processor partition = %d\n",universe->iworld);
      }
      if (logfile) {
	fprintf(logfile,"DSMC (%s)\n",universe->version);
	fprintf(logfile,"Processor partition = %d\n",universe->iworld);
      }
    }
  }

  // check datatype settings in dsmctype.h

  if (sizeof(smallint) != sizeof(int))
    error->all(FLERR,"Smallint setting in dsmctype.h is invalid");
  if (sizeof(bigint) < sizeof(smallint))
    error->all(FLERR,"Bigint setting in dsmctype.h is invalid");

  int mpisize;
  MPI_Type_size(MPI_DSMC_BIGINT,&mpisize);
  if (mpisize != sizeof(bigint))
      error->all(FLERR,
		 "MPI_DSMC_BIGINT and bigint in dsmctype.h are not compatible");

  if (sizeof(smallint) != 4 || sizeof(bigint) != 8)
    error->all(FLERR,"Small,big integers are not sized correctly");

  // allocate input class now that MPI is fully setup

  input = new Input(this,narg,arg);

  // allocate fundamental classes

  create();

  // other top-level classes

  //ranmaster = new RanMars(this);

  // if helpflag set, print help and exit

  if (helpflag) {
    if (universe->me == 0) print_styles();
    error->done();
  }
}

/* ----------------------------------------------------------------------
   shutdown DSMC
   delete top-level classes
   delete fundamental classes
   close screen and log files in world and universe
   output files were already closed in destroy()
   delete fundamental classes
------------------------------------------------------------------------- */

DSMC::~DSMC()
{
  //delete ranmaster;
  
  destroy();

  if (universe->nworlds == 1) {
    if (logfile) fclose(logfile);
  } else {
    if (screen && screen != stdout) fclose(screen);
    if (logfile) fclose(logfile);
    if (universe->ulogfile) fclose(universe->ulogfile);
  }

  if (world != universe->uworld) MPI_Comm_free(&world);

  delete [] suffix;

  delete input;
  delete universe;
  delete error;
  delete memory;
}

/* ----------------------------------------------------------------------
   allocate single instance of top-level classes
   fundamental classes are allocated in constructor
   some classes have package variants
------------------------------------------------------------------------- */

void DSMC::create()
{
  particle = new Particle(this);
  update = new Update(this);
  comm = new Comm(this);
  domain = new Domain(this);
  grid = new Grid(this);
  surf = new Surf(this);
  collide = NULL;
  output = new Output(this);
  timer = new Timer(this);
}

/* ----------------------------------------------------------------------
   initialize top-level classes
------------------------------------------------------------------------- */

void DSMC::init()
{
  particle->init();
  update->init();
  comm->init();
  domain->init();
  grid->init();
  surf->init();
  if (collide) collide->init();
  output->init();
  timer->init();
}

/* ----------------------------------------------------------------------
   delete single instance of top-level classes
   fundamental classes are deleted in destructor
------------------------------------------------------------------------- */

void DSMC::destroy()
{
  delete particle;
  delete update;
  delete comm;
  delete domain;
  delete grid;
  delete surf;
  delete collide;
  delete output;
  delete timer;
}

/* ----------------------------------------------------------------------
   for each style, print name of all child classes build into executable
------------------------------------------------------------------------- */

void DSMC::print_styles()
{
  printf("\nList of style options included in this executable:\n\n");

  printf("Collide styles:");
#define COLLIDE_CLASS
#define CollideStyle(key,Class) printf(" %s",#key);
#include "style_collide.h"
#undef CollideStyle
#undef COLLIDE_CLASS
  printf("\n");

  printf("Command styles (add-on input script commands):");
#define COMMAND_CLASS
#define CommandStyle(key,Class) printf(" %s",#key);
#include "style_command.h"
#undef CommandStyle
#undef COMMAND_CLASS
  printf("\n");
}
