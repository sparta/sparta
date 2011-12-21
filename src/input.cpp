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
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "ctype.h"
#include "unistd.h"
#include "sys/stat.h"
#include "input.h"
#include "style_command.h"
#include "style_collide.h"
#include "universe.h"
#include "variable.h"
#include "domain.h"
#include "grid.h"
#include "particle.h"
#include "update.h"
#include "collide.h"
#include "output.h"
#include "stats.h"
#include "math_extra.h"
#include "error.h"
#include "memory.h"

#ifdef _OPENMP
#include "omp.h"
#endif

using namespace DSMC_NS;

#define MAXLINE 2048
#define DELTA 4

/* ---------------------------------------------------------------------- */

Input::Input(DSMC *dsmc, int argc, char **argv) : Pointers(dsmc)
{
  MPI_Comm_rank(world,&me);

  line = new char[MAXLINE];
  copy = new char[MAXLINE];
  work = new char[MAXLINE];
  narg = maxarg = 0;
  arg = NULL;

  echo_screen = 0;
  echo_log = 1;

  label_active = 0;
  labelstr = NULL;
  jump_skip = 0;

  if (me == 0) {
    nfile = maxfile = 1;
    infiles = (FILE **) memory->smalloc(sizeof(FILE *),"input:infiles");
    infiles[0] = infile;
  } else infiles = NULL;

  variable = new Variable(dsmc);

  // process command-line args
  // check for args "-var" and "-echo"
  // caller has already checked that sufficient arguments exist

  int iarg = 0;
  while (iarg < argc) {
    if (strcmp(argv[iarg],"-var") == 0 || strcmp(argv[iarg],"-v") == 0) {
      int jarg = iarg+2;
      while (jarg < argc && argv[jarg][0] != '-') jarg++;
      variable->set(argv[iarg+1],jarg-iarg-2,&argv[iarg+2]);
      iarg = jarg;
    } else if (strcmp(argv[iarg],"-echo") == 0 || 
	       strcmp(argv[iarg],"-e") == 0) {
      narg = 1;
      char **tmp = arg;        // trick echo() into using argv instead of arg
      arg = &argv[iarg+1];
      echo();
      arg = tmp;
      iarg += 2;
    } else iarg++;
  }
}

/* ---------------------------------------------------------------------- */

Input::~Input()
{
  // don't free command and arg strings
  // they just point to other allocated memory

  delete variable;
  delete [] line;
  delete [] copy;
  delete [] work;
  if (labelstr) delete [] labelstr;
  memory->sfree(arg);
  memory->sfree(infiles);
}

/* ----------------------------------------------------------------------
   process all input from infile
   infile = stdin or file if command-line arg "-in" was used
------------------------------------------------------------------------- */

void Input::file()
{
  int m,n;

  while (1) {
    
    // read a line from input script
    // if line ends in continuation char '&', concatenate next line(s)
    // n = length of line including str terminator, 0 if end of file
    // m = position of last printable char in line or -1 if blank line

    if (me == 0) {
      m = 0;
      while (1) {
	if (fgets(&line[m],MAXLINE-m,infile) == NULL) n = 0;
	else n = strlen(line) + 1;
	if (n == 0) break;
	m = n-2;
	while (m >= 0 && isspace(line[m])) m--;
	if (m < 0 || line[m] != '&') break;
      }
    }

    // bcast the line
    // if n = 0, end-of-file
    // error if label_active is set, since label wasn't encountered
    // if original input file, code is done
    // else go back to previous input file

    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n == 0) {
      if (label_active) error->all(FLERR,"Label wasn't found in input script");
      if (me == 0) {
	if (infile != stdin) fclose(infile);
	nfile--;
      }
      MPI_Bcast(&nfile,1,MPI_INT,0,world);
      if (nfile == 0) break;
      if (me == 0) infile = infiles[nfile-1];
      continue;
    }

    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // if n = MAXLINE, line is too long

    if (n == MAXLINE) {
      char str[MAXLINE+32];
      sprintf(str,"Input line too long: %s",line);
      error->all(FLERR,str);
    }

    // echo the command unless scanning for label

    if (me == 0 && label_active == 0) {
      if (echo_screen && screen) fprintf(screen,"%s",line); 
      if (echo_log && logfile) fprintf(logfile,"%s",line);
    }

    // parse the line
    // if no command, skip to next line in input script

    parse();
    if (command == NULL) continue;

    // if scanning for label, skip command unless it's a label command

    if (label_active && strcmp(command,"label") != 0) continue;

    // execute the command

    if (execute_command()) {
      char str[MAXLINE];
      sprintf(str,"Unknown command: %s",line);
      error->all(FLERR,str);
    }
  }
}

/* ----------------------------------------------------------------------
   process all input from filename
------------------------------------------------------------------------- */

void Input::file(const char *filename)
{
  // error if another nested file still open
  // if single open file is not stdin, close it
  // open new filename and set infile, infiles[0]

  if (me == 0) {
    if (nfile > 1)
      error->one(FLERR,"Another input script is already being processed");
    if (infile != stdin) fclose(infile);
    infile = fopen(filename,"r");
    if (infile == NULL) {
      char str[128];
      sprintf(str,"Cannot open input script %s",filename);
      error->one(FLERR,str);
    }
    infiles[0] = infile;
  } else infile = NULL;

  file();
}

/* ----------------------------------------------------------------------
   parse the command in single and execute it
   return command name to caller
------------------------------------------------------------------------- */

char *Input::one(const char *single)
{
  strcpy(line,single);

  // echo the command unless scanning for label
  
  if (me == 0 && label_active == 0) {
    if (echo_screen && screen) fprintf(screen,"%s\n",line); 
    if (echo_log && logfile) fprintf(logfile,"%s\n",line);
  }

  // parse the line
  // if no command, just return NULL

  parse();
  if (command == NULL) return NULL;

  // if scanning for label, skip command unless it's a label command

  if (label_active && strcmp(command,"label") != 0) return NULL;

  // execute the command and return its name

  if (execute_command()) {
    char str[MAXLINE];
    sprintf(str,"Unknown command: %s",line);
    error->all(FLERR,str);
  }

  return command;
}

/* ----------------------------------------------------------------------
   parse copy of command line
   strip comment = all chars from # on
   replace all $ via variable substitution
   command = first word
   narg = # of args
   arg[] = individual args
   treat text between single/double quotes as one arg
------------------------------------------------------------------------- */

void Input::parse()
{
  // make a copy to work on

  strcpy(copy,line);

  // strip any # comment by resetting string terminator
  // do not strip # inside single/double quotes

  char quote = '\0';
  char *ptr = copy;
  while (*ptr) {
    if (*ptr == '#' && !quote) {
      *ptr = '\0';
      break;
    }
    if (*ptr == quote) quote = '\0';
    else if (*ptr == '"' || *ptr == '\'') quote = *ptr;
    ptr++;
  }

  // perform $ variable substitution (print changes)
  // except if searching for a label since earlier variable may not be defined

  if (!label_active) substitute(copy,1);

  // command = 1st arg

  command = strtok(copy," \t\n\r\f");
  if (command == NULL) return;

  // point arg[] at each subsequent arg
  // treat text between single/double quotes as one arg
  // insert string terminators in copy to delimit args

  quote = '\0';
  int iarg,argstart;

  narg = 0;
  while (1) {
    if (narg == maxarg) {
      maxarg += DELTA;
      arg = (char **) memory->srealloc(arg,maxarg*sizeof(char *),"input:arg");
    }
    arg[narg] = strtok(NULL," \t\n\r\f");
    if (!arg[narg]) break;
    if (!quote && (arg[narg][0] == '"' || arg[narg][0] == '\'')) {
      quote = arg[narg][0];
      argstart = narg;
      arg[narg] = &arg[narg][1];
    }
    if (quote && arg[narg][strlen(arg[narg])-1] == quote) {
      for (iarg = argstart; iarg < narg; iarg++)
	arg[iarg][strlen(arg[iarg])] = ' ';
      arg[narg][strlen(arg[narg])-1] = '\0';
      narg = argstart;
      quote = '\0';
    }
    narg++;
  }

  if (quote) error->all(FLERR,"Unbalanced quotes in input line");
}

/* ----------------------------------------------------------------------
   substitute for $ variables in str and return it
   str assumed to be long enough to hold expanded version
   print updated string if flag is set and not searching for label
------------------------------------------------------------------------- */

void Input::substitute(char *str, int flag)
{
  // use work[] as scratch space to expand str, then copy back to str
  // do not replace $ inside single/double quotes
  // var = pts at variable name, ended by NULL
  //   if $ is followed by '{', trailing '}' becomes NULL
  //   else $x becomes x followed by NULL
  // beyond = pts at text following variable

  char *var,*value,*beyond;
  char quote = '\0';
  char *ptr = str;

  while (*ptr) {
    if (*ptr == '$' && !quote) {
      if (*(ptr+1) == '{') {
	var = ptr+2;
	int i = 0;
	while (var[i] != '\0' && var[i] != '}') i++;
	if (var[i] == '\0') error->one(FLERR,"Invalid variable name");
	var[i] = '\0';
	beyond = ptr + strlen(var) + 3;
      } else {
	var = ptr;
	var[0] = var[1];
	var[1] = '\0';
	beyond = ptr + strlen(var) + 1;
      }
      value = variable->retrieve(var);
      if (value == NULL) error->one(FLERR,"Substitution for illegal variable");

      *ptr = '\0';
      strcpy(work,str);
      if (strlen(work)+strlen(value) >= MAXLINE)
	error->one(FLERR,"Input line too long after variable substitution");
      strcat(work,value);
      if (strlen(work)+strlen(beyond) >= MAXLINE)
	error->one(FLERR,"Input line too long after variable substitution");
      strcat(work,beyond);
      strcpy(str,work);
      ptr += strlen(value);
      if (flag && me == 0 && label_active == 0) {
	if (echo_screen && screen) fprintf(screen,"%s",str); 
	if (echo_log && logfile) fprintf(logfile,"%s",str);
      }
      continue;
    }
    if (*ptr == quote) quote = '\0';
    else if (*ptr == '"' || *ptr == '\'') quote = *ptr;
    ptr++;
  }
}

/* ----------------------------------------------------------------------
   process a single parsed command
   return 0 if successful, -1 if did not recognize command
------------------------------------------------------------------------- */

int Input::execute_command()
{
  int flag = 1;

  if (!strcmp(command,"clear")) clear();
  else if (!strcmp(command,"echo")) echo();
  else if (!strcmp(command,"if")) ifthenelse();
  else if (!strcmp(command,"include")) include();
  else if (!strcmp(command,"jump")) jump();
  else if (!strcmp(command,"label")) label();
  else if (!strcmp(command,"log")) log();
  else if (!strcmp(command,"next")) next_command();
  else if (!strcmp(command,"partition")) partition();
  else if (!strcmp(command,"print")) print();
  else if (!strcmp(command,"shell")) shell();
  else if (!strcmp(command,"variable")) variable_command();

  else if (!strcmp(command,"boundary")) boundary();
  else if (!strcmp(command,"collisions")) collisions();
  else if (!strcmp(command,"dimension")) dimension();
  else if (!strcmp(command,"global")) global();
  else if (!strcmp(command,"mixture")) mixture();
  else if (!strcmp(command,"species")) species();
  else if (!strcmp(command,"stats")) stats();
  else if (!strcmp(command,"stats_modify")) stats_modify();
  else if (!strcmp(command,"stats_style")) stats_style();
  else if (!strcmp(command,"timestep")) timestep();

  else flag = 0;

  // return if command was listed above

  if (flag) return 0;

  // check if command is added via style.h

  if (0) return 0;      // dummy line to enable else-if macro expansion

#define COMMAND_CLASS
#define CommandStyle(key,Class)         \
  else if (strcmp(command,#key) == 0) { \
    Class key(dsmc);			\
    key.command(narg,arg);              \
    return 0;                           \
  }
#include "style_command.h"
#undef CommandStyle
#undef COMMAND_CLASS

  // unrecognized command

  return -1;
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void Input::clear()
{
  if (narg > 0) error->all(FLERR,"Illegal clear command");
  dsmc->destroy();
  dsmc->create();
}

/* ---------------------------------------------------------------------- */

void Input::echo()
{
  if (narg != 1) error->all(FLERR,"Illegal echo command");

  if (strcmp(arg[0],"none") == 0) {
    echo_screen = 0;
    echo_log = 0;
  } else if (strcmp(arg[0],"screen") == 0) {
    echo_screen = 1;
    echo_log = 0;
  } else if (strcmp(arg[0],"log") == 0) {
    echo_screen = 0;
    echo_log = 1;
  } else if (strcmp(arg[0],"both") == 0) {
    echo_screen = 1;
    echo_log = 1;
  } else error->all(FLERR,"Illegal echo command");
}

/* ---------------------------------------------------------------------- */

void Input::ifthenelse()
{
  if (narg < 3) error->all(FLERR,"Illegal if command");

  // substitute for variables in Boolean expression for "if"
  // in case expression was enclosed in quotes
  // must substitute on copy of arg else will step on subsequent args

  char *scopy = new char[MAXLINE];
  strcpy(scopy,arg[0]);
  substitute(scopy,0);

  // evaluate Boolean expression for "if"

  double btest = variable->evaluate_boolean(scopy);

  // bound "then" commands

  if (strcmp(arg[1],"then") != 0) error->all(FLERR,"Illegal if command");

  int first = 2;
  int iarg = first;
  while (iarg < narg && 
	 (strcmp(arg[iarg],"elif") != 0 && strcmp(arg[iarg],"else") != 0))
    iarg++;
  int last = iarg-1;

  // execute "then" commands
  // make copies of all arg string commands
  // required because re-parsing a command via one() will wipe out args

  if (btest != 0.0) {
    int ncommands = last-first + 1;
    if (ncommands <= 0) error->all(FLERR,"Illegal if command");

    char **commands = new char*[ncommands];
    ncommands = 0;
    for (int i = first; i <= last; i++) {
      int n = strlen(arg[i]) + 1;
      if (n == 1) error->all(FLERR,"Illegal if command");
      commands[ncommands] = new char[n];
      strcpy(commands[ncommands],arg[i]);
      ncommands++;
    }
    
    for (int i = 0; i < ncommands; i++) input->one(commands[i]);
    
    for (int i = 0; i < ncommands; i++) delete [] commands[i];
    delete [] commands;
    delete [] scopy;

    return;
  }

  // done if no "elif" or "else"

  if (iarg == narg) {
    delete [] scopy;
    return;
  }

  // check "elif" or "else" until find commands to execute
  // substitute for variables and evaluate Boolean expression for "elif"
  // must substitute on copy of arg else will step on subsequent args
  // bound and execute "elif" or "else" commands

  while (1) {
    if (iarg+2 > narg) error->all(FLERR,"Illegal if command");
    if (strcmp(arg[iarg],"elif") == 0) {
      strcpy(scopy,arg[iarg+1]);
      substitute(scopy,0);
      btest = variable->evaluate_boolean(scopy);
      first = iarg+2;
    } else {
      btest = 1.0;
      first = iarg+1;
    }

    iarg = first;
    while (iarg < narg && 
	   (strcmp(arg[iarg],"elif") != 0 && strcmp(arg[iarg],"else") != 0))
      iarg++;
    last = iarg-1;

    if (btest == 0.0) continue;

    int ncommands = last-first + 1;
    if (ncommands <= 0) error->all(FLERR,"Illegal if command");

    char **commands = new char*[ncommands];
    ncommands = 0;
    for (int i = first; i <= last; i++) {
      int n = strlen(arg[i]) + 1;
      if (n == 1) error->all(FLERR,"Illegal if command");
      commands[ncommands] = new char[n];
      strcpy(commands[ncommands],arg[i]);
      ncommands++;
    }
    
    // execute the list of commands
    
    for (int i = 0; i < ncommands; i++) input->one(commands[i]);
    
    // clean up

    for (int i = 0; i < ncommands; i++) delete [] commands[i];
    delete [] commands;
    delete [] scopy;

    return;
  }
}

/* ---------------------------------------------------------------------- */

void Input::include()
{
  if (narg != 1) error->all(FLERR,"Illegal include command");

  if (me == 0) {
    if (nfile == maxfile) {
      maxfile++;
      infiles = (FILE **) 
        memory->srealloc(infiles,maxfile*sizeof(FILE *),"input:infiles");
    }
    infile = fopen(arg[0],"r");
    if (infile == NULL) {
      char str[128];
      sprintf(str,"Cannot open input script %s",arg[0]);
      error->one(FLERR,str);
    }
    infiles[nfile++] = infile;
  }
}

/* ---------------------------------------------------------------------- */

void Input::jump()
{
  if (narg < 1 || narg > 2) error->all(FLERR,"Illegal jump command");

  if (jump_skip) {
    jump_skip = 0;
    return;
  }

  if (me == 0) {
    if (strcmp(arg[0],"SELF") == 0) rewind(infile);
    else {
      if (infile != stdin) fclose(infile);
      infile = fopen(arg[0],"r");
      if (infile == NULL) {
	char str[128];
	sprintf(str,"Cannot open input script %s",arg[0]);
	error->one(FLERR,str);
      }
      infiles[nfile-1] = infile;
    }
  }

  if (narg == 2) {
    label_active = 1;
    if (labelstr) delete [] labelstr;
    int n = strlen(arg[1]) + 1;
    labelstr = new char[n];
    strcpy(labelstr,arg[1]);
  }
}

/* ---------------------------------------------------------------------- */

void Input::label()
{
  if (narg != 1) error->all(FLERR,"Illegal label command");
  if (label_active && strcmp(labelstr,arg[0]) == 0) label_active = 0;
}

/* ---------------------------------------------------------------------- */

void Input::log()
{
  if (narg != 1) error->all(FLERR,"Illegal log command");

  if (me == 0) {
    if (logfile) fclose(logfile);
    if (strcmp(arg[0],"none") == 0) logfile = NULL;
    else {
      logfile = fopen(arg[0],"w");
      if (logfile == NULL) {
	char str[128];
	sprintf(str,"Cannot open logfile %s",arg[0]);
	error->one(FLERR,str);
      }
    }
    if (universe->nworlds == 1) universe->ulogfile = logfile;
  }
}

/* ---------------------------------------------------------------------- */

void Input::next_command()
{
  if (variable->next(narg,arg)) jump_skip = 1;
}

/* ---------------------------------------------------------------------- */

void Input::partition()
{
  if (narg < 3) error->all(FLERR,"Illegal partition command");

  int yesflag;
  if (strcmp(arg[0],"yes") == 0) yesflag = 1;
  else if (strcmp(arg[0],"no") == 0) yesflag = 0;
  else error->all(FLERR,"Illegal partition command");

  int ilo,ihi;
  int flag = MathExtra::bounds(arg[1],universe->nworlds,ilo,ihi);
  if (flag) error->all(FLERR,"Partition numeric index is out of bounds");

  // copy original line to copy, since will use strtok() on it
  // ptr = start of 4th word

  strcpy(copy,line);
  copy[strlen(copy)-1] = '\0';
  char *ptr = strtok(copy," \t\n\r\f");
  ptr = strtok(NULL," \t\n\r\f");
  ptr = strtok(NULL," \t\n\r\f");
  ptr += strlen(ptr) + 1;
  ptr += strspn(ptr," \t\n\r\f");

  // execute the remaining command line on requested partitions

  if (yesflag) {
    if (universe->iworld+1 >= ilo && universe->iworld+1 <= ihi) one(ptr);
  } else {
    if (universe->iworld+1 < ilo || universe->iworld+1 > ihi) one(ptr);
  }
}

/* ---------------------------------------------------------------------- */

void Input::print()
{
  if (narg != 1) error->all(FLERR,"Illegal print command");

  // substitute for $ variables (no printing) and print arg

  substitute(arg[0],0);
  if (me == 0) {
    if (screen) fprintf(screen,"%s\n",arg[0]);
    if (logfile) fprintf(logfile,"%s\n",arg[0]);
  }
}

/* ---------------------------------------------------------------------- */

void Input::shell()
{
  if (narg < 1) error->all(FLERR,"Illegal shell command");

  if (strcmp(arg[0],"cd") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal shell command");
    chdir(arg[1]);

  } else if (strcmp(arg[0],"mkdir") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal shell command");
#if !defined(WINDOWS) && !defined(__MINGW32_VERSION) 
    if (me == 0)
      for (int i = 1; i < narg; i++)
	mkdir(arg[i], S_IRWXU | S_IRGRP | S_IXGRP);
#endif

  } else if (strcmp(arg[0],"mv") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal shell command");
    if (me == 0) rename(arg[1],arg[2]);

  } else if (strcmp(arg[0],"rm") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal shell command");
    if (me == 0)
      for (int i = 1; i < narg; i++) unlink(arg[i]);

  } else if (strcmp(arg[0],"rmdir") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal shell command");
    if (me == 0)
      for (int i = 1; i < narg; i++) rmdir(arg[i]);

  // use work to concat args back into one string separated by spaces
  // invoke string in shell via system()

  } else {
    strcpy(work,arg[0]);
    for (int i = 1; i < narg; i++) {
      strcat(work," ");
      strcat(work,arg[i]);
    }
    if (me == 0) system(work);
  }
}

/* ---------------------------------------------------------------------- */

void Input::variable_command()
{
  variable->set(narg,arg);
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   one function for each DSMC-specific input script command
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void Input::boundary()
{
  domain->set_boundary(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::collisions()
{
  if (narg < 1) error->all(FLERR,"Illegal collide command");

  if (collide) delete collide;

  if (strcmp(arg[0],"none") == 0) collide = NULL;

#define COLLIDE_CLASS
#define CollideStyle(key,Class) \
  else if (strcmp(arg[0],#key) == 0) \
    collide = new Class(dsmc,narg,arg);
#include "style_collide.h"
#undef CollideStyle
#undef COLLIDE_CLASS

  else error->all(FLERR,"Invalid collide style");
}

/* ---------------------------------------------------------------------- */

void Input::dimension()
{
  if (narg != 1) error->all(FLERR,"Illegal dimension command");
  if (domain->box_exist) 
    error->all(FLERR,"Dimension command after simulation box is defined");
  domain->dimension = atoi(arg[0]);
  if (domain->dimension != 2 && domain->dimension != 3)
    error->all(FLERR,"Illegal dimension command");
}

/* ---------------------------------------------------------------------- */

void Input::global()
{
  update->global(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::mixture()
{
  particle->add_mixture(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::species()
{
  particle->add_species(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::stats()
{
  if (narg != 1) error->all(FLERR,"Illegal stats command");
  int n = atoi(arg[0]);
  if (n < 0) error->all(FLERR,"Illegal stats command");
  output->stats_every = n;
}

/* ---------------------------------------------------------------------- */

void Input::stats_modify()
{
  output->stats->modify_params(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::stats_style()
{
  output->create_stats(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::timestep()
{
  if (narg != 1) error->all(FLERR,"Illegal timestep command");
  double dt = atof(arg[0]);
  if (dt <= 0.0) error->all(FLERR,"Illegal timestep command");
  update->dt = dt;
}
