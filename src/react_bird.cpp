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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "react_bird.h"
#include "collide.h"
#include "update.h"
#include "particle.h"
#include "comm.h"
#include "random_park.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{DISSOCIATION,EXCHANGE,IONIZATION,RECOMBINATION};  // other react files
enum{ARRHENIUS,QUANTUM};                               // other react files

#define MAXREACTANT 2
#define MAXPRODUCT 3
#define MAXCOEFF 6

#define MAXLINE 1024
#define DELTALIST 16

/* ---------------------------------------------------------------------- */

ReactBird::ReactBird(SPARTA *sparta, int narg, char **arg) :
  React(sparta, narg, arg)
{
  if (narg != 2) error->all(FLERR,"Illegal react tce command");

  nlist = maxlist = 0;
  rlist = NULL;
  readfile(arg[1]);

  reactions = NULL;
  indices = NULL;
}

/* ---------------------------------------------------------------------- */

ReactBird::~ReactBird()
{
  for (int i = 0; i < maxlist; i++) {
    for (int j = 0; j < rlist[i].nreactant; j++)
      delete [] rlist[i].id_reactants[j];
    for (int j = 0; j < rlist[i].nproduct; j++)
      delete [] rlist[i].id_products[j];
    delete [] rlist[i].id_reactants;
    delete [] rlist[i].id_products;
    delete [] rlist[i].reactants;
    delete [] rlist[i].products;
    delete [] rlist[i].coeff;
  }
  memory->destroy(rlist);

  memory->destroy(reactions);
  memory->destroy(indices);
}

/* ---------------------------------------------------------------------- */

void ReactBird::init() 
{
  // convert species IDs to species indices
  // flag reactions as active/inactive depending on whether all species exist

  for (int m = 0; m < nlist; m++) {
    OneReaction *r = &rlist[m];
    r->active = 1;
    for (int i = 0; i < r->nreactant; i++) {
      r->reactants[i] = particle->find_species(r->id_reactants[i]);
      if (r->reactants[i] < 0) {
        r->active = 0;
        break;
      }
    }
    for (int i = 0; i < r->nproduct; i++) {
      r->products[i] = particle->find_species(r->id_products[i]);
      if (r->products[i] < 0) {
        r->active = 0;
        break;
      }
    }
  }

  // count possible reactions for each species pair
  // include J,I reactions in I,J list and vice versa

  memory->destroy(reactions);
  int nspecies = particle->nspecies;
  reactions = memory->create(reactions,nspecies,nspecies,
                             "react/tce:reactions");

  for (int i = 0; i < nspecies; i++)
    for (int j = 0; j < nspecies; j++)
      reactions[i][j].n = 0;

  int n = 0;
  for (int m = 0; m < nlist; m++) {
    OneReaction *r = &rlist[m];
    if (!r->active) continue;
    int i = r->reactants[0];
    int j = r->reactants[1];
    reactions[i][j].n++;
    n++;
    if (i == j) continue;
    reactions[j][i].n++;
    n++;
  }

  // allocate indices = entire list of reactions for all IJ pairs

  memory->destroy(indices);
  memory->create(indices,n,"react/tce:indices");

  // reactions[i][j].list = offset into full indices vector

  int offset = 0;
  for (int i = 0; i < nspecies; i++)
    for (int j = 0; j < nspecies; j++) {
      reactions[i][j].list = &indices[offset];
      offset += reactions[i][j].n;
    }

  // reactions[i][j].list = indices of possible reactions for each species pair
  // include J,I reactions in I,J list and vice versa

  for (int i = 0; i < nspecies; i++)
    for (int j = 0; j < nspecies; j++)
      reactions[i][j].n = 0;
 
  for (int m = 0; m < nlist; m++) {
    OneReaction *r = &rlist[m];
    if (!r->active) continue;
    int i = r->reactants[0];
    int j = r->reactants[1];
    reactions[i][j].list[reactions[i][j].n++] = m;
    if (i == j) continue;
    reactions[j][i].list[reactions[j][i].n++] = m;
  }

  // check for ionizatoin and recombitation reactions, not yet supported

  for (int i = 0; i < nspecies; i++)
    for (int j = 0; j < nspecies; j++) {
      int *list = reactions[i][j].list;
      for (int n = 0; n < reactions[i][j].n; n++)
        if (rlist[list[n]].type == IONIZATION || 
            rlist[list[n]].type == RECOMBINATION)
          error->all(FLERR,"Ionization and recombination reactions are "
                     "not yet implemented");
    }

  // modify Arrhenius coefficients for TCE model
  // C1,C2 Bird 94, p 127
  // initflag logic insures only done once per reaction

  Particle::Species *species = particle->species;

  for (int m = 0; m < nlist; m++) {
    OneReaction *r = &rlist[m];
    if (!r->active) continue;
    if (r->initflag) continue;
    r->initflag = 1;

    int isp = r->reactants[0];
    int jsp = r->reactants[1];

    double mass_i = species[isp].mass;
    double mass_j = species[jsp].mass;

    // symmetry parameter

    double epsilon = 1.0;         
    if (isp == jsp) epsilon = 2.0;

    double idiam = collide->extract(isp,"diam");
    double jdiam = collide->extract(jsp,"diam");
    double iomega = collide->extract(isp,"omega");
    double jomega = collide->extract(jsp,"omega");
    double itref = collide->extract(isp,"tref");
    double jtref = collide->extract(jsp,"tref");

    double pre_vibdof_i = species[isp].vibdof;
    double pre_vibdof_j = species[jsp].vibdof;
    double pre_ave_vibdof = (species[isp].vibdof + species[jsp].vibdof)/2.;
    double diam = 0.5 * (idiam+jdiam);
    double omega = 0.5 * (iomega+jomega);
    double tref = 0.5 * (itref+jtref);
    double mr = species[isp].mass * species[jsp].mass /
        (species[isp].mass + species[jsp].mass);
    double sigma = MY_PI*diam*diam;

    // average DOFs participating in the reaction

    double z = 0.0;  
    if (r->type == EXCHANGE) z = 0.0;
    if (r->type == DISSOCIATION) z = 1.0;
    if (r->type == IONIZATION) z = 0.0;
    
    // add additional coeff for effective DOF

    double c1 = MY_PIS*epsilon*r->coeff[2]/(2.0*sigma) *
      sqrt(mr/(2.0*update->boltz*tref)) *
      pow(tref,1.0-omega)/pow(update->boltz,r->coeff[3]-1.0+omega) *
      tgamma(z+2.5-omega)/tgamma(z+r->coeff[3]+1.5);
    double c2 = r->coeff[3] - 1.0 + omega;

    r->coeff[2] = c1;
    r->coeff[3] = c2;
    r->coeff[5] = z + 1.5 - omega; 
  }
}

/* ---------------------------------------------------------------------- */

void ReactBird::readfile(char *fname) 
{
  int n,n1,n2,eof;
  char line1[MAXLINE],line2[MAXLINE];
  char *word;
  OneReaction *r;

  // proc 0 opens file

  if (comm->me == 0) {
    fp = fopen(fname,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open reaction file %s",fname);
      error->one(FLERR,str);
    }
  }

  // read reactions one at a time and store their info in rlist

  while (1) {
    if (comm->me == 0) eof = readone(line1,line2,n1,n2);
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    
    MPI_Bcast(&n1,1,MPI_INT,0,world);
    MPI_Bcast(&n2,1,MPI_INT,0,world);
    MPI_Bcast(line1,n1,MPI_CHAR,0,world);
    MPI_Bcast(line2,n2,MPI_CHAR,0,world);

    if (nlist == maxlist) {
      maxlist += DELTALIST;
      rlist = (OneReaction *) 
        memory->srealloc(rlist,maxlist*sizeof(OneReaction),"react/tce:rlist");
      for (int i = nlist; i < maxlist; i++) {
        r = &rlist[i];
        r->nreactant = r->nproduct = 0;
        r->id_reactants = new char*[MAXREACTANT];
        r->id_products = new char*[MAXPRODUCT];
        r->reactants = new int[MAXREACTANT];
        r->products = new int[MAXPRODUCT];
        r->coeff = new double[MAXCOEFF];
      }
    }

    r = &rlist[nlist];
    r->initflag = 0;

    int side = 0;
    int species = 1;

    word = strtok(line1," \t\n");

    while (1) {
      if (!word) {
        if (side == 0) error->all(FLERR,"Invalid reaction formula in file");
        if (species) error->all(FLERR,"Invalid reaction formula in file");
        break;
      }
      if (species) {
        species = 0;
        if (side == 0) {
          n = strlen(word) + 1;
          r->id_reactants[r->nreactant] = new char[n];
          strcpy(r->id_reactants[r->nreactant],word);
          r->nreactant++;
        } else {
          n = strlen(word) + 1;
          r->id_products[r->nproduct] = new char[n];
          strcpy(r->id_products[r->nproduct],word);
          r->nproduct++;
        }
      } else {
        species = 1;
        if (strcmp(word,"+") == 0) {
          word = strtok(NULL," \t\n");
          continue;
        }
        if (strcmp(word,"-->") != 0) 
          error->all(FLERR,"Invalid reaction formula in file");
        side = 1;
      }
      word = strtok(NULL," \t\n");
    }

    word = strtok(line2," \t\n");
    if (!word) error->all(FLERR,"Invalid reaction type in file");
    if (word[0] == 'D' || word[0] == 'd') r->type = DISSOCIATION;
    else if (word[0] == 'E' || word[0] == 'e') r->type = EXCHANGE;
    else if (word[0] == 'I' || word[0] == 'i') r->type = IONIZATION;
    else if (word[0] == 'R' || word[0] == 'r') r->type = RECOMBINATION;
    else error->all(FLERR,"Invalid reaction type in file");

    word = strtok(NULL," \t\n");
    if (!word) error->all(FLERR,"Invalid reaction style in file");
    if (word[0] == 'A' || word[0] == 'a') r->style = ARRHENIUS;
    else if (word[0] == 'Q' || word[0] == 'q') r->style = QUANTUM;
    else error->all(FLERR,"Invalid reaction style in file");

    if (r->style == ARRHENIUS || r->style == QUANTUM) r->ncoeff = 5;

    for (int i = 0; i < r->ncoeff; i++) {
      word = strtok(NULL," \t\n");
      if (!word) error->all(FLERR,"Invalid reaction coefficients in file");
      r->coeff[i] = atof(word);
    }

    word = strtok(NULL," \t\n");
    if (word) error->all(FLERR,"Invalid reaction coefficients in file");
    
    nlist++;
  }

  if (comm->me == 0) fclose(fp);
}

/* ----------------------------------------------------------------------
   read one reaction from file
   reaction = 2 lines
   return 1 if end-of-file, else return 0
------------------------------------------------------------------------- */

int ReactBird::readone(char *line1, char *line2, int &n1, int &n2) 
{
  char *eof;
  while (eof = fgets(line1,MAXLINE,fp)) {
    int pre = strspn(line1," \t\n");
    if (pre == strlen(line1) || line1[pre] == '#') continue;
    eof = fgets(line2,MAXLINE,fp);
    if (!eof) break;
    n1 = strlen(line1) + 1;
    n2 = strlen(line2) + 1;
    return 0;
  }

  return 1;
}
