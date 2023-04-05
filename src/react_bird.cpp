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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "react_bird.h"
#include "input.h"
#include "collide.h"
#include "update.h"
#include "particle.h"
#include "comm.h"
#include "modify.h"
#include "fix.h"
#include "fix_ambipolar.h"
#include "random_knuth.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{DISSOCIATION,EXCHANGE,IONIZATION,RECOMBINATION};  // other react files
enum{ARRHENIUS,QUANTUM};                               // other react files

#define MAXREACTANT 2
#define MAXPRODUCT 3
#define MAXCOEFF 7               // 5 in file, extra for pre-computation

#define MAXLINE 1024
#define DELTALIST 16

/* ---------------------------------------------------------------------- */

ReactBird::ReactBird(SPARTA *sparta, int narg, char **arg) :
  React(sparta, narg, arg)
{
  if (narg != 2) error->all(FLERR,"Illegal react tce or qk command");

  nlist = maxlist = 0;
  rlist = NULL;
  readfile(arg[1]);
  check_duplicate();

  tally_reactions = new bigint[nlist];
  tally_reactions_all = new bigint[nlist];
  tally_flag = 0;

  reactions = NULL;
  list_ij = NULL;
  sp2recomb_ij = NULL;
}

/* ---------------------------------------------------------------------- */

ReactBird::ReactBird(SPARTA *sparta) : React(sparta)
{
  rlist = NULL;
  reactions = NULL;
  list_ij = NULL;
  sp2recomb_ij = NULL;
}

/* ---------------------------------------------------------------------- */

ReactBird::~ReactBird()
{
  if (copy) return;

  delete [] tally_reactions;
  delete [] tally_reactions_all;

  if (rlist) {
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
      delete [] rlist[i].id;
    }
  }
  memory->destroy(rlist);

  memory->destroy(reactions);
  memory->destroy(list_ij);
  memory->destroy(sp2recomb_ij);
}

/* ---------------------------------------------------------------------- */

void ReactBird::init()
{
  tally_flag = 0;
  for (int i = 0; i < nlist; i++) tally_reactions[i] = 0;

  // convert species IDs to species indices
  // flag reactions as active/inactive depending on whether all species exist
  // mark recombination reactions inactive if recombflag_user = 0

  for (int m = 0; m < nlist; m++) {
    OneReaction *r = &rlist[m];
    r->active = 1;

    if (r->type == RECOMBINATION && recombflag_user == 0) {
      r->active = 0;
      continue;
    }

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

        // special case: recombination reaction with 2nd product = atom/mol

        if (r->type == RECOMBINATION && i == 1) {
          if (strcmp(r->id_products[i],"atom") == 0) {
            r->products[i] = -1;
            continue;
          } else if (strcmp(r->id_products[i],"mol") == 0) {
            r->products[i] = -2;
            continue;
          }
        }

        r->active = 0;
        break;
      }
    }
  }

  // count possible active reactions for each species pair
  // include J,I reactions in I,J list and vice versa
  // this allows collision pair I,J to be in either order in Collide

  memory->destroy(reactions);
  int nspecies = particle->nspecies;
  reactions = memory->create(reactions,nspecies,nspecies,
                             "react/bird:reactions");

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

  // allocate list_IJ = contiguous list of reactions for each IJ pair

  memory->destroy(list_ij);
  memory->create(list_ij,n,"react/bird:list_ij");

  // reactions[i][j].list = pointer into full list_ij vector

  int offset = 0;
  for (int i = 0; i < nspecies; i++)
    for (int j = 0; j < nspecies; j++) {
      reactions[i][j].list = &list_ij[offset];
      offset += reactions[i][j].n;
    }

  // reactions[i][j].list = indices of reactions for each species pair
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

    // symmetry parameter

    double epsilon = 1.0;
    if (isp == jsp) epsilon = 2.0;

    double diam = collide->extract(isp,jsp,"diam");
    double omega = collide->extract(isp,jsp,"omega");
    double tref = collide->extract(isp,jsp,"tref");

    // double pre_ave_vibdof = (species[isp].vibdof + species[jsp].vibdof)/2.0;
    double mr = species[isp].mass * species[jsp].mass /
        (species[isp].mass + species[jsp].mass);
    double sigma = MY_PI*diam*diam;

    // average DOFs participating in the reaction

    double z = r->coeff[0];

    // add additional coeff for effective DOF
    // added MAX() limit, 24Aug18

    double c1 = MY_PIS*epsilon*r->coeff[2]/(2.0*sigma) *
      sqrt(mr/(2.0*update->boltz*tref)) *
      pow(tref,1.0-omega)/pow(update->boltz,r->coeff[3]-1.0+omega) *
      tgamma(z+2.5-omega) / MAX(1.0e-6,tgamma(z+r->coeff[3]+1.5));
    double c2 = r->coeff[3] - 1.0 + omega;

    r->coeff[2] = c1;
    r->coeff[3] = c2;
    r->coeff[5] = z + 1.5 - omega;

    // add additional coeff for post-collision effective omega
    // mspec = post-collision species of the particle
    // aspec = post-collision species of the atom

    double momega,aomega;

    if (r->nproduct > 1) {
      int mspec = r->products[0];
      int aspec = r->products[1];

      if (species[mspec].rotdof < 2.0) {
        mspec = r->products[1];
        aspec = r->products[0];
      }

      int ncount = 0;
      if (mspec >= 0) {
        momega = collide->extract(mspec,mspec,"omega");
        ncount++;
      } else momega = 0.0;
      if (aspec >= 0) {
        aomega = collide->extract(aspec,aspec,"omega");
        ncount++;
      } else aomega = 0.0;

      r->coeff[6] = (momega+aomega) / ncount;

    } else {
      int mspec = r->products[0];
      momega = collide->extract(mspec,mspec,"omega");
      r->coeff[6] = momega;
    }
  }

  // set recombflag = 0/1 if any recombination reactions are defined & active
  // check for user disabling them is at top of this method

  recombflag = 0;
  for (int m = 0; m < nlist; m++) {
    if (!rlist[m].active) continue;
    if (rlist[m].type == RECOMBINATION) recombflag = 1;
  }

  if (!recombflag) return;

  // count how many IJ pairs have a recombination reaction
  // allocate sp2recomb_ij = contiguous list of reactions
  //   for all species for each IJ pair that has a recombination reaction

  OneReaction *r;

  int nij = 0;
  for (int i = 0; i < nspecies; i++)
    for (int j = 0; j < nspecies; j++) {
      int n = reactions[i][j].n;
      int *list = reactions[i][j].list;
      for (int m = 0; m < n; m++) {
        r = &rlist[list[m]];
        if (r->type == RECOMBINATION) {
          nij++;
          break;
        }
      }
    }

  memory->destroy(sp2recomb_ij);
  memory->create(sp2recomb_ij,nij*nspecies,"react/bird:sp2recomb_ij");

  // reactions[i][j].sp2recomb = pointer into full sp2recomb_ij vector

  offset = 0;
  for (int i = 0; i < nspecies; i++)
    for (int j = 0; j < nspecies; j++) {
      int n = reactions[i][j].n;
      int *list = reactions[i][j].list;
      if (offset < nij*nspecies)
        reactions[i][j].sp2recomb = &sp2recomb_ij[offset];
      else
        reactions[i][j].sp2recomb = NULL; // Needed for Kokkos
      for (int m = 0; m < n; m++) {
        r = &rlist[list[m]];
        if (r->type == RECOMBINATION) {
          offset += nspecies;
          break;
        }
      }
    }

  // loop over species K for each IJ pair
  // if the IJ pair has any recombination reactions,
  //   then fill in its reactions[i][j].sp2recomb entries,
  //   which are indices into rlist of specific recomb reaction
  //   for each 3rd particle species
  // if IJ pair has no recombination reactions, then do NOT set sp2recomb vec
  // matching recombination reaction is one most specific to species K
  // 4 levels of specificity from most to least, in 4 inner loops
  //   explicit K, K = atom/mol, any K, no match at all (sp2recomb = -1)

  int m;

  for (int i = 0; i < nspecies; i++)
    for (int j = 0; j < nspecies; j++) {
      int n = reactions[i][j].n;
      int *list = reactions[i][j].list;
      for (int k = 0; k < nspecies; k++) {
        for (m = 0; m < n; m++) {
          r = &rlist[list[m]];
          if (r->type != RECOMBINATION) continue;
          if (r->nproduct != 2 || r->products[1] < 0) continue;
          if (r->products[1] == k) {
            reactions[i][j].sp2recomb[k] = list[m];
            break;
          }
        }
        if (m < n) continue;

        for (m = 0; m < n; m++) {
          r = &rlist[list[m]];
          if (r->type != RECOMBINATION) continue;
          if (r->nproduct != 2 || r->products[1] >= 0) continue;
          if (r->products[1] == -1 && particle->species[k].vibdof == 0) {
            reactions[i][j].sp2recomb[k] = list[m];
            break;
          }
          if (r->products[1] == -2 && particle->species[k].vibdof > 0) {
            reactions[i][j].sp2recomb[k] = list[m];
            break;
          }
        }
        if (m < n) continue;

        for (m = 0; m < n; m++) {
          r = &rlist[list[m]];
          if (r->type != RECOMBINATION) continue;
          if (r->nproduct != 1) continue;
          reactions[i][j].sp2recomb[k] = list[m];
          break;
        }
        if (m < n) continue;

        for (m = 0; m < n; m++) {
          r = &rlist[list[m]];
          if (r->type != RECOMBINATION) continue;
          reactions[i][j].sp2recomb[k] = -1;
          break;
        }
      }
    }
}

/* ----------------------------------------------------------------------
   return 1 if any recombination reactions are defined for species pair ISP,JSP
   else return 0
   called from Collide::init(), after React::init() has been performed
------------------------------------------------------------------------- */

int ReactBird::recomb_exist(int isp, int jsp)
{
  if (reactions[isp][jsp].sp2recomb) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   check active reactions that include ambi ion or electron especies
   their format must be correct to work with ambi_reset()
   called after init() from collide::init()
------------------------------------------------------------------------- */

void ReactBird::ambi_check()
{
  int flag;
  OneReaction *r;

  // fix ambipolar must exist since collide caller extracted ambi vector/array

  int ifix;
  for (ifix = 0; ifix < modify->nfix; ifix++)
    if (strcmp(modify->fix[ifix]->style,"ambipolar") == 0) break;
  FixAmbipolar *afix = (FixAmbipolar *) modify->fix[ifix];
  int especies = afix->especies;
  int *ions = afix->ions;

  // loop over active reactions

  for (int i = 0; i < nlist; i++) {
    r = &rlist[i];
    if (!r->active) continue;

    // skip reaction if no ambipolar ions or electrons as reactant or product
    // r->products[j] can be < 0 for atom or mol

    flag = 0;
    for (int j = 0; j < r->nreactant; j++)
      if (r->reactants[j] == especies || ions[r->reactants[j]]) flag = 1;
    for (int j = 0; j < r->nproduct; j++) {
      if (r->products[j] < 0) continue;
      if (r->products[j] == especies || ions[r->products[j]]) flag = 1;
    }
    if (!flag) continue;

    // dissociation must match one of these orders
    // D: AB + e -> A + e + B
    // D: AB+ + e -> A+ + e + B

    flag = 1;

    if (r->type == DISSOCIATION) {
      if (r->nreactant == 2 && r->nproduct == 3) {
        if (ions[r->reactants[0]] == 0 && r->reactants[1] == especies &&
            ions[r->products[0]] == 0 && r->products[1] == especies &&
            ions[r->products[2]] == 0) flag = 0;
        else if (ions[r->reactants[0]] == 1 && r->reactants[1] == especies &&
                 ions[r->products[0]] == 1 && r->products[1] == especies &&
                 ions[r->products[2]] == 0) flag = 0;
      }
    }

    // ionization with 3 products must match this
    // I: A + e -> A+ + e + e

    else if (r->type == IONIZATION && r->nproduct == 3) {
      if (r->nreactant == 2 && r->nproduct == 3) {
        if (ions[r->reactants[0]] == 0 && r->reactants[1] == especies &&
            ions[r->products[0]] == 1 && r->products[1] == especies &&
            r->products[2] == especies) flag = 0;
      }
    }

    // ionization with 2 products must match this
    // I: A + B -> AB+ + e

    else if (r->type == IONIZATION && r->nproduct == 2) {
      if (r->nreactant == 2 && r->nproduct == 2) {
        if (ions[r->reactants[0]] == 0 && ions[r->reactants[1]] == 0 &&
            ions[r->products[0]] == 1 && r->products[1] == especies) flag = 0;
      }
    }

    // exchange must match one of these
    // E: AB+ + e -> A + B
    // E: AB+ + C -> A + BC+
    // E: C + AB+ -> A + BC+

    else if (r->type == EXCHANGE) {
      if (r->nreactant == 2 && r->nproduct == 2) {
        if (ions[r->reactants[0]] == 1 && r->reactants[1] == especies &&
            ions[r->products[0]] == 0 && ions[r->products[1]] == 0) flag = 0;
        else if (ions[r->reactants[0]] == 1 && ions[r->reactants[1]] == 0 &&
            ions[r->products[0]] == 0 && ions[r->products[1]] == 1) flag = 0;
        else if (ions[r->reactants[0]] == 0 && ions[r->reactants[1]] == 1 &&
            ions[r->products[0]] == 0 && ions[r->products[1]] == 1) flag = 0;
      }
    }

    // recombination must match one of these
    // R: A+ + e -> A
    // R: A + B+ -> AB+
    // R: A+ + B -> AB+

    else if (r->type == RECOMBINATION) {
      if (r->nreactant == 2 && r->nproduct == 1) {
        if (ions[r->reactants[0]] == 1 && r->reactants[1] == especies &&
            ions[r->products[0]] == 0) flag = 0;
        else if (ions[r->reactants[0]] == 0 && ions[r->reactants[1]] == 1 &&
            ions[r->products[0]] == 1) flag = 0;
        else if (ions[r->reactants[0]] == 1 && ions[r->reactants[1]] == 0 &&
            ions[r->products[0]] == 1) flag = 0;
      }
    }

    // flag = 1 means unrecognized reaction

    if (flag) {
      print_reaction_ambipolar(r);
      error->all(FLERR,"Invalid ambipolar reaction");
    }
  }
}

/* ---------------------------------------------------------------------- */

void ReactBird::readfile(char *fname)
{
  int n,n1,n2,eof;
  char line1[MAXLINE],line2[MAXLINE];
  char copy1[MAXLINE],copy2[MAXLINE];
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
        memory->srealloc(rlist,maxlist*sizeof(OneReaction),"react/bird:rlist");
      for (int i = nlist; i < maxlist; i++) {
        r = &rlist[i];
        r->nreactant = r->nproduct = 0;
        r->id_reactants = new char*[MAXREACTANT];
        r->id_products = new char*[MAXPRODUCT];
        r->reactants = new int[MAXREACTANT];
        r->products = new int[MAXPRODUCT];
        r->coeff = new double[MAXCOEFF];
        r->id = NULL;
      }
    }

    strcpy(copy1,line1);
    strcpy(copy2,line2);

    r = &rlist[nlist];
    r->initflag = 0;

    int side = 0;
    int species = 1;

    n = strlen(line1) - 1;
    r->id = new char[n+1];
    strncpy(r->id,line1,n);
    r->id[n] = '\0';

    word = strtok(line1," \t\n\r");

    while (1) {
      if (!word) {
        if (side == 0) {
          print_reaction(copy1,copy2);
          error->all(FLERR,"Invalid reaction formula in file");
        }
        if (species) {
          print_reaction(copy1,copy2);
          error->all(FLERR,"Invalid reaction formula in file");
        }
        break;
      }
      if (species) {
        species = 0;
        if (side == 0) {
          if (r->nreactant == MAXREACTANT) {
            print_reaction(copy1,copy2);
            error->all(FLERR,"Too many reactants in a reaction formula");
          }
          n = strlen(word) + 1;
          r->id_reactants[r->nreactant] = new char[n];
          strcpy(r->id_reactants[r->nreactant],word);
          r->nreactant++;
        } else {
          if (r->nreactant == MAXPRODUCT) {
            print_reaction(copy1,copy2);
            error->all(FLERR,"Too many products in a reaction formula");
          }
          n = strlen(word) + 1;
          r->id_products[r->nproduct] = new char[n];
          strcpy(r->id_products[r->nproduct],word);
          r->nproduct++;
        }
      } else {
        species = 1;
        if (strcmp(word,"+") == 0) {
          word = strtok(NULL," \t\n\r");
          continue;
        }
        if (strcmp(word,"-->") != 0) {
          print_reaction(copy1,copy2);
          error->all(FLERR,"Invalid reaction formula in file");
        }
        side = 1;
      }
      word = strtok(NULL," \t\n\r");
    }

    word = strtok(line2," \t\n\r");
    if (!word) {
      print_reaction(copy1,copy2);
      error->all(FLERR,"Invalid reaction type in file");
    }
    if (word[0] == 'D' || word[0] == 'd') r->type = DISSOCIATION;
    else if (word[0] == 'E' || word[0] == 'e') r->type = EXCHANGE;
    else if (word[0] == 'I' || word[0] == 'i') r->type = IONIZATION;
    else if (word[0] == 'R' || word[0] == 'r') r->type = RECOMBINATION;
    else {
      print_reaction(copy1,copy2);
      error->all(FLERR,"Invalid reaction type in file");
    }

    // check that reactant/product counts are consistent with type

    if (r->type == DISSOCIATION) {
      if (r->nreactant != 2 || r->nproduct != 3) {
        print_reaction(copy1,copy2);
        error->all(FLERR,"Invalid dissociation reaction");
      }
    } else if (r->type == EXCHANGE) {
      if (r->nreactant != 2 || r->nproduct != 2) {
        print_reaction(copy1,copy2);
        error->all(FLERR,"Invalid exchange reaction");
      }
    } else if (r->type == IONIZATION) {
      if (r->nreactant != 2 || (r->nproduct != 2 && r->nproduct != 3)) {
        print_reaction(copy1,copy2);
        error->all(FLERR,"Invalid ionization reaction");
      }
    } else if (r->type == RECOMBINATION) {
      if (r->nreactant != 2 || (r->nproduct != 1 && r->nproduct != 2)) {
        print_reaction(copy1,copy2);
        error->all(FLERR,"Invalid recombination reaction");
      }
    }

    word = strtok(NULL," \t\n\r");
    if (!word) {
      print_reaction(copy1,copy2);
      error->all(FLERR,"Invalid reaction style in file");
    }
    if (word[0] == 'A' || word[0] == 'a') r->style = ARRHENIUS;
    else if (word[0] == 'Q' || word[0] == 'q') r->style = QUANTUM;
    else {
      print_reaction(copy1,copy2);
      error->all(FLERR,"Invalid reaction style in file");
    }

    if (r->style == ARRHENIUS || r->style == QUANTUM) r->ncoeff = 5;

    for (int i = 0; i < r->ncoeff; i++) {
      word = strtok(NULL," \t\n\r");
      if (!word) {
        print_reaction(copy1,copy2);
        error->all(FLERR,"Invalid reaction coefficients in file");
      }
      r->coeff[i] = input->numeric(FLERR,word);
    }

    word = strtok(NULL," \t\n\r");
    if (word) {
      print_reaction(copy1,copy2);
      error->all(FLERR,"Too many coefficients in a reaction formula");
    }

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
  while ((eof = fgets(line1,MAXLINE,fp))) {
    size_t pre = strspn(line1," \t\n\r");
    if (pre == strlen(line1) || line1[pre] == '#') continue;
    eof = fgets(line2,MAXLINE,fp);
    if (!eof) break;
    n1 = strlen(line1) + 1;
    n2 = strlen(line2) + 1;
    return 0;
  }

  return 1;
}

/* ----------------------------------------------------------------------
   check for duplicates in list of reactions read from file
   error if any exist
------------------------------------------------------------------------- */

void ReactBird::check_duplicate()
{
  OneReaction *r,*s;

  for (int i = 0; i < nlist; i++) {
    r = &rlist[i];

    for (int j = i+1; j < nlist; j++) {
      s = &rlist[j];

      if (r->type != s->type) continue;
      if (r->style != s->style) continue;
      if (r->nreactant != s->nreactant) continue;
      if (r->nproduct != s->nproduct) continue;

      int reactant_match = 0;
      if (strcmp(r->id_reactants[0],s->id_reactants[0]) == 0 &&
          strcmp(r->id_reactants[1],s->id_reactants[1]) == 0)
        reactant_match = 1;
      else if (strcmp(r->id_reactants[0],s->id_reactants[1]) == 0 &&
               strcmp(r->id_reactants[1],s->id_reactants[0]) == 0)
        reactant_match = 2;
      if (!reactant_match) continue;

      int product_match = 0;
      if (r->nproduct == 1) {
        if (strcmp(r->id_products[0],s->id_products[0]) == 0)
          product_match = 1;
      } else if (r->nproduct >= 2) {
        if (strcmp(r->id_products[0],s->id_products[0]) == 0 &&
            strcmp(r->id_products[1],s->id_products[1]) == 0)
          product_match = 1;
        else if (strcmp(r->id_products[0],s->id_products[1]) == 0 &&
                 strcmp(r->id_products[1],s->id_products[0]) == 0)
          product_match = 2;
      }
      if (!product_match) continue;

      if (comm->me == 0) {
        printf("MATCH %d %d %d: %d\n",i,j,nlist,product_match);
        printf("MATCH %d %d %d %d\n",
               r->products[0],r->products[1],s->products[0],s->products[1]);
      }
      print_reaction(r);
      print_reaction(s);
      error->all(FLERR,"Duplicate reactions in reaction file");
    }
  }
}

/* ----------------------------------------------------------------------
   print reaction as read from file
   only proc 0 performs output
------------------------------------------------------------------------- */

void ReactBird::print_reaction(char *line1, char *line2)
{
  if (comm->me) return;
  printf("Bad reaction format:\n");
  printf("%s\n%s\n",line1,line2);
};

/* ----------------------------------------------------------------------
   print reaction as stored in rlist
   only proc 0 performs output
------------------------------------------------------------------------- */

void ReactBird::print_reaction(OneReaction *r)
{
  if (comm->me) return;
  printf("Bad reaction:\n");

  char type;
  if (r->type == DISSOCIATION) type = 'D';
  else if (r->type == EXCHANGE) type = 'E';
  else if (r->type == IONIZATION) type = 'I';
  else if (r->type == RECOMBINATION) type = 'R';

  char style;
  if (r->style == ARRHENIUS) style = 'A';
  else if (r->style == QUANTUM) style = 'Q';

  if (r->nproduct == 1)
    printf("  %c %c: %s + %s --> %s\n",type,style,
           r->id_reactants[0],r->id_reactants[1],
           r->id_products[0]);
  else if (r->nproduct == 2)
    printf("  %c %c: %s + %s --> %s %s\n",type,style,
           r->id_reactants[0],r->id_reactants[1],
           r->id_products[0],r->id_products[1]);
  else if (r->nproduct == 3)
    printf("  %c %c: %s + %s --> %s %s %s\n",type,style,
           r->id_reactants[0],r->id_reactants[1],
           r->id_products[0],r->id_products[1],r->id_products[2]);
};

/* ----------------------------------------------------------------------
   print reaction as stored in rlist
   only proc 0 performs output
------------------------------------------------------------------------- */

void ReactBird::print_reaction_ambipolar(OneReaction *r)
{
  if (comm->me) return;
  printf("Bad ambipolar reaction format:\n");
  printf("  type %d style %d\n",r->type,r->style);
  printf("  nreactant %d:",r->nreactant);
  for (int i = 0; i < r->nreactant; i++)
    printf(" %s",r->id_reactants[i]);
  printf("\n");
  printf("  nproduct %d:",r->nproduct);
  for (int i = 0; i < r->nproduct; i++)
    printf(" %s",r->id_products[i]);
  printf("\n");
  printf("  ncoeff %d:",r->ncoeff);
  for (int i = 0; i < r->ncoeff; i++)
    printf(" %g",r->coeff[i]);
  printf("\n");
};

/* ----------------------------------------------------------------------
   return reaction ID = chemical formula
------------------------------------------------------------------------- */

char *ReactBird::reactionID(int m)
{
  return rlist[m].id;
};

/* ----------------------------------------------------------------------
   return tally associated with a reaction
------------------------------------------------------------------------- */

double ReactBird::extract_tally(int m)
{
  if (!tally_flag) {
    tally_flag = 1;
    MPI_Allreduce(tally_reactions,tally_reactions_all,nlist,
                  MPI_SPARTA_BIGINT,MPI_SUM,world);
  }

  return 1.0*tally_reactions_all[m];
};
