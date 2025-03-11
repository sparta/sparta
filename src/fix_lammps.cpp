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

/* ----------------------------------------------------------------------
   Contributing author: Tim Linke (UC Davis, LLNL)
------------------------------------------------------------------------- */

#include <vector>
#include <cmath>
#include <limits>
#include <set>
#include <algorithm>

#include <string>
#include "string.h"
#include "fix_lammps.h"
#include "grid.h"
#include "surf.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "random_mars.h"
#include "random_knuth.h"

#include "compute.h"

// these are LAMMPS files
#include "lammps.h"
#include "../../lammps/src/input.h"
#include "../../lammps/src/atom.h"
#include "../../lammps/src/library.h"
#include "../../lammps/src/modify.h"
#include "../../lammps/src/fix.h"
#include "../../lammps/src/force.h"
#include "../../lammps/src/compute.h"


using namespace SPARTA_NS;

enum{INT,DOUBLE};

FixLAMMPS::FixLAMMPS(SPARTA *sparta, int narg, char **arg) : Fix(sparta, narg, arg)
{
  //call is "fix ID lammps group-ID Nevery c_nrho keywords ..."
  if (narg < 5) error->all(FLERR,"Illegal fix lammps command. Doublecheck input arguments.");

  if (surf->implicit)
    error->all(FLERR,"Cannot use fix lammps with implicit surfs");
  if (surf->distributed)
    error->all(FLERR,"Cannot use fix lammps with distributed surfs");

  int igroup = grid->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Could not find group ID of grid cells.");
  groupbit = grid->bitmask[igroup];

  nevery = atoi(arg[3]);
  if (nevery < 0) error->all(FLERR, "Please provide valid frequency for fix lammps");

  if (strncmp(arg[4],"c_",2) == 0) {
    int n = strlen(arg[4]);
    id_nrho = new char[n];
    strcpy(id_nrho,&arg[4][2]);
    // error checks
    icompute = modify->find_compute(id_nrho);
    c_nrho = modify->compute[icompute];
    if (icompute < 0) error->all(FLERR,"Could not find fix lammps compute ID");
  } else error->all(FLERR,"Invalid input for fix lammps compute ID");

  /*  code for custom keywords, such as filename input, computes, etc. */
  // int iarg = 5;
  // while (iarg < narg) {
  //   if (strncmp(arg[iarg],"filename",8) == 0) {
  //   // in case lammps input file should be an input to the fix:
  //   int n = strlen(arg[iarg]) + 1;
  //   filename = new char[n];
  //   strcpy(filename,arg[4]);
  //   }
  //   iarg++;
  // }

  // create per-surf temperature vector
  char *id_custom = new char[12];
  strcpy(id_custom,"temperature");
  tindex = surf->add_custom(id_custom,DOUBLE,0);
  delete [] id_custom;
  // trigger setup of list of owned surf elements belonging to surf group
  firstflag = 1;
  // initialize data structure
  tvector_me = NULL;
  // number of grid cells
  ncells = grid->nlocal;
  dimension = domain->dimension;
  nvalid = nextvalid();
  tvector_me = NULL;
}

/* ---------------------------------------------------------------------- */

FixLAMMPS::~FixLAMMPS()
{
  memory->destroy(tvector_me);
  surf->remove_custom(tindex);
  delete random;
}

/* ---------------------------------------------------------------------- */

int FixLAMMPS::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLAMMPS::init()
{
  // Initialize random number generator
  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,5,100);

  // one-time initialization of temperature for all surfs in custom vector
  if (!firstflag) return;
  firstflag = 0;
  int nlocal = surf->nlocal;
  tvector = surf->edvec[surf->ewhich[tindex]];
  for (int i = 0; i < nlocal; i++) {
    tvector[i] = 273.0;
  }
  
  // allocate per-surf vector for explicit all surfs
  memory->create(tvector_me,nlocal,"lammps:tvector_me");
}

void FixLAMMPS::end_of_step()
{
  // skip LAMMPS check if not on timestep that's multiple of Nevery
  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;

  // Duplicate communicator
  MPI_Comm comm_lammps;  // New communicator
  MPI_Comm_dup(MPI_COMM_WORLD, &comm_lammps);
  int rank,nprocs;
  MPI_Comm_rank(comm_lammps,&rank);
  MPI_Comm_size(comm_lammps,&nprocs);

  // Naming of MD input file
  sprintf(filename, "MDfiles/in.lammps.proc%d_", rank);

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  // Grid cell densities from compute
  double *gridrho;
  gridrho = c_nrho->vector_grid;
  if (gridrho == NULL)
  {
    error->all(FLERR, "Cannot find fix nrho info\n");
  }
  
  int nlocal = surf->nlocal;
  memset(tvector_me,0,nlocal*sizeof(double));

  // Set LAMMPS parameters
  double surftemp = 0.0;
  for (int i = 0; i < nlocal; i++) {
    surftemp += tvector[i];
  }
  surftemp /= nlocal;

  int coupleCell;
  int coupledProc;
  int local_fnum;
  int maxCells; // (all procs need to loop over the maximum number of cells of one proc to be able to assist in a coupling)
  MPI_Allreduce(&ncells, &maxCells, 1, MPI_INT, MPI_MAX, comm_lammps);

  for (int i = 0; i < maxCells; i++)
  {
    coupleCell = 0;
    coupledProc = 0;
    double xLength, yLength, yBox, xBox;
    local_fnum = 19; //design choice: MD atoms representing one DSMC particle

    if ((i < ncells) && (cinfo[i].mask & groupbit)) // apply coupling to selected cell group
    {
      coupleCell = 1;
      coupledProc = rank;
      fprintf(screen,"DSMC-MD coupling invoked for cell %d. Running LAMMPS ...\n", cells[i].id);

      // Check for DSMC particles in cell
      int particle_iterator = cinfo[i].first;
      if (particle_iterator == -1) {
        fprintf(screen, "No particles in cell %d, skipping.\n",cells[i].id);
        coupleCell = 2;
      }

      /****************************
      1. Lifting Operator: From SPARTA to LAMMPS, initialize microscopic properties derived from SPARTA
      2. Run LAMMPS simulation
      3. Restricting Operator: From LAMMPS to SPARTA, calculate relevant quantities to pass to SPARTA
      *****************************/

      // LIFTING OPERATOR -------------------------------------------------

      //Obtain SPARTA data
      double gridnrho = gridrho[cells[i].ilocal] * 1.0e28; // local number density
      numPar = (int) cinfo[i].count; // number of simulated particles in cell
      if (numPar > 27)
      {
        //LAMMPS is restricted to <32 groups
        fprintf(screen,"WARNING: Number of DSMC particles per cell is exceeding LAMMPS capabilities (more than 27 DSMC particles per cell).\n");
        coupleCell = 2;
      }
      int numAtoms = local_fnum * numPar;
      if (numAtoms > 1e5) fprintf(screen,"WARNING: High density detected in fix lammps. Number of MD atoms is %d. Reduce local_fnum\n",numAtoms);
      if (numAtoms < 10) fprintf(screen,"WARNING: Low density detected in fix lammps. Number of MD atoms is %d. Increase local_fnum\n",numAtoms);
      xLength = (cells + i)->hi[0] - (cells + i)->lo[0];
      yLength = (cells + i)->hi[1] - (cells + i)->lo[1];
      double AR = xLength/yLength; //aspect ratio of DSMC cell
      double gridnrhoA = gridnrho * 1e-30; //convert number density to 1/Angstrom^3

      // LAMMPS Simulation Parameters
      double zBox = 1.0; //in Angstrom
      yBox = sqrt(numAtoms / (AR * gridnrho * 1e-30 * zBox)); //in Angstrom
      double lattice_constant = 7.0; //in Angstrom
      double lattice_constant_surf = 2.55; // lattice constant of wall //3.52
      int lattice_length = (int) (yBox / lattice_constant_surf); //ensure MD length is an even multiple of lattice length
      yBox = lattice_length * lattice_constant_surf;
      xBox = AR * yBox; //in Angstrom
      
      // Obtain radius for sphere that will cover hexagonal lattice in atom initialization
      double radius = 0.0;
      radius = lattice_constant * sqrt( sqrt(3)*local_fnum / (2.0*M_PI) );

      int timesteps_equi = 10000;
      int timesteps = 10000;
      
      // create individual LAMMPS input file for every cell
      sprintf(filenameS, "%s%d", filename, cells[i].id);
      
      // fill input file with SPARTA data
      FILE *fp = fopen(filenameS,"w");
      if (fp == NULL) error->all(FLERR,"Could not open LAMMPS input script in fix lammps");
      if (screen == NULL) error->all(FLERR, "Screen output invalid in fix lammps");

      // Create LAMMPS input file
      fprintf(fp, "# LAMMPS input script created by SPARTA\n\n");
      fprintf(fp, "# Number of DSMC particles: %d\n",numPar);
      fprintf(fp, "# Number density nrho = %e\n", gridnrho);

      fprintf(fp,"units metal\n");
      fprintf(fp,"atom_style atomic\n\n");
      fprintf(fp,"dimension %d\n",dimension);
      fprintf(fp,"boundary f p p\n");
      fprintf(fp,"neighbor 3.0 bin\n\n");

      fprintf(fp, "region box block 0.0 %g 0.0 %g %g %g\n", xBox, yBox, -zBox*0.5, zBox*0.5);
      fprintf(fp, "create_box 2 box\n\n");
      fprintf(fp, "lattice hex %e\n\n", lattice_constant);

      // LIFTING OPERATOR: for every DSMC particle, create a batch of MD particles in space
      // MD particles are initialized in a hexagonal lattice around the DSMC particle's position
      double xMD, yMD, vxMD, vyMD;
      std::string regions;
      particle_iterator = cinfo[i].first;

      for (int p = 0; p < numPar; p++)
      {
        xMD = (particle->particles[particle_iterator].x[0] - (cells + i)->lo[0]) * (0.9*xBox / xLength); //moved to origin and scaled
        yMD = (particle->particles[particle_iterator].x[1] - (cells + i)->lo[1]) * (0.9*yBox / yLength);

        fprintf(fp, "region %d sphere %g %g %g %g units box\n", p+1, xMD, yMD, 0.0, radius);

        particle_iterator = particle->next[particle_iterator];
        regions += std::to_string(p+1);
        regions += " ";
      }

      fprintf(fp, "region flow union %d %s\n", numPar, regions.c_str());
      fprintf(fp, "create_atoms 2 region flow\n\n");

      fprintf(fp, "lattice hex %e\n", lattice_constant_surf);
      fprintf(fp, "region wall block %g %g %g %g %g %g units box\n", 0.9*xBox, xBox, 0.0, yBox, -zBox*0.5, zBox*0.5);
      fprintf(fp, "region sink block %g %g %g %g %g %g units box\n", 0.95*xBox, xBox, 0.0, yBox, -zBox*0.5, zBox*0.5);
      fprintf(fp, "region side block %g %g %g %g %g %g units box\n", 0.9*xBox, 0.95*xBox, 0.0, yBox, -zBox*0.5, zBox*0.5);
      fprintf(fp, "create_atoms 1 region wall\n\n");

      fprintf(fp, "# Material properties\n");
      fprintf(fp, "mass 1 58.710\n");
      fprintf(fp, "mass 2 28.96\n\n");

      fprintf(fp, "# Define potentials: EAM for wire, LJ for air\n");
      fprintf(fp, "pair_style hybrid eam/alloy lj/cut 2.5  # Hybrid potential for different interactions\n");
      fprintf(fp, "pair_coeff * * eam/alloy MDfiles/Ni_v6_2.0.eam Ni NULL\n");
      fprintf(fp, "pair_coeff 2 2 lj/cut 0.0091 3.55 # DOI: 10.1007/s11085-009-9180-z\n");
      fprintf(fp, "pair_coeff 1 2 lj/cut 0.3308 2.92 2.5 # Lorentz-Berthelot\n\n");

      fprintf(fp, "# Define groups\n");
      fprintf(fp, "group wire region wall\n");
      fprintf(fp, "group heatsink region sink\n");
      fprintf(fp, "group surfside region side\n");
      fprintf(fp, "group air region flow\n\n");
      for (int p = 0; p < numPar; p++) {
        fprintf(fp, "group %d region %d\n", p+1, p+1); // keep track of atoms per DSMC particle
      }

      fprintf(fp, "# Set wire velocity\n");
      fprintf(fp, "velocity wire create %g 12345 mom yes dist gaussian\n\n", surftemp);

      fprintf(fp, "timestep 5e-4\n");


      fprintf(fp, "# Equilibration of wire\n");
      fprintf(fp, "fix wall_boundary all wall/reflect xlo EDGE\n");
      fprintf(fp, "fix wall_lj wire wall/lj126 xhi EDGE 0.2314 2.49 6.5\n");
      fprintf(fp, "fix init wire nvt temp %g %g $(10.0*dt)\n\n", surftemp, surftemp);

      fprintf(fp, "# Compute wire temperature\n");
      fprintf(fp, "compute temp_wire wire temp\n\n");

      fprintf(fp, "# Output to file\n");
      fprintf(fp, "#fix wire_temp_output wire ave/time 1000 1 1000 c_temp_wire file MDfiles/MDwireTemp.txt\n\n");
      fprintf(fp, "fix wire_temp_ave wire ave/time %d %d %d c_temp_wire #file MDfiles/finalTemp.txt\n", 100, 200, timesteps_equi+timesteps); //sound averaging

      fprintf(fp, "thermo 0\n");
      fprintf(fp, "thermo_modify flush yes lost ignore\n");
      fprintf(fp, "run %d\n\n", timesteps_equi);

      fprintf(fp, "# Set DSMC velocities\n");
      particle_iterator = cinfo[i].first; //reset particle index
      for (int p = 0; p < numPar; p++)
      {
        vxMD = particle->particles[particle_iterator].v[0]*0.01; //converted from m/s to Angstrom/ps
        vyMD = particle->particles[particle_iterator].v[1]*0.01;
        fprintf(fp,"velocity %d set %g %g %g units box\n", p+1, vxMD, vyMD, 0.0);
        particle_iterator = particle->next[particle_iterator];
      }
      fprintf(fp, "\n");

      fprintf(fp, "timestep 1e-3\n");

      fprintf(fp, "# Fixes\n");
      fprintf(fp, "unfix init\n");
      fprintf(fp, "fix source heatsink nvt temp %g %g $(10.0*dt)\n", surftemp, surftemp);
      fprintf(fp, "fix cooling surfside nve\n");
      fprintf(fp, "fix intgl air nve\n");
      fprintf(fp,"fix 2d all enforce2d\n\n");

      fprintf(fp, "compute restricting air property/atom id x y vx vy\n");

      //fprintf(fp, "dump 2 all image 1000 MDfiles/lmp.image.*.ppm type type adiam 1.5 zoom 1.8\n");
      //fprintf(fp, "dump_modify  2 pad 7\n\n");

      fprintf(fp, "#Output settings\n");
      fprintf(fp, "thermo_style custom step time temp c_temp_wire\n");
      fprintf(fp, "thermo 4000  # Output every 100 steps\n");
      fprintf(fp, "thermo_modify flush yes lost ignore\n");
      fprintf(fp, "# Run simulation\n");

      fprintf(fp, "run %d\n", timesteps);


      fclose(fp);
    }

    //agree on coupled proc
    MPI_Barrier(comm_lammps);
    int maxProc;
    MPI_Allreduce(&coupledProc, &maxProc, 1, MPI_INT, MPI_MAX, comm_lammps);
    coupledProc = maxProc;
    MPI_Barrier(comm_lammps);
    MPI_Bcast(&coupleCell, 1, MPI_INT, coupledProc, comm_lammps); //communicate coupled cell
    MPI_Barrier(comm_lammps);
    
    if (coupleCell == 1)
    {
      // Inform all procs of LAMMPS input file
      int n = strlen(filenameS) + 1;
      MPI_Bcast(&n,1,MPI_INT,coupledProc,comm_lammps);
      MPI_Bcast(&filenameS, n, MPI_CHAR, coupledProc, comm_lammps);
      MPI_Comm_rank(comm_lammps, &rank); //ensure ranks are owned

      // RUN LAMMPS ----------------------------------------------------------------

      // instantiate LAMMPS
      LAMMPS_NS::LAMMPS *lmp;
      const char *lmpargv[] {"liblammps", "-screen", "MDfiles/MDlogfile.txt", "-log", "none"};
      int lmpargc = 5;

      lmp = new LAMMPS_NS::LAMMPS(lmpargc, (char **)lmpargv, comm_lammps);
      lmp->input->file(filenameS);

      // RESTRICTING OPERATOR ------------------------------------------------------
      // Sort atoms according to nearest neighbor criterion
      // sum up local atoms into a DSMC particle

      int natoms = static_cast<int> (lmp->atom->natoms); //total number of LAMMPS particles
      int nlocalLMP = lmp->atom->nlocal;
      auto lmp_compute = lmp->modify->get_compute_by_id("restricting");
      if (!lmp_compute)
      {
        error->all(FLERR,"Compute restricting cannot be found!");
      }
      if (!lmp_compute->peratom_flag) error->all(FLERR,"Compute restricting is not peratom!");
      lmp_compute->compute_peratom();
      double** lmp_props = lmp_compute->array_atom;
      if (!lmp_props)
      {
        error->all(FLERR,"Compute restricting has no info!");
      }

      int nAirAtoms_local = 0;
      for (int f = 0; f < nlocalLMP; f++) {
        if (lmp_props[f][1] != 0.0) {
            nAirAtoms_local++;
        }
      }

      int aindex = 0;
      int* IDs_local = new int[nAirAtoms_local];
      double* x_local = new double[nAirAtoms_local];
      double* y_local = new double[nAirAtoms_local];
      double* vx_local = new double[nAirAtoms_local];
      double* vy_local = new double[nAirAtoms_local];
      for (int f = 0; f < nlocalLMP; f++) {
        if (lmp_props[f][0] != 0.0) {
            IDs_local[aindex] = lmp_props[f][0];
            x_local[aindex] = lmp_props[f][1];
            y_local[aindex] = lmp_props[f][2];
            vx_local[aindex] = lmp_props[f][3];
            vy_local[aindex] = lmp_props[f][4];
            aindex++;
        }
      }
      int* nAirAtoms_all = new int[nprocs];
      MPI_Gather(&nAirAtoms_local, 1, MPI_INT, nAirAtoms_all, 1, MPI_INT, coupledProc, comm_lammps);

      int* displacements = new int[nprocs];
      int nAirAtoms = 0;
      if (rank == coupledProc) {
          for (int f = 0; f < nprocs; f++) {
              displacements[f] = nAirAtoms;
              nAirAtoms += nAirAtoms_all[f];
          }
      }

      int* IDs = nullptr;
      double* xLMP = nullptr;
      double* yLMP = nullptr;
      double* vxLMP = nullptr;
      double* vyLMP = nullptr;
      if (rank == coupledProc) {
          IDs = new int[nAirAtoms];
          xLMP = new double[nAirAtoms];
          yLMP = new double[nAirAtoms];
          vxLMP = new double[nAirAtoms];
          vyLMP = new double[nAirAtoms];
      }
      MPI_Gatherv(IDs_local, nAirAtoms_local, MPI_INT, IDs, nAirAtoms_all, displacements, MPI_INT, coupledProc, comm_lammps);
      MPI_Gatherv(x_local, nAirAtoms_local, MPI_DOUBLE, xLMP, nAirAtoms_all, displacements, MPI_DOUBLE, coupledProc, comm_lammps);
      MPI_Gatherv(y_local, nAirAtoms_local, MPI_DOUBLE, yLMP, nAirAtoms_all, displacements, MPI_DOUBLE, coupledProc, comm_lammps);
      MPI_Gatherv(vx_local, nAirAtoms_local, MPI_DOUBLE, vxLMP, nAirAtoms_all, displacements, MPI_DOUBLE, coupledProc, comm_lammps);
      MPI_Gatherv(vy_local, nAirAtoms_local, MPI_DOUBLE, vyLMP, nAirAtoms_all, displacements, MPI_DOUBLE, coupledProc, comm_lammps);
      delete[] displacements;
      delete[] IDs_local;
      delete[] x_local;
      delete[] y_local;
      delete[] vx_local;
      delete[] vy_local;
      delete[] nAirAtoms_all;

      auto lmp_fix = lmp->modify->get_fix_by_id("wire_temp_ave");
      if (!lmp_fix)
      {
        error->all(FLERR,"Fix wire_temp_ave cannot be found!");
      }

      double temp_lmp = (double) lmp_fix->compute_scalar();
      if (temp_lmp < 1e-5 || !temp_lmp)
      {
        fprintf(screen, "Error extracting data from fix lammps. Defaulting to calculated temperature.\n");
        temp_lmp = surftemp;
      }

      if(rank == coupledProc)
      {
        int particle_iterator = cinfo[i].first; //DSMC particle index
      
        // --------------------------------------------------------------------------------
        // Sort atoms according to their spatial distribution
        int nExpt = local_fnum * numPar;
        if (nAirAtoms < nExpt) 
        {
          local_fnum = (int) (nAirAtoms / numPar);
          nExpt = (local_fnum) * numPar;
        }
        bool* visited = new bool[nExpt]();
        int** groups = new int*[nExpt];
        int* groupSizes = new int[nExpt]();

        int groupCount = 0;

        for (int f = 0; f < nExpt; f++) 
        {
            if (visited[f]) continue;
            
            groups[groupCount] = new int[local_fnum];
            groups[groupCount][0] = f;
            visited[f] = true; 
            int groupIndex = 1;

            // Find the nearest neighbors
            std::pair<double, int>* distances = new std::pair<double, int>[nExpt];
            int distCount = 0;

            for (int j = 0; j < nExpt; j++) {
                if (!visited[j] && f != j) {
                    double dist = (xLMP[f] - xLMP[j]) * (xLMP[f] - xLMP[j]) + (yLMP[f] - yLMP[j]) * (yLMP[f] - yLMP[j]);
                    distances[distCount++] = std::make_pair(dist, j);
                }
            }

            // Sort by distance (ascending order)
            sort(distances, distances + distCount);

            // Add the closest particles to the group (up to 'local_fnum')
            for (int k = 0; k < local_fnum - 1 && k < distCount; k++) {
                int neighborIndex = distances[k].second;
                groups[groupCount][groupIndex++] = neighborIndex;
                visited[neighborIndex] = true;
            }

            groupSizes[groupCount] = groupIndex;
            groupCount++;

            delete[] distances;
        }

        // Reassign the xLMP array based on the groups
        int* new_IDs = new int[nExpt];
        double* new_xLMP = new double[nExpt];
        double* new_yLMP = new double[nExpt];
        double* new_vxLMP = new double[nExpt];
        double* new_vyLMP = new double[nExpt];
        int currentIndex = 0;

        for (int g = 0; g < groupCount; g++) {
            for (int j = 0; j < groupSizes[g]; j++) {
                int atomIndex = groups[g][j];
                // Reassign positions for each atom based on the group order
                new_IDs[currentIndex] = IDs[atomIndex];
                new_xLMP[currentIndex] = xLMP[atomIndex];
                new_yLMP[currentIndex] = yLMP[atomIndex];
                new_vxLMP[currentIndex] = vxLMP[atomIndex];
                new_vyLMP[currentIndex] = vyLMP[atomIndex];
                currentIndex++;
            }
            delete[] groups[g];
        }

        // Copy the new_xLMP back to xLMP
        for (int f = 0; f < nExpt; f++) {
            IDs[f] = new_IDs[f];
            xLMP[f] = new_xLMP[f];
            yLMP[f] = new_yLMP[f];
            vxLMP[f] = new_vxLMP[f];
            vyLMP[f] = new_vyLMP[f];
        }

        delete[] visited;
        delete[] groups;
        delete[] groupSizes;
        delete[] new_IDs;
        delete[] new_xLMP;
        delete[] new_yLMP;
        delete[] new_vxLMP;
        delete[] new_vyLMP;
        //--------------------------------------------------------------------------------
        double x = 0, y = 0, vx = 0, vy = 0;
        int count = local_fnum;

        // Per DSMC particle, find the center-of-mass average velocity of MD particles
        for (int a = 0; a < numPar; a++)
        {
          for (int b = (count - local_fnum); b < count; b++)
          {
            x += xLMP[b];
            y += yLMP[b];
            vx += vxLMP[b];
            vy += vyLMP[b];
          }
          count += local_fnum;
          // Assign averaged values to DSMC particles and scale back to DSMC domain
          particle->particles[particle_iterator].x[0] = ((x / local_fnum) * (xLength / (0.9*xBox))) + (cells + i)->lo[0];
          particle->particles[particle_iterator].x[1] = ((y / local_fnum) * (yLength / (0.9*yBox))) + (cells + i)->lo[1];
          particle->particles[particle_iterator].v[0] = vx / local_fnum * 100.0; // converted from Angstrom/ps to m/s
          particle->particles[particle_iterator].v[1] = vy / local_fnum * 100.0;
          x = y = vx = vy = 0; // reset averages
          particle_iterator = particle->next[particle_iterator];
        }
        // Assign surface temperature to DSMC surface
        for (int f = 0; f < nlocal; f++) {
          tvector_me[f] = temp_lmp;;
        }
      }
      MPI_Barrier(comm_lammps);
      if (IDs) {delete[] IDs; IDs = nullptr;}
      if (xLMP) {delete[] xLMP; xLMP = nullptr;}
      if (yLMP) {delete[] yLMP; yLMP = nullptr;}
      if (vxLMP) {delete[] vxLMP; vxLMP = nullptr;}
      if (vyLMP) {delete[] vyLMP; vyLMP = nullptr;}
      delete lmp;
    }
    else if (coupleCell == 2)
    {
      MPI_Comm_rank(comm_lammps, &rank); //ensure ranks are owned
      if(rank == coupledProc)
      {
        for (int f = 0; f < nlocal; f++) {
          tvector_me[f] = surftemp;
        }
      }
    }
    MPI_Barrier(comm_lammps);
  }
  nvalid += nevery;
  int size;
  MPI_Comm_size(comm_lammps, &size);
  MPI_Allreduce(tvector_me,tvector,nlocal,MPI_DOUBLE,MPI_SUM,comm_lammps);
  MPI_Comm_free(&comm_lammps);
  MPI_Barrier(MPI_COMM_WORLD);
}

bigint FixLAMMPS::nextvalid()
{
  bigint nvalid = (update->ntimestep/nevery)*nevery + nevery;
  if (nvalid < update->ntimestep) nvalid += nevery;
  return nvalid;
}