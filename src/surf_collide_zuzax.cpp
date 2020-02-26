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
#include "stdlib.h"
#include "string.h"
#include "surf_collide_zuzax.h"
#include "surf.h"
#include "surf_react.h"
#include "input.h"
#include "variable.h"
#include "particle.h"
#include "domain.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "random_mars.h"
#include "random_park.h"
#include "math_const.h"
#include "math_extra.h"
#include "error.h"

#ifdef USE_ZSURF

using namespace SPARTA_NS;
using namespace MathConst;

#include "zuzax/base/xml.h"
#include "zuzax/numerics/SolnMaterialDomain.h"
#include "zuzax/multiphase/PhaseList.h"
#include "zuzax/kinetics.h"
#include "zuzax/zeroD/ReactorNetDAE.h"
#include "zuzax/zerodim.h"
#include "zuzax/zeroD/IdealGasConstPressureReactor.h"
#include "zuzax/zeroD/RBdry.h"

#include "zuzax/zeroD/SubstrateElement.h"
#include "zuzax/zeroD/SurfPropagationSparta.h"

#include "surf_state.h"

/* ---------------------------------------------------------------------- */

SurfCollideZuzax::SurfCollideZuzax(SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal surf_collide zuzax command: id zuzax fileName twall acc");

  allowreact = 1;

  tstr = NULL;


  int n = strlen(arg[2]) + 1;
  inputConfigFile = new char[n];
  strcpy(inputConfigFile, arg[2]);
  
  if (strlen(inputConfigFile) < 2) {
      error->all(FLERR,"Expecting non-zero fileName");
  }
  FILE* fp = fopen(inputConfigFile, "r");
  if (!fp) {
     error->all(FLERR,"can't open inputConfigFile");
  }
  fclose(fp);

  twall = input->numeric(FLERR, arg[3]); 
  if (twall <= 0.0) {
    error->all(FLERR,"Surf_collide diffuse temp <= 0.0");
  }

  pwall = input->numeric(FLERR, arg[4]); 
  if (pwall <= 0.0) {
    error->all(FLERR,"Surf_collide diffuse press <= 0.0");
  }

  // Read the accomodation coefficient
  acc = input->numeric(FLERR,arg[5]); 
  if (acc < 0.0 || acc > 1.0) {
    error->all(FLERR,"Illegal surf_collide diffuse command");
  }
  vstream[0] = vstream[1] = vstream[2] = 0.0;

  // Initialize zuzax
  

  // initialize RNG

  random = new RanPark(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);
}

/* ---------------------------------------------------------------------- */

SurfCollideZuzax::~SurfCollideZuzax()
{
  if (copy) return;

  delete baseNet;
  delete [] tstr;
  delete random;
  delete inputConfigFile;
}

/* ---------------------------------------------------------------------- */

void SurfCollideZuzax::init()
{
  SurfCollide::init();

  // check variable

  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0) 
      error->all(FLERR,"Surf_collide diffuse variable name does not exist");
    if (!input->variable->equal_style(tvar))
      error->all(FLERR,"Surf_collide diffuse variable is invalid style");
  }

  initNetwork();
}
//==================================================================================================================================
void SurfCollideZuzax::initNetwork()
{
  Zuzax::XML_Node* xPhase =  Zuzax::get_XML_File(inputConfigFile);
  Zuzax::SolnMaterialDomain smd(*xPhase);
  Zuzax::PhaseList* pl = new Zuzax::PhaseList(smd);   

  Zuzax::IdealGasPhase* igp = dynamic_cast<Zuzax::IdealGasPhase*>(&(pl->volPhase(0)));
  if (!igp) {
    throw Zuzax::ZuzaxError("SurfCollideZuzax::initNetwork()", "First phase isn't an IdealGasPhase"); 
  }
  igp->setState_TP(twall, pwall);


  Zuzax::SubstrateElement* ablateFilm = new Zuzax::SubstrateElement("Ablator");


  /*
   *  Geometry of simulation (meters)
   */
  // This should be obtained from the geometry of the surface
  // TODO - connect
  double areaOfSurface = 1.0E-4;
  // The film depth will be needed eventually but is not significant now. We don't change the depth of the surface
  double filmDepth = 0.001;

  // Each bulk phase has a volume fraction associated with it. Here, we just assign the volume fraction
  // of the first phase to 1 and the rest to zero.
  double volumeFraction[10];
  volumeFraction[0] = 1.0;
  volumeFraction[1] = 0.0;
  volumeFraction[2] = 0.0;
  volumeFraction[3] = 0.0;


  
  size_t nVolPhases = pl->nVolPhases();
  for (size_t i = 1; i < nVolPhases; ++i) {
    double volPhase = areaOfSurface *filmDepth * volumeFraction[i];
    ablateFilm->addVolumePhase((pl->volPhase(i)), volPhase);
  }
  Zuzax::Reservoir* gasR = new Zuzax::Reservoir();
  gasR->setThermoMgr(*igp);


  ablateFilm->addLinkToAssociatedReactorVolume(gasR, false, true);
  ablateFilm->setState_TP(twall, pwall);

  // Create the object that will propagate the surface reactor forward in time
  /*
   *  This is built on top of the normal ODE solver that propagates the reactor in time using BDF methods
   */
  baseNet = new Zuzax::SurfPropagationSparta();
  // Add the surface reactor to the propagator
  /*
   * The surface reactor consists of the bulk volume phases and the interface
   */
  baseNet->addReactor(*ablateFilm, true);

  // Setup the interface between the reactor and the substrate where the reactions will occur
  // -> this will be a RBdry element which connects two reactors with an interfacial element
  Zuzax::RBdry* abFace_rrr = new Zuzax::RBdry("AirGraphite_rrr");

  // Define the area
  abFace_rrr->setArea(1.0E-7);

  // Turn on extended ROP analysis
  abFace_rrr->setKinBreakdownToggle(true);

  // Set up the InterfaceKinetics object that controls the surface reactions
  //   -> first we need a list of ThermoPhase objects
  std::vector<Zuzax::thermo_t_double*> thVec_abFace_rrr;
  thVec_abFace_rrr.push_back(igp);
  Zuzax::thermo_t_double* volPhase = ablateFilm->volPhasePtr(0);
  thVec_abFace_rrr.push_back(volPhase);
  Zuzax::InterfaceKinetics* abFaceKin_rrr = 
            (Zuzax::InterfaceKinetics*) Zuzax::newInterfaceKineticsMgrFromFile(thVec_abFace_rrr, "carbon_ablate.xml");

  Zuzax::thermo_t_double& sp = abFaceKin_rrr->reactionPhaseThermo();
  Zuzax::SurfPhase& surfPhase = dynamic_cast<Zuzax::SurfPhase&>(sp);
  abFace_rrr->setKinetics(abFaceKin_rrr);

  // Install the interface so that the equations reside in the SubstrateElement
  abFace_rrr->install(ablateFilm, gasR);

  double scl = 1;
  baseNet->setTMScaleFactor(scl);

  // The gas phase is added in as a reservoir
  baseNet->addReservoir(*gasR, true);

  // Add this so that abFace_rrr is deleted in the destructor
  baseNet->addUniqueConnector(abFace_rrr, true);

  size_t is_c = surfPhase.speciesIndex("(s)");
  size_t is_H = surfPhase.speciesIndex("H(s)");
  size_t is_O = surfPhase.speciesIndex("O(s)");
  size_t is_N = surfPhase.speciesIndex("N(s)");

  size_t ig_CO2 = igp->speciesIndex("CO2");
  size_t ig_CO = igp->speciesIndex("CO");

  double etchRate;
  double xmf[100];
  double theta[20];
  double press, uval, enth, vol, tt;
  FILE* FP = fopen("reactor_results.csv", "w");
  FILE* FP2 = fopen("step_results.csv", "w");
  baseNet->initialize();
  baseNet->printInitialSolution();

  baseNet->solveInitialConditions();

}
//==================================================================================================================================
/* ----------------------------------------------------------------------
   particle collision with surface with optional chemistry
   ip = particle with current x = collision pt, current v = incident v
   norm = surface normal unit vector
   isr = index of reaction model if >= 0, -1 for no chemistry
   ip = set to NULL if destroyed by chemsitry
   return jp = new particle if created by chemistry
   resets particle(s) to post-collision outward velocity
------------------------------------------------------------------------- */

Particle::OnePart *SurfCollideZuzax::
collide(Particle::OnePart *&ip, double *norm, double &, int isr, void* surfState)
{
  nsingle++;

  // if surface chemistry defined, attempt reaction
  // reaction = 1 if reaction took place

  Particle::OnePart iorig;
  Particle::OnePart *jp = NULL;
  int reaction = 0;

  if (isr >= 0) {
    // Save the original particle
    if (modify->n_surf_react) memcpy(&iorig,ip,sizeof(Particle::OnePart));
    // Call the surface reaction capability assigned to this surface
    reaction = surf->sr[isr]->react(ip,norm,jp);
    if (reaction) surf->nreact_one++;
  }

  // diffuse reflection for each particle
  // resets v, roteng, vibeng
  // if new particle J created, also need to trigger any fixes

  if (ip) diffuse(ip,norm);
  if (jp) {
    diffuse(jp,norm);
    if (modify->n_add_particle) {
      int j = jp - particle->particles;
      modify->add_particle(j,twall,twall,twall,vstream);
    }
  }

  // call any fixes with a surf_react() method
  // they may reset j to -1, e.g. fix ambipolar
  //   in which case newly created j is deleted

  if (reaction && modify->n_surf_react) {
    int i = -1;
    if (ip) i = ip - particle->particles;
    int j = -1;
    if (jp) j = jp - particle->particles;
    modify->surf_react(&iorig,i,j);
    if (jp && j < 0) {
      jp = NULL;
      particle->nlocal--;
    }
  }

  return jp;
}

/* ----------------------------------------------------------------------
   diffusive particle collision with surface
   p = particle with current x = collision pt, current v = incident v
   norm = surface normal unit vector
   resets particle(s) to post-collision outward velocity
------------------------------------------------------------------------- */

void SurfCollideZuzax::diffuse(Particle::OnePart *p, double *norm)
{
  // specular reflection
  // reflect incident v around norm

  if (random->uniform() > acc) {
    MathExtra::reflect3(p->v,norm);
    p->erot = particle->erot(p->ispecies,twall,random);
    p->evib = particle->evib(p->ispecies,twall,random);

  // diffuse reflection
  // vrm = most probable speed of species, eqns (4.1) and (4.7)
  // vperp = velocity component perpendicular to surface along norm, eqn (12.3)
  // vtan12 = 2 velocity components tangential to surface
  // tangent1 = component of particle v tangential to surface,
  //   check if tangent1 = 0 (normal collision), set randomly
  // tangent2 = norm x tangent1 = orthogonal tangential direction
  // tangent12 are both unit vectors

  } else {
    double tangent1[3],tangent2[3];
    Particle::Species *species = particle->species;
    int ispecies = p->ispecies;

    double vrm = sqrt(2.0*update->boltz * twall / species[ispecies].mass);
    double vperp = vrm * sqrt(-log(random->uniform()));

    double theta = MY_2PI * random->uniform();
    double vtangent = vrm * sqrt(-log(random->uniform()));
    double vtan1 = vtangent * sin(theta);
    double vtan2 = vtangent * cos(theta);

    double *v = p->v;
    double dot = MathExtra::dot3(v,norm);

    double beta_un,normalized_distbn_fn;

    tangent1[0] = v[0] - dot*norm[0];
    tangent1[1] = v[1] - dot*norm[1];
    tangent1[2] = v[2] - dot*norm[2];

    if (MathExtra::lensq3(tangent1) == 0.0) {
      tangent2[0] = random->uniform();
      tangent2[1] = random->uniform();
      tangent2[2] = random->uniform();
      MathExtra::cross3(norm,tangent2,tangent1);
    }

    MathExtra::norm3(tangent1);
    MathExtra::cross3(norm,tangent1,tangent2);

    // add in translation or rotation vector if specified
    // only keep portion of vector tangential to surface element

    if (trflag) {
      double vxdelta,vydelta,vzdelta;
      if (tflag) {
        vxdelta = vx; vydelta = vy; vzdelta = vz;
        double dot = vxdelta*norm[0] + vydelta*norm[1] + vzdelta*norm[2];
     
        if (fabs(dot) > 0.001) {
          dot /= vrm;
          do {
            do {
              beta_un = (6.0*random->uniform() - 3.0);
            } while (beta_un + dot < 0.0);
            normalized_distbn_fn = 2.0 * (beta_un + dot) /
              (dot + sqrt(dot*dot + 2.0)) *
              exp(0.5 + (0.5*dot)*(dot-sqrt(dot*dot + 2.0)) -
                  beta_un*beta_un);
          } while (normalized_distbn_fn < random->uniform());
          vperp = beta_un*vrm;
        }

      } else {
        double *x = p->x;
        vxdelta = wy*(x[2]-pz) - wz*(x[1]-py);
        vydelta = wz*(x[0]-px) - wx*(x[2]-pz);
        vzdelta = wx*(x[1]-py) - wy*(x[0]-px);
        double dot = vxdelta*norm[0] + vydelta*norm[1] + vzdelta*norm[2];
        vxdelta -= dot*norm[0];
        vydelta -= dot*norm[1];
        vzdelta -= dot*norm[2];
      }
      
      v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0] + vxdelta;
      v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1] + vydelta;
      v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2] + vzdelta;

    // no translation or rotation

    } else {
      v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0];
      v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1];
      v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2];
    }

    // initialize rot/vib energy

    p->erot = particle->erot(ispecies,twall,random);
    p->evib = particle->evib(ispecies,twall,random);
  }
}

/* ----------------------------------------------------------------------
   set current surface temperature
------------------------------------------------------------------------- */

void SurfCollideZuzax::dynamic()
{
  twall = input->variable->compute_equal(tvar);
  if (twall <= 0.0) error->all(FLERR,"Surf_collide diffuse temp <= 0.0");
}

/* ----------------------------------------------------------------------
    Provides a surface state object
------------------------------------------------------------------------- */

SurfState* SurfCollideZuzax::provideStateObject() const
{
    SurfState *ss = new SurfState();
    ss->net = new Zuzax::SurfPropagationSparta(*baseNet);
    return ss; 
}

#endif
