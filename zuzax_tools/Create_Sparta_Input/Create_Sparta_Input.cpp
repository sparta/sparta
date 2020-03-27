
#include <stdio.h>

#include "zuzax/thermo.h"
#include "zuzax/thermo/StatMech.h"




using namespace std;
using namespace Zuzax;

int CHECK_DEBUG_MODE = 0;

void printUsage()
{
    cout << "usage: HMW_test_1 " <<  endl;
    cout <<"                -> Everything is hardwired" << endl;
}
//==================================================================================================================================
void printline(const char* c, int n)
{
    for (int i = 0; i < n; i++) {
        printf("%s", c);
    }
    printf("\n");
}

//==================================================================================================================================
int main(int argc, char** argv)
{

    int retn = 0;

    try {
        double x[20];
        for (int k = 0; k < 20; k++) {
            x[k] = 0.0;
        }

        std::string zuzaxFileName = "airSM.xml";
        std::string zuzaxCompFileName = "air.xml";
#ifdef THERMOPHASE_TEMPLATED
        ThermoPhase* airSM = newPhase<doublevalue>(0.0, zuzaxFileName.c_str(), "");
#else
        ThermoPhase* airSM = newPhase(zuzaxFileName.c_str(), "");
#endif

#ifdef THERMOPHASE_TEMPLATED
        ThermoPhase* air = newPhase<doublevalue>(0.0, zuzaxCompFileName.c_str(), "");
#else
        ThermoPhase* air = newPhase(zuzaxCompFileName.c_str(), "");
#endif


        double enthl_RT[50], cp_R[50], s_R[50];
        double enthl_RT_SM[50], cp_R_SM[50], s_R_SM[50];

        size_t indexN_SM =airSM->speciesIndex("N2");
        size_t indexN = air->speciesIndex("N2");

        /*
         * set states
         */
        x[0] = 0.7;
        x[1] = 1.0 - x[0];
        double T = 273.15 + 352.;
        air->setState_TPX(T, OneBar, x);
        airSM->setState_TPX(T, OneBar, x);

        std::vector<std::string> speciesToInclude;
        speciesToInclude.push_back("N2");
        speciesToInclude.push_back("O2");
        for (size_t k = 0; k < airSM->nSpecies(); ++k) {
           std::string ss = airSM->speciesName(k);
           speciesToInclude.push_back(ss);
        }

        // Get the Stat Mech Struct. We'll need to reference it directly
        SpeciesThermo& st = airSM->speciesThermo();

        //
        //                    Write Header for air.species file to be written
        //
        FILE* FP = fopen("air.species", "w");

        // Dump N2 thermo to Sparta input file
        fprintf(FP,"#   File create by Create_Sparta_Input \n");
        fprintf(FP,"#   Original Zuzax file: %s\n", zuzaxFileName.c_str());
        fprintf(FP,"# \n");
        fprintf(FP,"#   SpeciesName      AMU           MW       RotDOF 1/RotRelax#      "
                "vibDOF     1/VibRelax#     VibTemp     speciesW     Charge \n");

        fprintf(FP,"# \n"); 
        //
        // Loop Over all species that are to be output to the Sparta input file
        //
        for  (size_t k = 0; k < speciesToInclude.size(); ++k) {
           std::string speciesName = speciesToInclude[k];
          
           indexN = airSM->speciesIndex(speciesName);
           size_t iSp_r = air->speciesIndex(speciesName);
           double hf298 = 0.0; 
           if (iSp_r != npos) {
              hf298 = air->Hf298SS(iSp_r);   
              printf("Hf298 ( %s ) (NASA) = %g J/kmol\n", speciesName.c_str(), hf298);
           }
           if (indexN == npos) {
               throw ZuzaxError("Create_Sparta_Input Error",
                                "Species name, %s, not found in ThermoPhase object", speciesName.c_str());
           }

           // Go get the pointer to the SpeciesThermoInterpType object for the current species
           const SpeciesThermoInterpType* stit = st.provideTempDepSTIT(indexN,  298.15);
           // Make sure it's a StatMech type -> won't work if it isn't
           const StatMech* gsm = dynamic_cast<const StatMech*>(stit);
           if (!gsm) {
               throw ZuzaxError("main", "SpeciesThermo for species, %s, isn't a StatMech object", speciesName.c_str());
           }

           // Get the StatMech Input file for the species
           const StatMech::speciesStatMechInput& si = gsm->statMechInput();

           if (si.theta.size() > 0) {
           double theta0 = si.theta[0]; 
           double freq0 = Boltzmann * theta0 / (lightSpeed * 100. * Planck);


           printf("theta = %g K\n", theta0);
           printf("freq0 = %g cm-1\n", freq0);
           }
           double rotRelN =  si.rotRelaxationNumber;

           // Write out the species name
           fprintf(FP, "%16s ", speciesName.c_str());

           // Write out the atomic mass units of the species
           //  -> Here we assume that the molecular weight is the average of the isotope compositions
           double moleW = airSM->molecularWeight(indexN);
           fprintf(FP," % 15.8E ", moleW);

           // Write out the mass of a single molecule
           //  -> Here we assume that the molecular weight is the average of the isotope compositions
           double moleW_molecule = moleW/Zuzax::Avogadro;
           fprintf(FP," % 15.8E ", moleW_molecule);

           // Write out the number of rotational degrees of freedom
           int ndofRot = 0;
           if (si.geom == "linear") {
               ndofRot = 2;
           } else if (si.geom == "nonlinear") {
               ndofRot = 3;
           } else if (si.geom == "atom") {
               ndofRot = 0;
           } else {
               throw ZuzaxError("main", "Unknown geom keyword: %s", si.geom.c_str());
           }
           fprintf(FP," % 5d   ", ndofRot); 

           // Write out the inverse of the rotational relaxation rate
           // for single atoms 
           if (si.geom == "atom") {
               fprintf(FP, "  % 9.3E  " , 0.0);
           } else {
               fprintf(FP, "  % 9.3E  " , 1.0 / rotRelN);
           }
           // write out the vibrational DOFS x 2 
           int nn = 2 * si.nvib;
           fprintf(FP, "  % 5d  " , nn);
 
           double vibRelN =  si.vibRelaxationNumber;
           if (si.nvib > 0) {
               fprintf(FP, "  % 12.3E  " , 1.0 / vibRelN);
           } else {
               fprintf(FP, "  % 12.6E  " , 0.0);
           }

           // Print out the first vibrational temperature (doesn't seem to be room for multiple temperature)
           if (si.nvib > 0) {
               fprintf(FP, " % 12.6E  ", si.theta[0]);
           } else {
               fprintf(FP, "  % 12.6E  " , 0.0);
           }

           // Print out the species weighting factor = 1.0
           fprintf(FP, " 1.0   ");

           // Print out the charge of the species
           double ch = airSM->charge(indexN);
           fprintf(FP, "  % 12.6E  " , ch);

           fprintf(FP, "\n");

        }

        //---------------------------------------------------------------------------------------------------------------
        double Hf298_N_air = air->Hf298SS(indexN);
        double Hf298_N_airSM = airSM->Hf298SS(indexN_SM);
 
        printf("   Comparison of N2 Atom Heat of Formation (GRI vs StatMech):\n");
        printf(" Hf298:   | % 24.15E % 24.15E |\n", Hf298_N_air, Hf298_N_airSM);

        //---------------------------------------------------------------------------------------------------------------




        printline("-", 80);
        printf("    Comparison of Enthalpy Calculation of N2 Atom (GRI vs StatMech):\n");
        printline("-", 80);
        printf("  Temp    | Enth_reg(\"N2\") Enth_SM(\"N2\")|    DIFF     |   DIFF/RT   |\n");
        printline("-", 80);
        for (size_t i = 0; i < 20; i++) {
            if (i == 0) {
                T = 298.15;
            } else if (i == 1) {
                T = 300.0;
            } else {
                T += 100.;
            }
            air->setState_TPX(T, OneBar, x);
            airSM->setState_TPX(T, OneBar, x);
            airSM->getEnthalpy_RT(enthl_RT_SM);
            air->getEnthalpy_RT(enthl_RT);
            double RT = GasConstant * T;
            printf(" %7g  |", T);
            printf(" %12.5E %12.5E | ", enthl_RT[indexN] * RT, enthl_RT_SM[indexN_SM] * RT);
            printf(" % 10.3E | ", (enthl_RT_SM[indexN_SM] - enthl_RT[indexN]) * RT);
            printf(" % 10.3E | ", (enthl_RT_SM[indexN_SM] - enthl_RT[indexN]));
            printf("\n");
        }
        printline("-", 80);

        //---------------------------------------------------------------------------------------------------------------
        printf("    Comparison of Cp Calculation of O Atom (GRI vs StatMech):\n");
        printline("-", 80);
        printf("  Temp    |  Cp_reg(\"N\")  Cp_SM(\"N\")  |    DIFF     |    DIFF/R   |\n");
        printline("-", 80);
        for (size_t i = 0; i < 20; i++) {
            if (i == 0) {
                T = 298.15;
            } else if (i == 1) {
                T = 300.0;
            } else {
                T += 100.;
            }
            air->setState_TPX(T, OneBar, x);
            airSM->setState_TPX(T, OneBar, x);
            airSM->getCp_R(cp_R_SM);
            air->getCp_R(cp_R);
            printf(" %7g  |", T);
            printf(" %12.5E %12.5E | ", cp_R[indexN] * GasConstant, cp_R_SM[indexN_SM] * GasConstant);
            printf(" % 10.3E | ", (cp_R_SM[indexN_SM] - cp_R[indexN]) * GasConstant);
            printf(" % 10.3E | ", (cp_R_SM[indexN_SM] - cp_R[indexN]));
            printf("\n");
        }
        printline("-", 80);

        //---------------------------------------------------------------------------------------------------------------
        printf("    Comparison of Entropy Calculation of N2 Atom (GRI vs StatMech):\n");
        printline("-", 80);
        printf("  Temp    |   S_reg(\"N2\")   S_SM(\"N2\")  |    DIFF     |    DIFF/R   |\n");
        printline("-", 80);
        for (size_t i = 0; i < 20; i++) {
            if (i == 0) {
                T = 298.15;
            } else if (i == 1) {
                T = 300.0;
            } else {
                T += 100.;
            }
            air->setState_TPX(T, OneBar, x);
            airSM->setState_TPX(T, OneBar, x);
            airSM->getEntropy_R(s_R_SM);
            air->getEntropy_R(s_R);
            printf(" %7g  |", T);
            printf(" %12.5E %12.5E | ", s_R[indexN] * GasConstant, s_R_SM[indexN_SM] * GasConstant);
            printf(" % 10.3E | ", (s_R_SM[indexN_SM] - s_R[indexN]) * GasConstant);
            printf(" % 10.3E | ", (s_R_SM[indexN_SM] - s_R[indexN]));
            printf("\n");
        }
        printline("-", 80);
        //---------------------------------------------------------------------------------------------------------------
        printf("    Comparison of IntEng Calculation of N2 Atom (GRI vs StatMech):\n");
        printline("-", 80);
        printf("  Temp    | IntE_reg(\"N2\")  IntE_SM(\"N2\") |    DIFF     |    DIFF/RT  |\n");
        printline("-", 80);
        for (size_t i = 0; i < 20; i++) {
            if (i == 0) {
                T = 298.15;
            } else if (i == 1) {
                T = 300.0;
            } else {
                T += 100.;
            }
            air->setState_TPX(T, OneBar, x);
            airSM->setState_TPX(T, OneBar, x);
            airSM->getIntEnergy_RT(s_R_SM);
            air->getIntEnergy_RT(s_R);
            printf(" %7g  |", T);
            printf(" %12.5E %12.5E | ", s_R[indexN] * GasConstant * T, s_R_SM[indexN_SM] * GasConstant * T);
            printf(" % 10.3E | ", (s_R_SM[indexN_SM] - s_R[indexN]) * GasConstant * T);
            printf(" % 10.3E | ", (s_R_SM[indexN_SM] - s_R[indexN]));
            printf("\n");
        }
        printline("-", 80);
        //---------------------------------------------------------------------------------------------------------------

        delete air;
        delete airSM;

        appdelete();

        return retn;

    } catch (ZuzaxError) {

        showErrors();
        return -1;
    }
}
