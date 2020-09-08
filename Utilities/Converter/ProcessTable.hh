// File use to convert the process and sub process number into a process name
// The numbers can be checked in the include files of G4

// /cvmfs/larsoft.opensciencegrid.org/products/geant4/v4_10_3_p03e/Linux64bit+3.10-2.17-e19-prof/include/Geant4/G4DecayProcessType.hh
// /cvmfs/larsoft.opensciencegrid.org/products/geant4/v4_10_3_p03e/Linux64bit+3.10-2.17-e19-prof/include/Geant4/G4EmProcessSubType.hh
// /cvmfs/larsoft.opensciencegrid.org/products/geant4/v4_10_3_p03e/Linux64bit+3.10-2.17-e19-prof/include/Geant4/G4FastSimulationProcessType.hh
// /cvmfs/larsoft.opensciencegrid.org/products/geant4/v4_10_3_p03e/Linux64bit+3.10-2.17-e19-prof/include/Geant4/G4HadronicProcessType.hh
// /cvmfs/larsoft.opensciencegrid.org/products/geant4/v4_10_3_p03e/Linux64bit+3.10-2.17-e19-prof/include/Geant4/G4OpProcessSubType.hh
// /cvmfs/larsoft.opensciencegrid.org/products/geant4/v4_10_3_p03e/Linux64bit+3.10-2.17-e19-prof/include/Geant4/G4ProcessType.hh
// /cvmfs/larsoft.opensciencegrid.org/products/geant4/v4_10_3_p03e/Linux64bit+3.10-2.17-e19-prof/include/Geant4/G4TransportationProcessType.hh
// /cvmfs/larsoft.opensciencegrid.org/products/geant4/v4_10_3_p03e/Linux64bit+3.10-2.17-e19-prof/include/Geant4/G4UCNProcessSubType.hh

// The convention is Process:SubProcess such as 2:14 -> EM:GammaConversion or 4:111 -> HAD:HadronElastic

#include <unordered_map>
#include <string>

namespace gar {
    namespace util {

        inline std::unordered_map<int, std::unordered_map<int, std::string> > &processTable()
        {
            static std::unordered_map<int, std::unordered_map<int, std::string> > m_processTable;
            static bool firstTime = true;
            if(firstTime)
            {
                firstTime = false;
                //Not defined

                //Transportation
                m_processTable[1][91] = "Transportation:Transportation";
                m_processTable[1][92] = "Transportation:CoupledTransportation";

                //Electromagnetic
                m_processTable[2][1] = "EM:CoulombScattering";
                m_processTable[2][2] = "EM:Ionisation";
                m_processTable[2][3] = "EM:Bremsstrahlung";
                m_processTable[2][4] = "EM:PairProductionByCharged";
                m_processTable[2][5] = "EM:Annihilation";
                m_processTable[2][6] = "EM:AnnihilationToMuMu";
                m_processTable[2][7] = "EM:AnnihilationToHadrons";
                m_processTable[2][8] = "EM:NuclearStopping";
                m_processTable[2][10] = "EM:MultipleScattering";
                m_processTable[2][11] = "EM:Rayleigh";
                m_processTable[2][12] = "EM:PhotoElectricEffect";
                m_processTable[2][13] = "EM:ComptonScattering";
                m_processTable[2][14] = "EM:GammaConversion";
                m_processTable[2][15] = "EM:GammaConversionToMuMu";
                m_processTable[2][21] = "EM:Cerenkov";
                m_processTable[2][22] = "EM:Scintillation";
                m_processTable[2][23] = "EM:SynchrotronRadiation";
                m_processTable[2][24] = "EM:TransitionRadiation";

                //Optical
                m_processTable[4][31] = "OPT:OpAbsorption";
                m_processTable[4][32] = "OPT:OpBoundary";
                m_processTable[4][33] = "OPT:OpRayleigh";
                m_processTable[4][34] = "OPT:OpWLS";
                m_processTable[4][35] = "OPT:OpMieHG";

                //Hadronic
                m_processTable[4][111] = "HAD:HadronElastic";
                m_processTable[4][121] = "HAD:HadronInelastic";
                m_processTable[4][131] = "HAD:Capture";
                m_processTable[4][141] = "HAD:Fission";
                m_processTable[4][151] = "HAD:HadronAtRest";
                m_processTable[4][152] = "HAD:LeptonAtRest";
                m_processTable[4][161] = "HAD:ChargeExchange";
                m_processTable[4][210] = "HAD:RadioactiveDecay";

                //PhotoLepton_Hadron

                //Decay
                m_processTable[6][201] = "DEC:Decay";
                m_processTable[6][202] = "DEC:DecayWithSpin";
                m_processTable[6][203] = "DEC:DecayPionMakeSpin";
                m_processTable[6][210] = "DEC:Radioactive";
                m_processTable[6][211] = "DEC:Unknown";
                m_processTable[6][231] = "DEC:External";

                //General
                m_processTable[7][401] = "General:StepLimiter";
                m_processTable[7][402] = "General:UserSpecialCuts";
                m_processTable[7][403] = "General:NeutronKiller";

                //Parametrisation

                //UserDefined

                //Parallel

                //Photon

                //UCN
                m_processTable[12][41] = "UCN:UCNLoss";
                m_processTable[12][42] = "UCN:UCNAbsorption";
                m_processTable[12][43] = "UCN:UCNBoundary";
                m_processTable[12][44] = "UCN:UCNMultiScattering";
            }
            return m_processTable;
        }

        inline std::string FindProcessName(int process, int subprocess);
    }
}

inline std::string gar::util::FindProcessName(int process, int subprocess)
{
    std::string process_name = "unknown";
    if(processTable().find(process) != processTable().end()) {
        if(processTable()[process].find(subprocess) != processTable()[process].end()) {
            process_name = processTable()[process][subprocess];
        }
    }
    return process_name;
}
