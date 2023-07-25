include_directories( $ENV{GEANT4_FQ_DIR}/include )
include_directories( $ENV{NUTOOLS_DIR}/include   )
include_directories( $ENV{CLHEP_INCLUDE_DIR}     )
link_directories( $ENV{NUG4_FQ_DIR}/lib )

art_make(
	BASENAME_ONLY
	LIB_LIBRARIES garana::GaranaProducts
		      ReconstructionDataProducts
		      nusimdata_SimulationBase
                    ROOT::Core
                    ROOT::Hist
                    ROOT::Tree              
        MODULE_LIBRARIES Geometry
                       ReadoutSimulation
                       SimulationDataProducts
                       RawDataProducts
                       ReconstructionDataProducts
                       MCCheater_BackTracker_service
                       nusimdata_SimulationBase
                       nurandom_RandomUtils_NuRandomService_service
                       nug4_MagneticFieldServices_MagneticFieldServiceStandard_service
                       nugen_EventGeneratorBase_GENIE
                       ${ART_FRAMEWORK_CORE}
                       ${ART_FRAMEWORK_IO_SOURCES}
                       ${ART_FRAMEWORK_PRINCIPAL}
                       art_Persistency_Common
                       art_Persistency_Provenance
                       art_Utilities
                       canvas
                       ${ART_FRAMEWORK_SERVICES_REGISTRY}
                       ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                       ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                       ${MF_MESSAGELOGGER}
                       ${FHICLCPP}
                       ${CETLIB_EXCEPT}
                    ROOT::Core
                    ROOT::Hist
                    ROOT::Tree              
                       ${ROOT_EG}
                       garana::GaranaProducts # use when setting up GArAna from UPS (default behavior)
                       garsoft_StandardRecord
         )

install_headers()
install_fhicl()
install_source()

add_subdirectory(ParamSim)
add_subdirectory(gntpc_dune)
add_subdirectory(TPCperformance)
add_subdirectory(HeavyNeutralLeptonAnalysis)
