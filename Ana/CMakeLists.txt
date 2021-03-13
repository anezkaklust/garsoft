include_directories( $ENV{GEANT4_FQ_DIR}/include )
include_directories( $ENV{NUTOOLS_DIR}/include   )
include_directories( $ENV{CLHEP_INCLUDE_DIR}     )

art_make(BASENAME_ONLY
      MODULE_LIBRARIES Geometry
                       ReadoutSimulation
                       SimulationDataProducts
                       RawDataProducts
                       ReconstructionDataProducts
                       AnalysisDataProducts
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
                       ${ROOT_BASIC_LIB_LIST}
                       ${ROOT_EG}
                       GaranaProducts
         )

install_headers()
install_fhicl()
install_source()

add_subdirectory(ParamSim)
add_subdirectory(DSTproduction)
add_subdirectory(gntpc_dune)
add_subdirectory(TPCperformance)
add_subdirectory(HeavyNeutralLeptonAnalysis)
