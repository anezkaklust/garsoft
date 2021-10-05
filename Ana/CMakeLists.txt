include_directories( $ENV{GEANT4_FQ_DIR}/include )
include_directories( $ENV{NUTOOLS_DIR}/include   )
include_directories( $ENV{GARANA_DIR}/include )
include_directories( $ENV{CLHEP_INCLUDE_DIR}     )

# N.B. due to some issues setting up garana libraries, there are two ways to 
# to set that library paths, depending on whether or not you are developing
# side-by-side with GArAna or not (two lines below)

art_make(
	BASENAME_ONLY
	LIB_LIBRARIES #$ENV{GARANA_LIB}/libGaranaProducts.so # use when developing side-by-side with GArAna
                      GaranaProducts # use when setting up GArAna from UPS (default behavior)
		      ReconstructionDataProducts
		      nusimdata_SimulationBase
		      ROOT_BASIC_LIB_LIST
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
                       ROOT_BASIC_LIB_LIST
                       ${ROOT_EG}
                       #$ENV{GARANA_LIB}/libGaranaProducts.so # use when developing side-by-side with GArAna
                       GaranaProducts # use when setting up GArAna from UPS (default behavior)
                       garsoft_StandardRecord
         )

install_headers()
install_fhicl()
install_source()

add_subdirectory(ParamSim)
add_subdirectory(DSTproduction)
add_subdirectory(gntpc_dune)
add_subdirectory(TPCperformance)
add_subdirectory(HeavyNeutralLeptonAnalysis)
