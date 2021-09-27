include_directories( $ENV{GEANT4_FQ_DIR}/include )
include_directories( $ENV{NUTOOLS_DIR}/include   )
include_directories( $ENV{CLHEP_INCLUDE_DIR}     )
include_directories( $ENV{GSL_FQ_DIR}/include     )

cet_find_library(GSGSL NAMES gsl PATHS ENV GSL_FQ_DIR PATH_SUFFIXES lib NO_DEFAULT_PATH)
cet_find_library(GSGSLCBLAS NAMES gslcblas PATHS ENV GSL_FQ_DIR PATH_SUFFIXES lib NO_DEFAULT_PATH)


cet_make( LIBRARIES
                       Geometry
                       ${ART_FRAMEWORK_PRINCIPAL}
                       ${ART_FRAMEWORK_SERVICES_REGISTRY}
                       ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                       ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                       ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                       art_Persistency_Common
                       art_Persistency_Provenance
                       art_Utilities
	                    canvas
                       ${MF_MESSAGELOGGER}
                       ${FHICLCPP}
                       ${CETLIB_EXCEPT}
                       ${ROOT_BASIC_LIB_LIST}
                       ${CLHEP}
                       ${Boost_SYSTEM_LIBRARY}
		       GSGSL
		       GSGSLCBLAS
                       nug4_MagneticFieldServices_MagneticFieldServiceStandard_service
                       NO_DICTIONARY )

# Simulation base is included in making the dictionary for associations
# between the MC truth and reconstructed quantities when cheating

art_dictionary( DICTIONARY_LIBRARIES ReconstructionDataProducts
                                     RawDataProducts
                                     nug4_MagneticFieldServices_MagneticFieldServiceStandard_service
                                     nusimdata_SimulationBase)

install_headers()
install_source()
