
include_directories( $ENV{NUTOOLS_DIR}/include )

art_make_library( LIBRARY_NAME MCCheater
                  SOURCE       BackTrackerCore.cxx
                  LIBRARIES    ReconstructionDataProducts
                               SimulationDataProducts
                               Geometry
                               Geometry_Geometry_service
                               nusimdata_SimulationBase
                               nutools_ParticleNavigation
                               ${ART_FRAMEWORK_SERVICES_REGISTRY}
                               art_Persistency_Common
                               art_Persistency_Provenance
                               art_Utilities
                               canvas
                               ${MF_MESSAGELOGGER}
                               
                               ${FHICLCPP}
                               ${CETLIB_EXCEPT}
                               ${Boost_SYSTEM_LIBRARY}
                )

simple_plugin(BackTracker "service"
              MCCheater
	          Geometry
	          Geometry_Geometry_service
              ReconstructionDataProducts
              RawDataProducts
              SimulationDataProducts
              ${ART_FRAMEWORK_CORE}
			        ${ART_FRAMEWORK_PRINCIPAL}
			        ${ART_FRAMEWORK_SERVICES_REGISTRY}
			        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			        art_Persistency_Common
			        art_Persistency_Provenance
			        art_Utilities
              canvas
			        ${MF_MESSAGELOGGER}			        
			        ${FHICLCPP}
			        ${CETLIB_EXCEPT}
			        ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()
