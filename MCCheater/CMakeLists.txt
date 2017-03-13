
include_directories( $ENV{NUTOOLS_DIR}/include )

art_make_library( LIBRARY_NAME MCCheater
                  SOURCE       BackTrackerCore.cxx
                  LIBRARIES    ReconstructionDataProducts
                               SimulationDataProducts
                               nusimdata_SimulationBase
                               nutools_ParticleNavigation
                               art_Persistency_Common     canvas_Persistency_Common
                               art_Persistency_Provenance canvas_Persistency_Provenance
                               art_Utilities canvas_Utilities
                               ${MF_MESSAGELOGGER}
                               ${MF_UTILITIES}
                               ${FHICLCPP}
                               ${CETLIB}
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
			        art_Persistency_Common     canvas_Persistency_Common
			        art_Persistency_Provenance canvas_Persistency_Provenance
			        art_Utilities              canvas_Utilities
			        ${MF_MESSAGELOGGER}
			        ${MF_UTILITIES}
			        ${FHICLCPP}
			        ${CETLIB}
			        ${ROOT_BASIC_LIB_LIST}
)

simple_plugin(BackTrackerLoader "module"
              MCCheater_BackTracker_service
        			${ART_FRAMEWORK_SERVICES_REGISTRY}              
              
)

install_headers()
install_fhicl()
install_source()
