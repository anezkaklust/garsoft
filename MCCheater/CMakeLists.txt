
include_directories( $ENV{NUTOOLS_DIR}/include )

art_make_library( LIBRARY_NAME MCCheater
                  SOURCE       BackTrackerCore.cxx
                  LIBRARIES    ReconstructionDataProducts
                               SimulationDataProducts
                               Geometry
                               Geometry_GeometryGAr_service
                               nusimdata_SimulationBase
                               nug4_ParticleNavigation
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
	          Geometry_GeometryGAr_service
              ReconstructionDataProducts
              RawDataProducts
              SimulationDataProducts
              ${ART_FRAMEWORK_CORE}
			        ${ART_FRAMEWORK_PRINCIPAL}
			        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                                ${ART_ROOT_IO_TFILESERVICE_SERVICE}
			        art_Persistency_Common
			        art_Persistency_Provenance
			        art_Utilities
                                canvas
			        ${MF_MESSAGELOGGER}			        
			        ${FHICLCPP}
			        ${CETLIB_EXCEPT}
                                ROOT::Core
                                ROOT::Hist
                                ROOT::Tree              
)

install_headers()
install_fhicl()
install_source()
