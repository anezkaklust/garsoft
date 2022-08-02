include_directories( $ENV{NUGEN_INC} )

art_make(MODULE_LIBRARIES Geometry
                           Geometry_GeometryGAr_service
                           SummaryDataProducts
                           nusimdata_SimulationBase
                           nurandom_RandomUtils_NuRandomService_service
                           ${ART_FRAMEWORK_CORE}
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
                           ${CETLIB}
                           ${CETLIB_EXCEPT}
                           ${CLHEP}
                           ROOT::Core
                           ROOT::Hist
                           ROOT::Tree              
                           ROOT::EG
        )

install_headers()
install_fhicl()
install_source()

add_subdirectory(CRY)
add_subdirectory(GENIE)
