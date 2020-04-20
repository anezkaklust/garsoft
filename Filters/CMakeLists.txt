include_directories( $ENV{NUTOOLS_DIR}/include   )

art_make(LIB_LIBRARIES
           nusimdata_SimulationBase
           Geometry
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           canvas
           ${MF_MESSAGELOGGER}
           ${FHICLCPP}
           cetlib_except
           ${ART_FRAMEWORK_PRINCIPAL}
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
           art_Persistency_Common
           art_Persistency_Provenance
           art_Utilities
           ${ROOT_BASIC_LIB_LIST}
         MODULE_LIBRARIES
           nusimdata_SimulationBase
           ${ART_FRAMEWORK_CORE}
           ${ART_FRAMEWORK_IO_SOURCES}
           ${ART_FRAMEWORK_PRINCIPAL}
           art_Persistency_Common
           art_Persistency_Provenance
           art_Utilities
           ${MF_MESSAGELOGGER}
           ${ROOT_BASIC_LIB_LIST}
           Geometry
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
         )

install_headers()
install_fhicl()
install_source()
