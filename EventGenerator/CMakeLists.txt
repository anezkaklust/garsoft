
art_make(MODULE_LIBRARIES Geometry
                           Geometry_Geometry_service
                           SummaryDataProducts
                           nusimdata_SimulationBase
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                           art_Persistency_Common     canvas_Persistency_Common
                           art_Persistency_Provenance canvas_Persistency_Provenance
                           art_Utilities              canvas_Utilities
                           ${MF_MESSAGELOGGER}
                           ${MF_UTILITIES}
                           ${FHICLCPP}
                           ${CETLIB}
                           ${CLHEP}
                           ${ROOT_BASIC_LIB_LIST}
                           ${ROOT_EG}
         EXCLUDE GeneratedEventTimestamp_plugin.cc
        )

simple_plugin( GeneratedEventTimestamp "plugin"
               ${ART_FRAMEWORK_CORE}
               art_Persistency_Provenance canvas_Persistency_Provenance
               ${MF_MESSAGELOGGER}
               ${MF_UTILITIES}
               ${CETLIB}
               ${Boost_FILESYSTEM_LIBRARY}
				       ${Boost_SYSTEM_LIBRARY}
             )

install_headers()
install_fhicl()
install_source()

add_subdirectory(CRY)
add_subdirectory(GENIE)
