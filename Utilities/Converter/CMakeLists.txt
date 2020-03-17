add_subdirectory(edep-io)

include_directories( $ENV{CLHEP_INCLUDE_DIR}     )

art_make(LIB_LIBRARIES      Utilities_Converter_edep-io
                            SimulationDataProducts
                            Geometry
                            ${ART_FRAMEWORK_PRINCIPAL}
                            ${ART_FRAMEWORK_SERVICES_REGISTRY}
                            ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                            ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                            ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                            art_Persistency_Common
                            art_Persistency_Provenance
                            art_Utilities
                            canvas
                            ${MF_MESSAGELOGGER}
                            ${FHICLCPP}
                            ${CETLIB}
                            ${CETLIB_EXCEPT}
                            ${ROOT_EG}
                            ${ROOT_TREEPLAYER}
                            ${ROOT_FFTW}
                            ${ROOT_REFLEX}
                            ${ROOT_EGPYTHIA6}
                            ${ROOT_GUI}
                            ${ROOT_BASIC_LIB_LIST}
                            ${ROOT_GEOM}
                            ${CLHEP}
                            ${Boost_SYSTEM_LIBRARY}
       MODULE_LIBRARIES     Utilities_Converter_edep-io
                            SimulationDataProducts
                            Geometry
                            MCCheater_BackTracker_service
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
                            ${CETLIB}
                            ${CETLIB_EXCEPT}
                            ${ROOT_BASIC_LIB_LIST}
)


install_headers()
install_fhicl()
install_source()
