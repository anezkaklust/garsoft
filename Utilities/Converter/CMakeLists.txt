add_subdirectory(edep-io)

include_directories( $ENV{CLHEP_INCLUDE_DIR} )
include_directories( $ENV{GENIE_INC}/GENIE )

execute_process(COMMAND genie-config --libs OUTPUT_VARIABLE GENIE_LIBS_DIR_WITH_LFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)

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
                            nutools_EventGeneratorBase_GENIE
                            nutools_RandomUtils_NuRandomService_service
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
                            ${ROOT_EGPYTHIA6}
                            ${GENIE_LIBS_DIR_WITH_LFLAGS}
)


install_headers()
install_fhicl()
install_source()
