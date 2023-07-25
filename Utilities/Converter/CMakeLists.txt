include_directories( $ENV{CLHEP_INCLUDE_DIR} )
include_directories( $ENV{GENIE_INC}/GENIE )
include_directories( $ENV{GENIE_INC} )
include_directories( $ENV{EDEPSIM_INC}/EDepSim )
link_directories( $ENV{NUG4_FQ_DIR}/lib )

execute_process(COMMAND genie-config --libs OUTPUT_VARIABLE GENIE_LIBS_DIR_WITH_LFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)

art_make(LIB_LIBRARIES      SimulationDataProducts
                            SummaryDataProducts
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
                            ${CETLIB}
                            ${CETLIB_EXCEPT}
                            ${ROOT_EG}
                            ${ROOT_TREEPLAYER}
                            ${ROOT_FFTW}
                            ${ROOT_REFLEX}
                            ${ROOT_EGPYTHIA6}
                            ${ROOT_GUI}
                    ROOT::Core
                    ROOT::Hist
                    ROOT::Tree              
                            ${ROOT_GEOM}
                            ${CLHEP}
                            ${Boost_SYSTEM_LIBRARY}
			    EDepSim::edepsim_io
       MODULE_LIBRARIES     SimulationDataProducts
                            SummaryDataProducts
                            Geometry
                            MCCheater_BackTracker_service
                            nugen_EventGeneratorBase_GENIE
                            nurandom_RandomUtils_NuRandomService_service
                            ${ART_FRAMEWORK_CORE}
                            ${ART_FRAMEWORK_IO_SOURCES}
                            ${ART_FRAMEWORK_PRINCIPAL}
                            art_Persistency_Common
                            art_Persistency_Provenance
                            art_Utilities
                            canvas
                            ${ART_FRAMEWORK_SERVICES_REGISTRY}
                            ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                            ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                            ${MF_MESSAGELOGGER}
                            ${FHICLCPP}
                            ${CETLIB}
                            ${CETLIB_EXCEPT}
                    ROOT::Core
                    ROOT::Hist
                    ROOT::Tree              
                            ${ROOT_EGPYTHIA6}
                            ${GENIE_LIBS_DIR_WITH_LFLAGS}
			    EDepSim::edepsim_io
)


install_headers()
install_fhicl()
install_source()
