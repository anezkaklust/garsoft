
include_directories ( $ENV{GENIE_INC}/GENIE )

# genie
execute_process(COMMAND genie-config --libs OUTPUT_VARIABLE GENIE_LIBS_DIR_WITH_LFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)

art_make( MODULE_LIBRARIES SimulationDataProducts
                           nusimdata_SimulationBase
                           SummaryDataProducts
                           Geometry
                           Geometry_Geometry_service
                           nutools_EventGeneratorBase_GENIE
                           nutools_RandomUtils_NuRandomService_service
                           ${ART_FRAMEWORK_CORE}
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
                           ${CLHEP}
    			                 ${GENIE_LIBS_DIR_WITH_LFLAGS}
                           ${ROOT_EGPYTHIA6}    # FIXME!!! - resolving genie run time reference
                           ${ROOT_BASIC_LIB_LIST}
                           ${ROOT_EG}
        )

install_headers()
install_fhicl()
install_source()


# install xml files
set( ${product}_xml_dir "${product}/${version}/xml" )
file(GLOB xml_file *.xml )
install( FILES ${xml_file} DESTINATION ${${product}_xml_dir} )
