include_directories( $ENV{GEANT4_FQ_DIR}/include )
include_directories( $ENV{NUTOOLS_DIR}/include   )
include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{LCPANDORACONTENT_INC} )

art_make(BASENAME_ONLY
          LIBRARY_NAME GArPandora
          LIB_LIBRARIES
          GArG4
          nusimdata_SimulationBase
          Geometry
          SimulationDataProducts
          RawDataProducts
          ReconstructionDataProducts
          nutools_MagneticField_MagneticField_service
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
          ${CETLIB_EXCEPT}

          ${PANDORASDK}
          ${PANDORAMONITORING}
          ${LCPANDORACONTENT}

          MODULE_LIBRARIES 
          Geometry
          GArG4
          SimulationDataProducts
          RawDataProducts
          ReconstructionDataProducts
          nusimdata_SimulationBase
          nutools_RandomUtils_NuRandomService_service
          nutools_MagneticField_MagneticField_service
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
          ${CETLIB_EXCEPT}
          GArPandora
)

install_headers()
install_fhicl()
install_source()

# install xml files
set( ${product}_xml_dir "${product}/${version}/xml" )
file(GLOB xml_file xml/*.xml )
install( FILES ${xml_file} DESTINATION ${${product}_xml_dir} )