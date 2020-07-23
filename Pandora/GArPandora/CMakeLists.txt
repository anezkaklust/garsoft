include_directories( $ENV{GEANT4_FQ_DIR}/include )
include_directories( $ENV{NUTOOLS_DIR}/include   )
include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{LCPANDORACONTENT_INC} )
include_directories( $ENV{CLHEP_INCLUDE_DIR}     )
include_directories( ../GArContent/include  )

art_make(BASENAME_ONLY
          LIBRARY_NAME GArPandora
          LIB_LIBRARIES
          GArG4
          nusimdata_SimulationBase
          Geometry
          SimulationDataProducts
          RawDataProducts
          ReconstructionDataProducts
          nug4_MagneticField_MagneticField_service
          GArRecoAlg
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
          ${CETLIB_EXCEPT}
          MCCheater_BackTracker_service
          ${PANDORASDK}
          ${PANDORAMONITORING}
          ${LCPANDORACONTENT}
          ${CLHEP}
          GArContent

          MODULE_LIBRARIES
          Geometry
          GArG4
          SimulationDataProducts
          RawDataProducts
          ReconstructionDataProducts
          nusimdata_SimulationBase
          nurandom_RandomUtils_NuRandomService_service
          nug4_MagneticField_MagneticField_service
          GArRecoAlg
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
          ${CETLIB_EXCEPT}
          GArPandora
          MCCheater_BackTracker_service
)

install_headers()
install_fhicl()
install_source()

# install xml files
set( ${product}_xml_dir "${product}/${version}/xml" )
file(GLOB xml_file xml/*.xml )
install( FILES ${xml_file} DESTINATION ${${product}_xml_dir} )
