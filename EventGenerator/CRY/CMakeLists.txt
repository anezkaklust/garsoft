
include_directories ( $ENV{CRYHOME}/src )

art_make(MODULE_LIBRARIES Geometry
                          Geometry_GeometryGAr_service
                          SummaryDataProducts
                          nutools_EventGeneratorBase_CRY
                          nusimdata_SimulationBase
                          nurandom_RandomUtils_NuRandomService_service
                          ${ART_FRAMEWORK_CORE}
                          ${ART_FRAMEWORK_PRINCIPAL}
                          ${ART_FRAMEWORK_SERVICES_REGISTRY}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                          ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                          ${ART_ROOT_IO_TFILE_SUPPORT}
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
        )

install_headers()
install_fhicl()
install_source()

