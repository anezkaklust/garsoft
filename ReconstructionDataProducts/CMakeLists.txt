
cet_make( LIBRARIES
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
                       ${CETLIB_EXCEPT}
                    ROOT::Core
                    ROOT::Hist
                    ROOT::Tree              
                       ${CLHEP}
                       ${Boost_SYSTEM_LIBRARY}
                       nug4_MagneticFieldServices_MagneticFieldServiceStandard_service
                       NO_DICTIONARY )

# Simulation base is included in making the dictionary for associations
# between the MC truth and reconstructed quantities when cheating

art_dictionary( DICTIONARY_LIBRARIES ReconstructionDataProducts
                                     RawDataProducts
                                     nug4_MagneticFieldServices_MagneticFieldServiceStandard_service
                                     nusimdata_SimulationBase)

install_headers()
install_source()
