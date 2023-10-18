
add_subdirectory(scripts)
add_subdirectory(Dumpers)
add_subdirectory(Metadata)

art_make(
  LIB_LIBRARIES
    ${ROOT_EG}
    ${ROOT_CORE}
    nusimdata_SimulationBase
    ReconstructionDataProducts
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_PRINCIPAL}
    art_Persistency_Common
    art_Persistency_Provenance
    art_Utilities
    canvas
    ${MF_MESSAGELOGGER}
    ${ROOT_EG}
                    ROOT::Core
                    ROOT::Hist
                    ROOT::Tree              
    ${CETLIB_EXCEPT}
    ${FHICLCPP}
    ${Boost_SYSTEM_LIBRARY}
    ${CLHEP}

  MODULE_LIBRARIES
    Utilities
    nusimdata_SimulationBase
    ReconstructionDataProducts
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    art_Persistency_Common
    art_Persistency_Provenance
    art_Utilities
    canvas
    ${MF_MESSAGELOGGER}
    ${ROOT_EG}
                    ROOT::Core
                    ROOT::Hist
                    ROOT::Tree              
    ${FHICLCPP}
    ${CETLIB_EXCEPT}
  )

add_subdirectory(Converter)
add_subdirectory(RootGraphics)

install_headers()
install_fhicl()
install_source()
