
art_make(LIB_LIBRARIES      ${MF_MESSAGELOGGER}
                            ${MF_UTILITIES}
                            ${FHICLCPP}
                            ${CETLIB}
                            ${CETLIB_EXCEPT}
                            ${CLHEP}
                            ${ROOT_BASIC_LIB_LIST}
                            ${ROOT_GEOM}
                            ${ROOT_XMLIO}
                            ${ROOT_GDML}
                            ${ROOT_EG}
          SERVICE_LIBRARIES Geometry
                            ${ART_FRAMEWORK_CORE}
                            ${ART_FRAMEWORK_IO_SOURCES}
                            ${ART_FRAMEWORK_PRINCIPAL}
                            art_Persistency_Common
                            art_Persistency_Provenance
                            art_Utilities
                            canvas
                            ${ART_FRAMEWORK_SERVICES_REGISTRY}
                            ${MF_MESSAGELOGGER}
                            ${MF_UTILITIES}
                            ${FHICLCPP}
                            ${CETLIB}
                            ${CETLIB_EXCEPT}

          )

install_headers()
install_fhicl()
install_source()

add_subdirectory(gdml)
