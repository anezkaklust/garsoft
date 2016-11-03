
include_directories( $ENV{GEANT4_FQ_DIR}/include )
include_directories( $ENV{NUTOOLS_DIR}/include   )
include_directories( $ENV{CLHEP_INCLUDE_DIR}     )

art_make(LIB_LIBRARIES      ${MF_MESSAGELOGGER}
                            ${MF_UTILITIES}
                            ${FHICLCPP}
                            ${CETLIB}
                            ${CLHEP}
                            ${ROOT_BASIC_LIB_LIST}
                            ${ROOT_GEOM}
                            ${ROOT_XMLIO}
                            ${ROOT_GDML}
                            ${ROOT_EG}
                            ${G4FR}
                            ${G4GMOCREN}
                            ${G4OPENGL}
                            ${G4RAYTRACER}
                            ${G4TREE}
                            ${G4VRML}
                            ${G4ANALYSIS}
                            ${G4DIGITS_HITS}
                            ${G4ERROR_PROPAGATION}
                            ${G4EVENT}
                            ${G4GEOMETRY}
                            ${G4GL2PS}
                            ${G4GLOBAL}
                            ${G4GRAPHICS_REPS}
                            ${G4INTERCOMS}
                            ${G4INTERFACES}
                            ${G4MATERIALS}
                            ${G4MODELING}
                            ${G4PARMODELS}
                            ${G4PARTICLES}
                            ${G4PERSISTENCY}
                            ${G4PHYSICSLISTS}
                            ${G4PROCESSES}
                            ${G4READOUT}
                            ${G4RUN}
                            ${G4ERROR_PROPAGATION}
                            ${G4EVENT}
                            ${G4GEOMETRY}
                            ${G4GL2PS}
                            ${G4GLOBAL}
                            ${G4GRAPHICS_REPS}
                            ${G4INTERCOMS}
                            ${G4INTERFACES}
                            ${G4MATERIALS}
                            ${G4MODELING}
                            ${G4PARMODELS}
                            ${G4PARTICLES}
                            ${G4PERSISTENCY}
                            ${G4PHYSICSLISTS}
                            ${G4PROCESSES}
                            ${G4READOUT}
                            ${G4RUN}
                            ${G4TRACK}
                            ${G4TRACKING}
                            ${G4VISHEPREP}
                            ${G4VISXXX}
                            ${G4VIS_MANAGEMENT}
                            ${G4ZLIB}
           MODULE_LIBRARIES Geometry
                            ${ART_FRAMEWORK_CORE}
                            ${ART_FRAMEWORK_IO_SOURCES}
                            ${ART_FRAMEWORK_PRINCIPAL}
                            art_Persistency_Common canvas_Persistency_Common
                            art_Persistency_Provenance canvas_Persistency_Provenance
                            art_Utilities canvas_Utilities
                            ${ART_FRAMEWORK_SERVICES_REGISTRY}
                            ${MF_MESSAGELOGGER}
                            ${MF_UTILITIES}
                            ${FHICLCPP}
                            ${CETLIB}
)


install_headers()
install_fhicl()
install_source()
