
include_directories( $ENV{GEANT4_FQ_DIR}/include )
include_directories( $ENV{NUTOOLS_DIR}/include   )
include_directories( $ENV{CLHEP_INCLUDE_DIR}     )
include_directories( $ENV{NUG4_INC} )
link_directories( $ENV{NUG4_FQ_DIR}/lib )

set(G4LIBS 
             Geant4::G4digits_hits
  Geant4::G4event
  Geant4::G4geometry
  Geant4::G4global
  Geant4::G4materials
  Geant4::G4particles
  Geant4::G4physicslists
  Geant4::G4processes
  Geant4::G4run
  Geant4::G4track
  Geant4::G4vis_management
  Geant4::G4visHepRep
  Geant4::G4VRML
  Geant4::G4GMocren
  Geant4::G4RayTracer
  Geant4::G4Tree
    )

art_make(LIB_LIBRARIES      nug4_ParticleNavigation
                            nug4_G4Base
                            nusimdata_SimulationBase
                            Geometry
                            GArChannelMapAlg
                            SimulationDataProducts
                            ${ART_FRAMEWORK_PRINCIPAL}
                            ${ART_FRAMEWORK_SERVICES_REGISTRY}
                            ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                            ${ART_ROOT_IO_TFILE_SUPPORT}
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
                            ${G4LIBS}
                            ${Boost_SYSTEM_LIBRARY}
           MODULE_LIBRARIES Geometry
                            GArG4
                            SimulationDataProducts
                            MCCheater_BackTracker_service
                            nusimdata_SimulationBase
                            nug4_G4Base
                            nug4_ParticleNavigation
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
                            ${G4LIBS}
)


install_headers()
install_fhicl()
install_source()

# install macro files in lib
file(GLOB g4_macros *.mac )
install_fw( LIST ${g4_macros} )
