
cet_make( LIBRARIES nusimdata_SimulationBase
                    canvas_Persistency_Provenance
                    canvas_Utilities
                    MF_MessageLogger
                    MF_Utilities
                    cetlib
                    ${ROOT_BASIC_LIB_LIST}
          NO_DICTIONARY )

art_dictionary( DICTIONARY_LIBRARIES larsimobj_Simulation )

install_headers()
install_fhicl()
install_source()
