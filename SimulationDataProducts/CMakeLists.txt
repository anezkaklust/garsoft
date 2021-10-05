
cet_make( LIBRARIES nusimdata_SimulationBase
                    canvas
                    MF_MessageLogger
                    
                    cetlib
                    ROOT_BASIC_LIB_LIST
          NO_DICTIONARY )

art_dictionary( DICTIONARY_LIBRARIES SimulationDataProducts )

install_headers()
install_fhicl()
install_source()
