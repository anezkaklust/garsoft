
cet_make( LIBRARIES canvas
                    MF_MessageLogger
                    cetlib cetlib_except
                    ROOT::Core
                    ROOT::Hist
                    ROOT::Tree              
                    ${ROOT_EG}
                    NO_DICTIONARY )

art_dictionary( DICTIONARY_LIBRARIES RawDataProducts )

install_headers()
install_source()
install_fhicl()
