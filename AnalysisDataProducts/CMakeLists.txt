
cet_make( LIBRARIES nusimdata_SimulationBase
                    canvas
                    MF_MessageLogger
                    cetlib
                    ${ROOT_BASIC_LIB_LIST}
                    GaranaProducts
          NO_DICTIONARY )

art_dictionary( DICTIONARY_LIBRARIES AnalysisDataProducts 
                                     ReconstructionDataProducts )

install_headers()
install_fhicl()
install_source()
