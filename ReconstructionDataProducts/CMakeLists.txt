

cet_make( LIBRARIES canvas_Persistency_Common
                    ${MF_MESSAGELOGGER}
          NO_DICTIONARY )

# Simulation base is included in making the dictionary for associations
# between the MC truth and reconstructed quantities when cheating

art_dictionary( DICTIONARY_LIBRARIES ReconstructionDataProducts
                                     RawDataProducts
                                     nusimdata_SimulationBase)

install_headers()
install_source()
