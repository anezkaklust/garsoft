
include_directories( $ENV{NUSIMDATA_INC} )

cet_make( LIBRARIES canvas
		    MF_MessageLogger
		    MF_Utilities
		    cetlib cetlib_except
		    ${ROOT_BASIC_LIB_LIST}
                    ${ROOT_EG}
                    NO_DICTIONARY )

art_dictionary( DICTIONARY_LIBRARIES RawDataProducts )

install_headers()
install_source()
install_fhicl()
