include_directories( $ENV{CETLIB_EXCEPT_INC} )
include_directories( $ENV{ART_INC} )

art_make(LIB_LIBRARIES ${CETLIB_EXCEPT} canvas
        )


install_headers()
install_source()

