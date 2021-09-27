include_directories( include/ )
include_directories( $ENV{PANDORA_INC} )

cet_make_library(LIBRARY_NAME GArContent
  SOURCE
  src/GArObjects/Helix.cc
  LIBRARIES
  PANDORASDK
  PANDORAMONITORING
  LCPANDORACONTENT
  ${CLHEP}
)

install_source( SUBDIRS src/GArObjects/ )
install_headers( SUBDIRS include/GArObjects/ )
