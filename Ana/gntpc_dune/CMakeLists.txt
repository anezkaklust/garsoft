include_directories ( $ENV{GENIE_INC}/GENIE )
include_directories ( $ENV{GENIE_INC} )
include_directories ( $ENV{MESSAGEFACILITY_INC} )

execute_process(COMMAND genie-config --libs OUTPUT_VARIABLE GENIE_LIBS_DIR_WITH_LFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)

include_directories("$ENV{LIBXML2_INC}")

cet_make_exec(gntpc_dune
  SOURCE
    main/gNtpConvDUNE.cxx
  LIBRARIES
                    ROOT::Core
                    ROOT::Hist
                    ROOT::Tree              
    ROOT::Eve
    ROOT::Core
    ${ROOT_EVE_LIB_LIST}
    ${ROOT_EGPYTHIA6}
    ROOT::EG
    ROOT::MathCore
    ${ROOT_GENVECTOR}
    ${GENIE_LIBS_DIR_WITH_LFLAGS}
    ${LIBXML2_LIBRARIES}
    messagefacility::MF_MessageLogger
)

install_source()
install_headers()
