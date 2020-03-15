include_directories ( $ENV{GENIE_INC}/GENIE )

execute_process(COMMAND genie-config --libs OUTPUT_VARIABLE GENIE_LIBS_DIR_WITH_LFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)

cet_find_library(log4cpp NAMES log4cpp PATHS $ENV{LOG4CPP_LIB} PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
include_directories("$ENV{LOG4CPP_INC}")

cet_find_library(pythia6 NAMES Pythia6 pythia6 PATHS $ENV{PYTHIA6_LIBRARY} PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
include_directories("$ENV{PYTHIA_INC}")

cet_find_library(libxml2 NAMES libxml2 xml2 PATHS $ENV{LIBXML2_FQ_DIR} PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
include_directories("$ENV{LIBXML2_INC}")

cet_make_exec(gntpc_dune
  SOURCE
    main/gNtpConvDUNE.cxx
  LIBRARIES
    ${ROOT_BASIC_LIB_LIST}
    ${ROOT_EVE}
    ${ROOT_EVE_LIB_LIST}
    ${ROOT_EGPYTHIA6}
    ${ROOT_EG}
    ${ROOT_GENVECTOR}
    ${GENIE_LIBS_DIR_WITH_LFLAGS}
    ${log4cpp}
    ${pythia6}
    ${libxml2}
)

install_source()
install_headers()
