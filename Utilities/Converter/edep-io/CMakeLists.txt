# Set the default verbosity of the compilation...
# set(CMAKE_VERBOSE_MAKEFILE TRUE)

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

# Compile the base library with private I/O fields.
add_definitions(-DEDEPSIM_FORCE_PRIVATE_FIELDS)

cet_make( LIBRARIES
		    ${ROOT_BASIC_LIB_LIST}
        ${ROOT_EG}
        NO_DICTIONARY )

art_dictionary(DICTIONARY_LIBRARIES Utilities_Converter_edep-io NO_CHECK_CLASS_VERSION)

install_source()
