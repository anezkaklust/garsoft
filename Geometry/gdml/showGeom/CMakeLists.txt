
cet_make_exec(showGeometry
  SOURCE
    showGeometry.cpp
  LIBRARIES
    ROOT_BASIC_LIB_LIST
    ${ROOT_EG}
    ${ROOT_EVE}
    ${ROOT_EVE_LIB_LIST}
)

install_source()
