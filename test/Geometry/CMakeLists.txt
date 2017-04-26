cet_enable_asserts()

# test libraries
art_make_library( LIBRARY_NAME GeometryTestLib
                  SOURCE GeometryTestAlg.cxx
                  LIBRARIES Geometry
                            TestUtils
                            art_Utilities canvas
                            ${MF_MESSAGELOGGER}
                            ${MF_UTILITIES}
                            ${FHICLCPP}
                            ${CETLIB_EXCEPT}
                            ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
          )


# geometry test on "standard" geometry

# unit test (use the hard-coded channel mapping for "standard" GArTPCdetector)
cet_test(geometry_test
         SOURCES   geometry_test.cxx
         DATAFILES test_geometry.fcl
         TEST_ARGS ./test_geometry.fcl
         LIBRARIES Geometry
                   GeometryTestLib
                   ${MF_MESSAGELOGGER}
                   ${MF_UTILITIES}
                   ${FHICLCPP}
                   ${CETLIB_EXCEPT}
         USE_BOOST_UNIT
)

# set FHICL_FILE_PATH and FW_SEARCH_PATH for the case where building without mrb
# this is not generic, we take advantage of knowing the subdirectories for the fcl and gdml files
set( mrb_build_dir $ENV{MRB_BUILDDIR} )
if( NOT mrb_build_dir )
  #message(STATUS "MRB_BUILDDIR is undefined, set FHICL_FILE_PATH and FW_SEARCH_PATH")
  set_tests_properties( geometry geometry_iterator 
                        PROPERTIES ENVIRONMENT 
                        "FHICL_FILE_PATH=.:${CMAKE_BINARY_DIR}/fcl;FW_SEARCH_PATH=${CMAKE_BINARY_DIR}/gdml" )
endif()



install_headers()
install_fhicl()
install_source()
