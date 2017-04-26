# ======================================================================
#
# Testing
#
# ======================================================================

cet_test( printBacktrace_test
          LIBRARIES CoreUtils
  )

cet_test(ServiceUtil_test
        LIBRARIES canvas
                  ${MF_MESSAGELOGGER}
                  ${MF_UTILITIES}
                  ${FHICLCPP}
                  ${CETLIB}
                  ${CETLIB_EXCEPT}
        USE_BOOST_UNIT
  )

cet_test(UncopiableAndUnmovableClass_test
        USE_BOOST_UNIT
  )

cet_test( ProviderPack_test USE_BOOST_UNIT )

