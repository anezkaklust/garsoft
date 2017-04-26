
art_make( NO_PLUGINS
          LIB_LIBRARIES Geometry
	                   		CoreUtils
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${CETLIB}
                        ${CETLIB_EXCEPT}
                  			${FHICLCPP}
                        ${ROOT_BASIC_LIB_LIST}
                  			${Boost_SYSTEM_LIBRARY}
        )


simple_plugin(DetectorClocksServiceStandard "service"
              DetectorInfo
              RawDataProducts
              ${ART_FRAMEWORK_CORE}
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_PERSISTENCY_ROOTDB}
              art_Utilities
              canvas
              art_Persistency_Common
              art_Persistency_Provenance 
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${FHICLCPP}
              ${CETLIB_EXCEPT}
              ${SQLITE3}
              ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(DetectorPropertiesServiceStandard "service"
              DetectorInfo
              Geometry_Geometry_service
              Geometry
              ${ART_PERSISTENCY_ROOTDB}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              art_Utilities canvas
              ${SQLITE3}
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${CETLIB_EXCEPT}
              ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(GArPropertiesServiceStandard "service"
              DetectorInfo
              art_Utilities canvas
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${CETLIB_EXCEPT}
             )


install_headers()
install_fhicl()
install_source()

#add_subdirectory(test)

