
include(CetRootCint)

set( PACKAGE EventDisplay )
FILE( GLOB src_files *.cxx )
LIST( REMOVE_ITEM src_files "Display3DPad.cxx")
LIST( REMOVE_ITEM src_files "Display3DView.cxx")

cet_rootcint( ${PACKAGE} )

art_make_library(SOURCE    ${src_files} ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE}Cint.cc
                 LIBRARIES EventDisplay_EVD_ColorDrawingOptions_service
		 	                     EventDisplay_EVD_RawDrawingOptions_service
		                             EventDisplay_EVD_RecoDrawingOptions_service
			                     EventDisplay_EVD_SimulationDrawingOptions_service
 			                     SimulationDataProducts
 			                     RawDataProducts
	                                     DetectorInfo_DetectorPropertiesServiceStandardGAr_service
                                             Geometry_GeometryGAr_service
			                     ReconstructionDataProducts
                                             nusimdata_SimulationBase
                                             nuevdb_EventDisplayBase
			                     nuevdb_EventDisplayBase_Colors_service
                                             ${ART_FRAMEWORK_CORE}
			                     ${ART_FRAMEWORK_PRINCIPAL}
                                             ${ART_FRAMEWORK_IO}
			                     art_Persistency_Common
                                             canvas
		 	                     art_Utilities
 			                     ${MF_MESSAGELOGGER}
                                             ${FHICLCPP}
			                     ${CETLIB}
                           ROOT::Eve
                           ${ROOT_EVE_LIB_LIST}
                           ROOT::EG
                           ${ROOT_RGL}
                           ROOT::X3d
                )

simple_plugin (EVD module
               EventDisplay_EVD
			         nuevdb_EventDisplayBase
               ${ART_FRAMEWORK_CORE}
			         ${ART_FRAMEWORK_PRINCIPAL}
			         ${ART_FRAMEWORK_SERVICES_REGISTRY}
               ${ART_FRAMEWORK_SERVICES_BASIC}
			         ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			         art_Utilities
               canvas
         			 ${MF_MESSAGELOGGER}
	             ${FHICLCPP}
			         ${CETLIB}
     )

simple_plugin (AnalysisDrawingOptions service
			         Geometry_GeometryGAr_service
			         Geometry
			         ReconstructionDataProducts
               nuevdb_EventDisplayBase
               ${ART_FRAMEWORK_CORE}
			         ${ART_FRAMEWORK_PRINCIPAL}
			         ${ART_FRAMEWORK_SERVICES_REGISTRY}
			         ${ART_FRAMEWORK_SERVICES_BASIC}
			         ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
			         art_Persistency_Common
			         art_Persistency_Provenance
			         art_Utilities
               canvas
			         ${MF_MESSAGELOGGER}
	             ${FHICLCPP}
			         ${CETLIB}
                    ROOT::Core
                    ROOT::Hist
                    ROOT::Tree
     )

simple_plugin (ColorDrawingOptions service
			         nuevdb_EventDisplayBase
               ${ART_FRAMEWORK_CORE}
			         ${ART_FRAMEWORK_PRINCIPAL}
			         ${ART_FRAMEWORK_SERVICES_REGISTRY}
			         ${ART_FRAMEWORK_SERVICES_BASIC}
			         ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
			         art_Persistency_Common
			         art_Persistency_Provenance
			         art_Utilities
               canvas
		         	 ${MF_MESSAGELOGGER}
	             ${FHICLCPP}
			         ${CETLIB}
                    ROOT::Core
                    ROOT::Hist
                    ROOT::Tree              

     )

simple_plugin (EvdLayoutOptions service
			         nuevdb_EventDisplayBase
               ${ART_FRAMEWORK_CORE}
			         ${ART_FRAMEWORK_PRINCIPAL}
			         ${ART_FRAMEWORK_SERVICES_REGISTRY}
			         ${ART_FRAMEWORK_SERVICES_BASIC}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
			         art_Persistency_Common
			         art_Persistency_Provenance
			         art_Utilities
               canvas
        			 ${MF_MESSAGELOGGER}
	             ${FHICLCPP}
			         ${CETLIB}
                    ROOT::Core
                    ROOT::Hist
                    ROOT::Tree              
     )

simple_plugin (RawDrawingOptions service
			         nuevdb_EventDisplayBase
               ${ART_FRAMEWORK_CORE}
			         ${ART_FRAMEWORK_PRINCIPAL}
			         ${ART_FRAMEWORK_SERVICES_REGISTRY}
			         ${ART_FRAMEWORK_SERVICES_BASIC}
			         ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
			         art_Persistency_Common
			         art_Persistency_Provenance
			         art_Utilities
               canvas
        			 ${MF_MESSAGELOGGER}
	             ${FHICLCPP}
			         ${CETLIB}
                    ROOT::Core
                    ROOT::Hist
                    ROOT::Tree              
     )

simple_plugin (RecoDrawingOptions service
			         nuevdb_EventDisplayBase
               ${ART_FRAMEWORK_CORE}
         			 ${ART_FRAMEWORK_PRINCIPAL}
         			 ${ART_FRAMEWORK_SERVICES_REGISTRY}
         			 ${ART_FRAMEWORK_SERVICES_BASIC}
         			 ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
         			 art_Persistency_Common
	         		 art_Persistency_Provenance
	         		 art_Utilities
               canvas
	         		 ${MF_MESSAGELOGGER}
               ${FHICLCPP}
		         	 ${CETLIB}
                    ROOT::Core
                    ROOT::Hist
                    ROOT::Tree              
     )

simple_plugin (SimulationDrawingOptions service
               nuevdb_EventDisplayBase
               ${ART_FRAMEWORK_CORE}
               ${ART_FRAMEWORK_PRINCIPAL}
         			 ${ART_FRAMEWORK_SERVICES_REGISTRY}
         			 ${ART_FRAMEWORK_SERVICES_BASIC}
	         		 ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
	         		 art_Persistency_Common
	         		 art_Persistency_Provenance
		         	 art_Utilities
               canvas
	         		 ${MF_MESSAGELOGGER}
               ${FHICLCPP}
		         	 ${CETLIB}
                    ROOT::Core
                    ROOT::Hist
                    ROOT::Tree              
     )



install_headers()
install_source()
