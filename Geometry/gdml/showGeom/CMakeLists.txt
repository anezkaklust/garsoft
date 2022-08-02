
cet_make_exec(showGeometry
  SOURCE
    showGeometry.cpp
  LIBRARIES
                    ROOT::Core
                    ROOT::Hist
                    ROOT::Tree              
		    ROOT::Geom
		    ROOT::EG
		    ROOT::Eve
)

install_source()
