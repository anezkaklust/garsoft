# install schema files
file(GLOB gdml_file *.gdml)
install( FILES ${gdml_file} DESTINATION ${${product}_gdml_dir}/MPD_Standalone )
