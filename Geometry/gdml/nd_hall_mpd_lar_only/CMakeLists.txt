# install schema files
file(GLOB gdml_file *.gdml)
install( FILES ${gdml_file} DESTINATION ${${product}_gdml_dir}/nd_hall_mpd_lar_only )
