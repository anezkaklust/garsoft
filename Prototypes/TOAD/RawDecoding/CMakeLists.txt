simple_plugin(TOADInput "source"
  garsoft::dunedaqhdf5utils2
  art::Framework_Services_Registry
)

add_subdirectory(dunedaqhdf5utils2)
install_headers()
install_fhicl()
install_source()
