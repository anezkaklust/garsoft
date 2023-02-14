include_directories ( $ENV{NUSIMDATA_INC} )

cet_make_exec(makeDST.exe
  SOURCE makeDST.cxx
  LIBRARIES
                    ROOT::Core
                    ROOT::Hist
                    ROOT::Tree              
    cetlib cetlib_except
)

install_source()
install_headers()
