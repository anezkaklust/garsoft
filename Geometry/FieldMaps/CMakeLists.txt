# install macro files in lib
file(GLOB txt_file *.txt )
list(FILTER txt_file EXCLUDE REGEX ".*Lists.txt$") 
install_fw( LIST ${txt_file} )