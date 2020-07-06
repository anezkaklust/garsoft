source /grid/fermiapp/products/dune/setup_dune.sh
setup cmake v3_12_1
setup gcc v6_4_0
setup jobsub_client
setup pycurl

setup dk2nu        v01_05_01b -q e15:prof
setup genie        v2_12_10c  -q e15:prof
setup genie_xsec   v2_12_10   -q DefaultPlusValenciaMEC
setup genie_phyopt v2_12_10   -q dkcharmtau
setup geant4 v4_10_3_p01b -q e15:prof

# we need something called TBB for ROOT
export TBB_DIR=/cvmfs/larsoft.opensciencegrid.org/products/tbb/v2018_1
export TBB_LIB=/cvmfs/larsoft.opensciencegrid.org/products/tbb/v2018_1/Linux64bit+2.6-2.12-e15-prof/lib
export TBB_INC=/cvmfs/larsoft.opensciencegrid.org/products/tbb/v2018_1/Linux64bit+2.6-2.12-e15-prof/include

G4_cmake_file=`find ${GEANT4_FQ_DIR}/lib64 -name 'Geant4Config.cmake'`
export Geant4_DIR=`dirname $G4_cmake_file`

export PATH=$PATH:$GEANT4_FQ_DIR/bin

export GXMLPATH=$GXMLPATH:/nashome/m/mtanaz/HPgTPC_sim/

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}/edep-sim/lib

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}/nusystematics/build/Linux/lib:${PWD}/nusystematics/artless

export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$GENIE_INC/GENIE
export NUSYST=`pwd`/nusystematics
export LD_LIBRARY_PATH=$NUSYST/build/Linux/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$NUSYST/build/nusystematics/artless:$LD_LIBRARY_PATH

echo "cmake ../ -DUSEART=0 -DLIBXML2_LIB=/cvmfs/larsoft.opensciencegrid.org/products/libxml2/v2_9_5/Linux64bit+2.6-2.12-prof/lib/ -DLIBXML2_INC=/cvmfs/larsoft.opensciencegrid.org/products/libxml2/v2_9_5/Linux64bit+2.6-2.12-prof/include/libxml2 -DPYTHIA6=/cvmfs/larsoft.opensciencegrid.org/products/pythia/v6_4_28i/Linux64bit+2.6-2.12-gcc640-prof/lib"
