# garsoft: Software for ND-GAr and TOAD

## Building GArSoft from Source -- all other dependencies in CVMFS

To set up a new test release for doing development in your own area, execute the following commands (or put them in a shell script) to set up a new development area.
```
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
mkdir <new empty directory>
cd <new emtpy directory>
export MRB_PROJECT=garsoft
export COMPILER=e20
export BUILDTYPE=prof
mrb newDev -v develop -q ${COMPILER}:${BUILDTYPE}
source localProducts*/setup
mkdir work
cd srcs
mrb g garsoft # git clone OG garsoft and include in the CMakeLists.txt
rm -r garsoft # remove OG garsoft src
git clone https://github.com/anezkaklust/garsoft.git # replace by toad garsoft
cd $MRB_BUILDDIR
mrbsetenv
mrb i -j4
mrbslp
```

The COMPILER variable above is a UPS qualifier indicating which version of which compiler to use. The head of develop has been upgraded to art v3_09_03, which uses a compiler qualifier "e20", which corresponds to GCC v9.3.0: https://cdcvs.fnal.gov/redmine/projects/cet-is-public/wiki/AboutQualifiers

The BUILDTYPE variable above is either "prof" or "debug". Both kinds of builds include debug symbols, but "prof" turns on optimization, which can make using a debugger more challenging, but will make the code run faster.

The -j4 argument on the mrb i command refers to the number of concurrent build processes (like g++ compiling a single source file).  A good rule of thum is to match it to the number of CPUs available on the machine you are working.  DUNE gpvms have four cores and a commensurate amount of memory, so use -j4 on those.  DUNE build nodes have 16 cores each and enough memory to support running the compiler on each one.

Each time you log in and want to work with your garsoft test release, execute (or write a script), do the following to set up a build environment:

```
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
export MRB_PROJECT=garsoft
source <directory you created from the creation step above>/localProducts*/setup
cd ${MRB_BUILDDIR}
mrbsetenv
```

Note: mrbsetenv sets up a build environment, and environment variables will point to the build directory. mrbslp will set up installed local products for running. Sometimes this makes a difference -- if a build fails, the build directory can have an incomplete set of things in it. It may also be missing some items (pandora XML files is a common one in far detector sim/reco that only get picked up with mrbslp).

See the DUNE computing [training](https://dune.github.io/computing-basics/) for a tutorial on how to develop software.  Some older examples from dunetpc on how to work with git, check out code, develop on branches, and pushing code: [dunetpc git/mrb tutorial](https://cdcvs.fnal.gov/redmine/projects/dunetpc/wiki/_Tutorial_)

## Run a 2 sample events for TOAD

Job after building or setting up code:

```
# event generation
art -n 2 -c toadtfgen.fcl
# readout simulation
art -c readoutsimjob.fcl text_gen.root 
# event display
art -c evd3D.fcl readoutsim.root
```

Reconstrution job for TOAD consisting of hit clustering can be run as:

```
art -c recojob.fcl readoutsim.root   
```

An analysis root tree with can be made
```
art -c anajob.fcl reco.root
```

There is an example root macro that reads the anatree.root file made by anajob.fcl. It's in garsoft/Ana/ExampleAnalysisScripts/garana.C and garana.h. There's no CMakeLists.txt file in there so the script doesn't get installed anywhere by default, so it just lives in the git repository. Copy it to your working directory along with the header file and run it to make some plots.


## Event Display tips
The event display starts a ROOT TApplication, and reads rootlogon.C. Users may have a rootlogon.C which is incompatible with the event display -- the symptom is a segfault starting up the event display. It works fine with a blank rootlogon.C, which you can put in the directory you're running the event display from so as not to clobber your regular rootlogon.C. 

For graphics visualization, one can set up VNC port.

## Computing tips
Lots of computing how-to information is here: https://wiki.dunescience.org/wiki/DUNE_Computing/Computing_How-To_Documentation

## Copyright and Licensing
Copyright © 2023 FERMI NATIONAL ACCELERATOR LABORATORY for the benefit of the DUNE Collaboration.

This repository, and all software contained within, except where noted within the individual source files, is licensed under
the Apache License, Version 2.0 (the "License"); you may not use this
file except in compliance with the License. You may obtain a copy of
the License at

    http://www.apache.org/licenses/LICENSE-2.0

Copyright is granted to FERMI NATIONAL ACCELERATOR LABORATORY on behalf
of the Deep Underground Neutrino Experiment (DUNE). Unless required by
applicable law or agreed to in writing, software distributed under the
License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for
the specific language governing permissions and limitations under the
License.
