#setup_allpix_cvmfs.sh
# Created on: 15 Feb 2017
#     Author: etahirov
#
# Setup allpix with dependencies from cvmfs
# and locally compiled GEANT4 with some visualisation driver 
# (OpenInventor, OpenGL or else).
#

setupATLAS

lsetup "cmake 3.3.2"

lsetup "root 6.04.02-x86_64-slc6-gcc48-opt"

# Setup install directory for allpix executable
export G4WORKDIR=/home/etahirov/ITk/Allpix/allpix-install/bin
export PATH=$PATH:$G4WORKDIR

# Setup local GEANT4
source /home/etahirov/opt/share/Geant4-10.2.2/geant4make/geant4make.sh

# Working with EUTELESCOPE
export EUTELESCOPE=1

echo 
echo "To compile & install:"
echo "cd /home/etahirov/ITk/Allpix/allpix-build"
echo cmake "-DCMAKE_INSTALL_PREFIX=/home/etahirov/ITk/Allpix/allpix-install -D CMAKE_ECLIPSE_VERSION=\"4.4\" -G\"Eclipse CDT4 - Unix Makefiles\" -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_FLAGS=\"-g2\" -DCMAKE_CXX_FLAGS=\"-g2\" ../allpix
echo
echo "For debugging symbols add manually to CMakeFiles/allpix.dir/flags.make -g -gdwarf-2"
echo "Or: 
   -DCMAKE_BUILD_TYPE:STRING=Debug -DCMAKE_C_FLAGS=\"-g2\" -DCMAKE_CXX_FLAGS=\"-g2\""
echo "make -jN"
echo "make install"

echo "For EUTelescope conversion"
echo "... unfinished"
