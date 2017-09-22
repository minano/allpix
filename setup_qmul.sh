#source /afs/cern.ch/sw/lcg/external/gcc/4.8/x86_64-slc6/setup.sh
#source /home/mbenoit/workspace/GEANT4/root-6.04.00/bin/thisroot.sh 
#source /opt/root_v6.06.00/bin/thisroot.sh
source /opt/root_v6.08.02/bin/thisroot.sh

#source /home/mbenoit/workspace/GEANT4/geant4.10.1-install/bin/geant4.sh
source /users/tahirovic/opt/geant4-install/bin/geant4.sh

# this is already defined in my .bashrc
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/cern.ch/eng/clic/software/Pixel_TestBeam_Software/xerces/lib

export G4WORKDIR=$HOME/ITk/Allpix/allpix-install/bin
export PATH=$G4WORKDIR:$PATH
#export LD_LIBRARY_PATH=/home/mbenoit/miniconda/lib:$LD_LIBRARY_PATH

#echo "$HOME/opt/cmake-3.3.2-Linux-x86_64/bin/cmake"
echo "$HOME/opt/cmake-3.3.2-Linux-x86_64/bin/cmake -DCMAKE_INSTALL_PREFIX=/users/tahirovic/ITk/Allpix/allpix-install -D CMAKE_ECLIPSE_VERSION=\"4.4\" -G\"Eclipse CDT4 - Unix Makefiles\" -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_FLAGS=\"-g2\" -DCMAKE_CXX_FLAGS=\"-g2\" ../allpix"

# Working with EUTELESCOPE
export EUTELESCOPE=1
export _EUTELESCOPE=1
