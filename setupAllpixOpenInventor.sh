

#gcc4.8
source /afs/cern.ch/sw/lcg/external/gcc/4.8/x86_64-slc6/setup.sh
#source /afs/cern.ch/sw/lcg/contrib/gcc/4.9/x86_64-slc6/setup.sh


#python2.7.4
export PATH="/afs/cern.ch/sw/lcg/external/Python/2.7.4/x86_64-slc6-gcc48-opt/bin:$PATH"
export LD_LIBRARY_PATH="/afs/cern.ch/sw/lcg/external/Python/2.7.4/x86_64-slc6-gcc48-opt/lib:$LD_LIBRARY_PATH" 


#numpy/scipy/sympy

export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/eng/clic/TBData/software/ROOT6_gcc48_python2.7/pytools/numpy/lib/python2.7/site-packages
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/eng/clic/TBData/software/ROOT6_gcc48_python2.7/pytools/scipy/lib/python2.7/site-packages
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/eng/clic/TBData/software/ROOT6_gcc48_python2.7/pytools/sympy/lib/python2.7/site-packages
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/eng/clic/TBData/software/ROOT6_gcc48_python2.7/pytools/cython/lib/python2.7/site-packages
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/eng/clic/TBData/software/ROOT6_gcc48_python2.7/pytools/fastcluster/lib/python2.7/site-packages


#root6
# Use CliC ROOT because of rootcling. If too slow, move ont to gcc4.9.3 and cvmfs ROOT
source /afs/cern.ch/eng/clic/TBData/software/ROOT6_gcc48_python2.7/geant4/root-6.04.00/bin/thisroot.sh
#source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/6.04.14-x86_64-slc6-gcc49-opt/bin/thisroot.sh
#source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/6.04.10-x86_64-slc6-gcc48-opt/bin/thisroot.sh
#source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.18/x86_64-centos7-gcc49-opt

export G4WORKDIR=/home/etahirov/ITk/Allpix/allpix-install/bin
export PATH=$PATH:$G4WORKDIR

source /home/etahirov/opt/share/Geant4-10.2.2/geant4make/geant4make.sh

# Work with EUTELESCOPE
export EUTELESCOPE=1



echo "To compile & install:"
echo "cd /home/etahirov/ITk/Allpix/allpix-build"
echo "cmake -DCMAKE_INSTALL_PREFIX=/home/etahirov/ITk/Allpix/allpix-install -D CMAKE_ECLIPSE_VERSION="4.4" -G\"Eclipse CDT4 - Unix Makefiles\" -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_FLAGS=\"-g2\" -DCMAKE_CXX_FLAGS=\"-g2\" ../allpix"
#echo "/home/etahirov/bin/cmake -DCMAKE_INSTALL_PREFIX=/home/etahirov/allpix/allpix-install -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug ../allpix"
echo "For debugging symbols add manually to CMakeFiles/allpix.dir/flags.make -g -gdwarf-2"
echo "Or: 
   -DCMAKE_BUILD_TYPE:STRING=Debug -DCMAKE_C_FLAGS=\"-g2\" -DCMAKE_CXX_FLAGS=\"-g2\""
echo "make -jN"
echo "make install"

echo "To setup pyLCIO"
echo "source share/setup_pyLCIO.py"
#echo "source /home/etahirov/EUTelescope/v01-17-05/Eutelescope/trunk/build_env.sh" 
#echo "export PYTHONPATH=${LCIO}/src/python:${LCIO}/examples/python:${PYTHONPATH}"

