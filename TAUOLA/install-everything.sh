if test -e documentation; then
  echo ""
  echo "Install-everything.sh has to be moved outside TAUOLA directory first"
  echo "consider mv install-everything.sh ../.  before executing"
  echo "also check if packages to be installed exist at the predefined links"
  echo "After its execution TAUOLA directory will be moved to (freshly created)"
  echo " external/TAUOLA/TAUOLA"
  echo ""
  sleep 15
  exit
fi


if test -e external; then
  echo ""
  echo "Directory 'external' exists. Remove it then re-run the script"
  echo ""
  sleep 15
  exit
fi
# ####
#  0 # 
# ####
#      This script and its command are expected to be in higher level
#      directory. The directory external and  subdirectories will be created
#      as external/TAUOLA, external/LHAPDF, etc.

# ####
#  1 #  
# ####
#      Get all packages TAUOLA and TauSpinner or their examples use
#
# POSSIBLE VARIANTS:
# If you have some of them installed already, you may skip some of the
# lines below, but still there may be cross dependencies which need to
# match. 
#
mkdir external
cd external/
mkdir  TAUOLA MC-TESTER LHAPDF HEPMC PYTHIA

cd TAUOLA/
#wget https://tauolapp.web.cern.ch/tauolapp/resources/TAUOLA.1.1.6b/TAUOLA.1.1.6b.tar.gz
#tar -xvzf TAUOLA.1.1.6b.tar.gz
#rm TAUOLA.1.1.6b.tar.gz
mv ../../TAUOLA .
cd ../MC-TESTER/
wget https://mc-tester.web.cern.ch/MC-TESTER/MC-TESTER-1.25.0.tar.gz
tar -xvzf MC-TESTER-1.25.0.tar.gz
rm MC-TESTER-1.25.0.tar.gz
cd ../HEPMC/
wget http://lcgapp.cern.ch/project/simu/HepMC/download/HepMC-2.06.09.tar.gz
tar -xvzf HepMC-2.06.09.tar.gz
rm HepMC-2.06.09.tar.gz
cd ../LHAPDF/
wget http://www.hepforge.org/archive/lhapdf/lhapdf-5.9.1.tar.gz
tar -xvzf lhapdf-5.9.1.tar.gz
rm lhapdf-5.9.1.tar.gz
cd ../PYTHIA/
wget http://home.thep.lu.se/~torbjorn/pythia8/pythia8201.tgz
tar -xvzf pythia8201.tgz
rm pythia8201.tgz

##############
#  2 Builds: #
##############

#
# Build pythia
#
cd pythia8201/
./configure
make -j4

#
# Build HepMC
#
cd ../../HEPMC/HepMC-2.06.09/
./configure --with-momentum=GEV --with-length=MM --prefix=`pwd`/../HepMC
make install -j4

#
# Build MC-TESTER
#
cd ../../MC-TESTER/MC-TESTER/
./configure --with-HepMC=`pwd`/../../HEPMC/HepMC
make -j4

#
# Build LHAPDF
#
cd ../../LHAPDF/lhapdf-5.9.1/
./configure --prefix=`pwd`/../lhapdf
make install -j4
cd ../lhapdf/share/lhapdf/
../../bin/lhapdf-getdata MSTW2008nnlo90cl.LHgrid
../../bin/lhapdf-getdata cteq6ll.LHpdf

#
# Build/make TAUOLA TauSpinner, TauSpinner examples build must be done after Tauola/TauSpinner build.
#
cd ../../../../TAUOLA/TAUOLA/
./configure --with-tau-spinner --with-lhapdf=`pwd`/../../LHAPDF/lhapdf --with-mc-tester=`pwd`/../../MC-TESTER/MC-TESTER --with-hepmc=`pwd`/../../HEPMC/HepMC --with-pythia8=`pwd`/../../PYTHIA/pythia8201
# NOTE, you may need PYTHIA8DATA= for older versions of Pythia
make
source configure.paths.sh

cd TauSpinner/examples/
./configure --with-mc-tester=$MCTESTERLOCATION --with-hepmc=$HEPMCLOCATION

make
./tau-reweight-test.exe

######
# 3  # Build and run of user analysis program sequence may look as follows
######
  echo "Finally you may like to add to setup.sh"
  echo "export SPIN2LIB=`pwd`/../../../TAUOLA/TauSpinner/examples/example-VBF/SPIN2/lib"
  echo "or link to your own directory for calculation of alternative tautau2j matrix elements"
# cd ../..
# source setup.sh
# cd ../../..
# cd myAnalysis
# make

