#!/bin/bash

source setup_slc6.sh
source myGeant4.sh 
export CXX=`which g++`; export CC=`which gcc`

cd build 
rm -rf *
cmake ../   

#cd build 
cp ../vis.mac ./ 

make
./shower 50