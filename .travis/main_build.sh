#!/bin/bash

set -e

mkdir build
cd build

### Cmake
echo "Using C++ = $CXXCOMPILER"
cmake ..  $BTYPE -DCMAKE_CXX_COMPILER=$CXXCOMPILER -DCMAKE_C_COMPILER=$CCOMPILER 

if [ $DEC = "true"];
    then
        echo "Compile Dec not in //";
        cd examples; 
        make exampleDiscreteExteriorCalculusChladni;
        make exampleDiscreteExteriorCalculusSolve;
        make exampleDECSurface;
        make examplePropagation;
        cd ../tests;
        make testDiscreteExteriorCalculusExtended;        
fi
    
### DGtal  build
make -j 4

echo "NeedExample $NEEDEXAMPLESANDTESTS"
### DGtal Examples and Examples
if [ $NEEDEXAMPLESANDTESTS = "true" ];
then
   if [ $DEC = "true"];
    then
        echo "Compile Dec not in //";
        cd examples; 
        make exampleDiscreteExteriorCalculusChladni;
        make exampleDiscreteExteriorCalculusSolve;
        make exampleDECSurface;
        make examplePropagation;
        cd ../tests;
        make testDiscreteExteriorCalculusExtended;        
   fi
    
    cd examples ; make  -j 3
    cd ../tests ;  make -j 3

    if [ -f io/writers/testMagickWriter ]; then
      io/writers/testMagickWriter -s
    fi
    if [ -f io/readers/testMagickReader ]; then
      io/readers/testMagickReader
    fi
    if [ $DEC = "true"];
    then
        echo "Compile Dec not in //";
        make testDiscreteExteriorCalculusExtended;
        make exampleDiscreteExteriorCalculusChladni;
        make exampleDiscreteExteriorCalculusSolve;
        make exampleDECSurface;
        make examplePropagation;
    fi
    ctest -j 3--output-on-failure


fi

### DGtal doc
if [ $BUILD_DOC = "true" ];
then
    echo "Building the doc"
    make doc > buildDoc.log
fi
