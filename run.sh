#!/bin/bash

cd ../pyscan/build
rm -rf *
cmake ..
make -j8
cd ../../trajectory_paper_code

cp -rf ../pyscan/build/libpyscan.so .
cp -rf ../pyscan/build/pyscan.py .
