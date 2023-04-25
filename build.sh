#!/bin/bash

cmake -B _build/ -S . -DCMAKE_INSTALL_PREFIX=$PWD #-DUSE_OPENMP=ON #-DUSE_MPI=ON 

cmake --build _build --target install $*