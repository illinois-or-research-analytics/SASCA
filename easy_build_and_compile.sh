#!/bin/sh
rm ./abm
rm -r build
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
# cmake -DCMAKE_BUILD_TYPE=Debug ..
cp compile_commands.json ../
make -j 1
ln -s $(readlink -f ./bin/abm) ../abm
