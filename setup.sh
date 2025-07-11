#!/bin/sh
# # to avoid linking shared
# export CXXFLAGS="$CXXFLAGS -fPIC"


external_libs_full_path=$(readlink -f ./external_libs)
rm -rf ${external_libs_full_path}
# rm -rf ${external_libs_full_path}/lib64
mkdir ${external_libs_full_path}
cd ${external_libs_full_path}

wget "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz"
tar -xzvf "./eigen-3.4.0.tar.gz"
cd eigen-3.4.0
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=${external_libs_full_path}
cmake --build .
cmake --install .


cd ${external_libs_full_path}
git clone "https://github.com/imneme/pcg-cpp.git"
# rm -rf ${external_libs_full_path}/libleidenalg
# rm -rf ${external_libs_full_path}/lib
# rm -rf ${external_libs_full_path}/lib64
# mkdir ${external_libs_full_path}
# cd ${external_libs_full_path}
# install igraph to external_libs
# echo "installing igraph to ${external_libs_full_path}"
# git clone https://github.com/igraph/igraph.git
# cd igraph
# mkdir build
# cd build
# cmake .. -DCMAKE_INSTALL_PREFIX=${external_libs_full_path} -DIGRAPH_ENABLE_TLS=ON -DIGRAPH_ENABLE_LTO=ON
# cmake --build .
# cmake --install .

# install libleidenalg to external_libs
# cd ${external_libs_full_path}
# echo "installing libleidenalg to ${external_libs_full_path}"
# git clone https://github.com/vtraag/libleidenalg.git
# cd libleidenalg
# mkdir build
# cd build
# cmake .. -DCMAKE_INSTALL_PREFIX=${external_libs_full_path} -DCMAKE_PREFIX_PATH=${external_libs_full_path} -DBUILD_SHARED_LIBS=FALSE
# cmake --build .
# cmake --build . --target install
