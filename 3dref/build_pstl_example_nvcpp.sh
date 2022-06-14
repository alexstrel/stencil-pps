#!/usr/bin/bash
###
#nvc++ -O0 -g -std=c++17 -stdpar -o base_pstl.exe base_pstl.cpp
#nvc++ -O3 -std=c++17 -stdpar -I./include -o nv_hk_pstl_reference.exe hk_pstl_reference.cpp

export GCC_HOME=/home/astrel/install/spack-0.16.3/opt/spack/linux-ubuntu20.04-skylake/gcc-9.3.0/gcc-11.2.0-wby54stmjhxbuxrgsxwjz3zzr4edo6o4

nvc++ -O2 -std=c++20 --gcc-toolchain=${GCC_HOME} -stdpar=gpu  -gpu=cc86  -gpu=fma -gpu=fastmath -gpu=managed -gpu=unroll -gpu=autocollapse -gpu=loadcache:L1  -I./include -o nv_3dhk_pstl_reference.exe hk_pstl_reference.cpp

nvc++ -O2 -std=c++20 --gcc-toolchain=${GCC_HOME} -stdpar=multicore -I./include -o x86_3dhk_pstl_reference.exe hk_pstl_reference.cpp
