#!/usr/bin/bash
###
export PSTL_USAGE_WARNINGS=1
export ONEDPL_USE_DPCPP_BACKEND=1

GCC_HOME=/opt/spack-0.16.2/opt/spack/linux-ubuntu20.04-icelake/gcc-9.3.0/gcc-11.2.0-r5szlxwm6grctbo6q7mx47gito2mtnva

dpcpp -std=c++20 -O3 --gcc-toolchain= -I./include ./hk_dpl_reference.cpp -o test-dpl.exe 

