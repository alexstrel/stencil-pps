#!/usr/bin/bash
###
#nvc++ -O0 -g -std=c++17 -stdpar -o base_pstl.exe base_pstl.cpp
#nvc++ -O3 -std=c++17 -stdpar -I./include -o nv_hk_pstl_reference.exe hk_pstl_reference.cpp

nvc++ -O2 -std=c++17  -stdpar  -gpu=cc75  -gpu=fma -gpu=fastmath -I./include -o nv_hk_pstl_reference.exe hk_pstl_reference.cpp
