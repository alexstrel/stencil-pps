
#nvc++ -std=c++20  -stdpar=gpu -gpu=cc75 -gpu=managed -gpu=fma -gpu=fastmath -gpu=autocollapse -gpu=loadcache:L1 -gpu=unroll -I./include -o test.exe U1D2-stdpar.cpp
#
nvc++ -std=c++20  -stdpar=gpu -gpu=cc75 -gpu=managed -gpu=fma -gpu=fastmath -gpu=autocollapse -gpu=loadcache:L1 -gpu=unroll -I./include -o test2.exe U1D2-stdpar-2.cpp
