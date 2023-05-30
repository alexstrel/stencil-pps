nvc++ -std=c++23  -stdpar=gpu -gpu=cc75 -gpu=managed -gpu=fma -gpu=fastmath -gpu=autocollapse -gpu=loadcache:L1 -gpu=unroll -I./include -o test.exe U1D2-stdpar.cpp

