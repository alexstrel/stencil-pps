CXX = g++
NVCXX = nvc++

NVARCH = cc70

#note -O2 option led to a substantial performance degradation for (gnu) x86 backend 
CXXCOPT_SKL  = -O3 -mavx512f -std=c++17 -fopenmp 
CXXCOPT_HSW  = -O3 -mavx2 -mfma -std=c++17 -fopenmp
LDFLAG =  -lm  -fopenmp -L/opt/intel/tbb-gnu9.3/lib -ltbb

NVCXXOPT = -O2 -std=c++17 -stdpar -gpu=$(NVARCH) -I.

INC = -I./include

all: pstl_stencil_test.avx512 pstl_stencil_test.avx2 pstl_stencil_test.$(NVARCH)

pstl_stencil_test.avx512: hk_pstl_reference_avx512.o
	$(CXX) -o  $@  $^  $(LDFLAG)

hk_pstl_reference_avx512.o : hk_pstl_reference.cpp
	$(CXX) $(CXXCOPT_SKL) -c -o  $@  $?  $(INC)

pstl_stencil_test.avx2: hk_pstl_reference_avx2.o
	$(CXX) -o  $@  $^  $(LDFLAG)

hk_pstl_reference_avx2.o : hk_pstl_reference.cpp
	$(CXX) $(CXXCOPT_HSW) -c -o  $@  $?  $(INC)

pstl_stencil_test.$(NVARCH): hk_pstl_reference_$(NVARCH).o
	$(NVCXX) -o  $@  $^

hk_pstl_reference_$(NVARCH).o : hk_pstl_reference.cpp
	$(NVCXX) $(NVCXXOPT) -c -o  $@  $?  $(INC)

clean:
	rm *.o *.avx512 *.avx2 *.$(NVARCH)

.PHONY:	clean
