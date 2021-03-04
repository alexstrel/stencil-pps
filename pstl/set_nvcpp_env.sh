
#export JITIFY_OPTIONS="-include complex"

#nvc++
export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/20.9/compilers/bin:$PATH
export LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/20.9/compilers/lib:$LD_LIBRARY_PATH

CDIR=`pwd`

echo $CDIR

makelocalrc -x -d $CDIR  -gcc /opt/gcc-10.2.0/bin/gcc -gpp /opt/gcc-10.2.0/bin/g++ -g77 /opt/gcc-10.2.0/bin/gfortran

export NVLOCALRC=$CDIR/localrc

