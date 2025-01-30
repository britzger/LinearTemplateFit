setupATLAS
lsetup "views LCG_104b x86_64-centos7-gcc11-dbg"


function compile(){
    echo "Start compilation..."
    make -S . -B build
    cmake --build build
}
