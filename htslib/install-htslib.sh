#!/bin/bash

pdir=$(pwd)/htslib

FILE=${pdir}/htslib-1.19.1.tar.bz2
SRC=${pdir}/htslib-1.19.1
dir=${pdir}/htslib

mkdir -p ${dir}
mkdir -p ${SRC}

echo ${FILE}
echo ${SRC}
echo ${dir}

tar xjvf ${FILE} -C ${pdir}
if [ $? -ne 0 ]; then
    echo "fail at uncompressing ${FILE}"
    exit 1
fi

curr_dir=$(pwd)
echo ${SRC}
cd ${SRC} && ./configure --disable-lzma --disable-libcurl --prefix=${dir} && make && make install && cd ${curr_dir}
