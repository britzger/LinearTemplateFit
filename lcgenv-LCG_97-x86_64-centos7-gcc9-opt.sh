# setup new env
export PLATFORM=x86_64-centos7-gcc9-opt
# Using globaly provided gcc from /cvmfs/sft.cern.ch/lcg/contrib
# found package Davix-22716
# found package Boost-d1982
# found package Python-e553a
# found package sqlite-472f2
# found package numpy-4438d
# found package blas-bb5ea
# found package Python-e553a
# found package lapack-bb169
# found package blas-bb5ea
# found package setuptools-43c0e
# found package Python-e553a
# found package libxml2-3501e
# found package srm_ifce-be254
# found package blas-bb5ea
# found package tbb-bbcac
# found package vdt-992df
# found package Vc-d25d4
# found package xrootd-4c7b6
# found package zlib-da225
# found package Python-e553a
# found package libxml2-3501e
# found package fftw-102c2
# found package gl2ps-30f0b
# found package zlib-da225
# found package png-9c2fe
# found package mysql-cf3b2
# found package lz4-9bdfe
# found package zeromq-e6a84
# found package pkg_config-c6baf
# found package libsodium-0b20d
# found package Boost-d1982
# found package libxml2-3501e
# found package msgpackc-c4011
# found package zlib-da225
# found package jemalloc-8154a
# found package bison-9b377
# found package m4-b8d0d
# found package libxml2-3501e
# found package oracle-1837a
# found package zlib-da225
# found package gfal-6fc75
# found package Python-e553a
# found package dcap-cdd28
# found package numpy-4438d
# found package GSL-32fc5
source /cvmfs/sft.cern.ch/lcg/contrib/gcc/9/x86_64-centos7/setup.sh;
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/ROOT/v6.20.02/x86_64-centos7-gcc9-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/ROOT/v6.20.02/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/ROOT/v6.20.02/x86_64-centos7-gcc9-opt/lib/JupyROOT:$LD_LIBRARY_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/ROOT/v6.20.02/x86_64-centos7-gcc9-opt/lib/JsMVA:$LD_LIBRARY_PATH";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Davix/0.7.3/x86_64-centos7-gcc9-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Davix/0.7.3/x86_64-centos7-gcc9-opt/lib64:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Davix/0.7.3/x86_64-centos7-gcc9-opt/lib64/pkgconfig:$PKG_CONFIG_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Boost/1.72.0/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Boost/1.72.0/x86_64-centos7-gcc9-opt/lib/cmake:$LD_LIBRARY_PATH";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Python/2.7.16/x86_64-centos7-gcc9-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Python/2.7.16/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Python/2.7.16/x86_64-centos7-gcc9-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export PYTHONPATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Python/2.7.16/x86_64-centos7-gcc9-opt/lib/python2.7/site-packages:$PYTHONPATH";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/sqlite/3280000/x86_64-centos7-gcc9-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/sqlite/3280000/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/sqlite/3280000/x86_64-centos7-gcc9-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export SQLITE__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/sqlite/3280000/x86_64-centos7-gcc9-opt";
export PYTHON__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Python/2.7.16/x86_64-centos7-gcc9-opt";
cd "/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Python/2.7.16/x86_64-centos7-gcc9-opt"
export PYTHONHOME="${PYTHON__HOME}"
export PYTHONVERSION=`python -c "from __future__ import print_function; import sys; print ('%d.%d' % (sys.version_info.major, sys.version_info.minor))"`
cd - 1>/dev/null # from /cvmfs/sft.cern.ch/lcg/releases/LCG_97/Python/2.7.16/x86_64-centos7-gcc9-opt
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/numpy/1.16.4/x86_64-centos7-gcc9-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/numpy/1.16.4/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export PYTHONPATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/numpy/1.16.4/x86_64-centos7-gcc9-opt/lib/python2.7/site-packages:$PYTHONPATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/blas/0.3.5.openblas/x86_64-centos7-gcc9-opt/lib64:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/blas/0.3.5.openblas/x86_64-centos7-gcc9-opt/lib64/pkgconfig:$PKG_CONFIG_PATH";
export BLAS__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/blas/0.3.5.openblas/x86_64-centos7-gcc9-opt";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/lapack/3.8.0/x86_64-centos7-gcc9-opt/lib64:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/lapack/3.8.0/x86_64-centos7-gcc9-opt/lib64/pkgconfig:$PKG_CONFIG_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/lapack/3.8.0/x86_64-centos7-gcc9-opt/lib64/cmake:$LD_LIBRARY_PATH";
export LAPACK__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/lapack/3.8.0/x86_64-centos7-gcc9-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/setuptools/41.0.0/x86_64-centos7-gcc9-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/setuptools/41.0.0/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export PYTHONPATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/setuptools/41.0.0/x86_64-centos7-gcc9-opt/lib/python2.7/site-packages:$PYTHONPATH";
export SETUPTOOLS__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/setuptools/41.0.0/x86_64-centos7-gcc9-opt";
export NUMPY__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/numpy/1.16.4/x86_64-centos7-gcc9-opt";
export BOOST__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Boost/1.72.0/x86_64-centos7-gcc9-opt";
cd "/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Boost/1.72.0/x86_64-centos7-gcc9-opt"
export CPLUS_INCLUDE_PATH=${BOOST__HOME}/include:$CPLUS_INCLUDE_PATH
export C_INCLUDE_PATH=${BOOST__HOME}/include:$C_INCLUDE_PATH
cd - 1>/dev/null # from /cvmfs/sft.cern.ch/lcg/releases/LCG_97/Boost/1.72.0/x86_64-centos7-gcc9-opt
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/libxml2/2.9.9/x86_64-centos7-gcc9-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/libxml2/2.9.9/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/libxml2/2.9.9/x86_64-centos7-gcc9-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/libxml2/2.9.9/x86_64-centos7-gcc9-opt/lib/cmake:$LD_LIBRARY_PATH";
export LIBXML2__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/libxml2/2.9.9/x86_64-centos7-gcc9-opt";
export DAVIX__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Davix/0.7.3/x86_64-centos7-gcc9-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Grid/srm-ifce/1.13.0-0/x86_64-centos7-gcc9-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Grid/srm-ifce/1.13.0-0/x86_64-centos7-gcc9-opt/lib64:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Grid/srm-ifce/1.13.0-0/x86_64-centos7-gcc9-opt/lib64/pkgconfig:$PKG_CONFIG_PATH";
export SRM_IFCE__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Grid/srm-ifce/1.13.0-0/x86_64-centos7-gcc9-opt";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/tbb/2020_U1/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export TBB__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/tbb/2020_U1/x86_64-centos7-gcc9-opt";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/vdt/0.4.3/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export VDT__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/vdt/0.4.3/x86_64-centos7-gcc9-opt";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Vc/1.4.1/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Vc/1.4.1/x86_64-centos7-gcc9-opt/lib/cmake:$LD_LIBRARY_PATH";
export VC__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Vc/1.4.1/x86_64-centos7-gcc9-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/xrootd/4.11.2/x86_64-centos7-gcc9-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/xrootd/4.11.2/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export PYTHONPATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/xrootd/4.11.2/x86_64-centos7-gcc9-opt/lib/python2.7/site-packages:$PYTHONPATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/xrootd/4.11.2/x86_64-centos7-gcc9-opt/lib64:$LD_LIBRARY_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/zlib/1.2.11/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/zlib/1.2.11/x86_64-centos7-gcc9-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export ZLIB__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/zlib/1.2.11/x86_64-centos7-gcc9-opt";
export XROOTD__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/xrootd/4.11.2/x86_64-centos7-gcc9-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/fftw3/3.3.8/x86_64-centos7-gcc9-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/fftw3/3.3.8/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/fftw3/3.3.8/x86_64-centos7-gcc9-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/fftw3/3.3.8/x86_64-centos7-gcc9-opt/lib/cmake:$LD_LIBRARY_PATH";
export FFTW__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/fftw3/3.3.8/x86_64-centos7-gcc9-opt";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/gl2ps/1.4.0/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/png/1.6.37/x86_64-centos7-gcc9-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/png/1.6.37/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/png/1.6.37/x86_64-centos7-gcc9-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export PNG__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/png/1.6.37/x86_64-centos7-gcc9-opt";
export GL2PS__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/gl2ps/1.4.0/x86_64-centos7-gcc9-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/mysql/10.4.12/x86_64-centos7-gcc9-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/mysql/10.4.12/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/mysql/10.4.12/x86_64-centos7-gcc9-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/mysql/10.4.12/x86_64-centos7-gcc9-opt/lib/plugin:$LD_LIBRARY_PATH";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/lz4/1.9.2/x86_64-centos7-gcc9-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/lz4/1.9.2/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/lz4/1.9.2/x86_64-centos7-gcc9-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export LZ4__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/lz4/1.9.2/x86_64-centos7-gcc9-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/zeromq/4.3.2/x86_64-centos7-gcc9-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/zeromq/4.3.2/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/zeromq/4.3.2/x86_64-centos7-gcc9-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/pkg_config/0.29.2/x86_64-centos7-gcc9-opt/bin:$PATH";
export PKG_CONFIG__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/pkg_config/0.29.2/x86_64-centos7-gcc9-opt";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/libsodium/1.0.18/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/libsodium/1.0.18/x86_64-centos7-gcc9-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export LIBSODIUM__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/libsodium/1.0.18/x86_64-centos7-gcc9-opt";
export ZEROMQ__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/zeromq/4.3.2/x86_64-centos7-gcc9-opt";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/msgpackc/3.2.0/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/msgpackc/3.2.0/x86_64-centos7-gcc9-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/msgpackc/3.2.0/x86_64-centos7-gcc9-opt/lib/cmake:$LD_LIBRARY_PATH";
export MSGPACKC__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/msgpackc/3.2.0/x86_64-centos7-gcc9-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/jemalloc/5.2.1/x86_64-centos7-gcc9-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/jemalloc/5.2.1/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/jemalloc/5.2.1/x86_64-centos7-gcc9-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export JEMALLOC__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/jemalloc/5.2.1/x86_64-centos7-gcc9-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/bison/3.3.2/x86_64-centos7-gcc9-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/bison/3.3.2/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/m4/1.4.18/x86_64-centos7-gcc9-opt/bin:$PATH";
export M4__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/m4/1.4.18/x86_64-centos7-gcc9-opt";
cd "/cvmfs/sft.cern.ch/lcg/releases/LCG_97/m4/1.4.18/x86_64-centos7-gcc9-opt"
export M4=${M4__HOME}/bin/m4
cd - 1>/dev/null # from /cvmfs/sft.cern.ch/lcg/releases/LCG_97/m4/1.4.18/x86_64-centos7-gcc9-opt
export BISON__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/bison/3.3.2/x86_64-centos7-gcc9-opt";
export MYSQL__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/mysql/10.4.12/x86_64-centos7-gcc9-opt";
cd "/cvmfs/sft.cern.ch/lcg/releases/LCG_97/mysql/10.4.12/x86_64-centos7-gcc9-opt"
export MYSQL_HOME=${MYSQL_HOME}
cd - 1>/dev/null # from /cvmfs/sft.cern.ch/lcg/releases/LCG_97/mysql/10.4.12/x86_64-centos7-gcc9-opt
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/oracle/19.3.0.0.0/x86_64-centos7-gcc9-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/oracle/19.3.0.0.0/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/oracle/19.3.0.0.0/x86_64-centos7-gcc9-opt/lib/network:$LD_LIBRARY_PATH";
export ORACLE__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/oracle/19.3.0.0.0/x86_64-centos7-gcc9-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Grid/gfal/1.13.0-0/x86_64-centos7-gcc9-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Grid/gfal/1.13.0-0/x86_64-centos7-gcc9-opt/lib64:$LD_LIBRARY_PATH";
export PYTHONPATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Grid/gfal/1.13.0-0/x86_64-centos7-gcc9-opt/lib64/python2.6/site-packages:$PYTHONPATH";
export GFAL__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Grid/gfal/1.13.0-0/x86_64-centos7-gcc9-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Grid/dcap/2.47.7-1/x86_64-centos7-gcc9-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Grid/dcap/2.47.7-1/x86_64-centos7-gcc9-opt/lib64:$LD_LIBRARY_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Grid/dcap/2.47.7-1/x86_64-centos7-gcc9-opt/lib64/dcap:$LD_LIBRARY_PATH";
export DCAP__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/Grid/dcap/2.47.7-1/x86_64-centos7-gcc9-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/GSL/2.5/x86_64-centos7-gcc9-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/GSL/2.5/x86_64-centos7-gcc9-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/GSL/2.5/x86_64-centos7-gcc9-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export GSL__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/GSL/2.5/x86_64-centos7-gcc9-opt";
export ROOT__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/ROOT/v6.20.02/x86_64-centos7-gcc9-opt";
cd "/cvmfs/sft.cern.ch/lcg/releases/LCG_97/ROOT/v6.20.02/x86_64-centos7-gcc9-opt"
test -s $ROOT__HOME/bin/thisroot.sh && source $ROOT__HOME/bin/thisroot.sh
cd - 1>/dev/null # from /cvmfs/sft.cern.ch/lcg/releases/LCG_97/ROOT/v6.20.02/x86_64-centos7-gcc9-opt
# Using globaly provided gcc from /cvmfs/sft.cern.ch/lcg/contrib
source /cvmfs/sft.cern.ch/lcg/contrib/gcc/9/x86_64-centos7/setup.sh;
export EIGEN__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_97/eigen/3.3.7/x86_64-centos7-gcc9-opt";
