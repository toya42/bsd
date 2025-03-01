cmake_minimum_required(VERSION 3.16)

enable_language(Fortran)

project(bayesian-spectral-deconvolution-post Fortran)

set(EXECUTABLE post.exe)

if(CMAKE_Fortran_COMPILER_ID MATCHES IntelLLVM)
   set(CMAKE_Fortran_FLAGS_DEBUG   "-O0                -warn all -standf08     \
                                    -warn declarations -CB       -check uninit \
                                    -check all         -fpe0     -traceback    \
                                    -g                 -fpp      -Ddebug -lpthread -ldl")
   set(CMAKE_Fortran_FLAGS_PROFILE "-O3 -g -ipo -static -xHost")
   set(CMAKE_Fortran_FLAGS_FAST    "-O3 -fpp -march=native -mtune=native  -lpthread -ldl ")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
   set(CMAKE_Fortran_FLAGS_DEBUG   "-O0            -Wall -pedantic  -std=f2008                        \
                                    -fbounds-check -Wuninitialized  -ffpe-trap=invalid,zero,overflow  \
                                    -fbacktrace    -g -Ddebug -lpthread -ldl -Wl,--no-as-needed")
   set(CMAKE_Fortran_FLAGS_PROFILE "-O3 -march=narive -mtune=native -p -lpthread -ldl -Wl,--no-as-needed")
   set(CMAKE_Fortran_FLAGS_FAST    "-O3 -march=native -mtune=native -lpthread -ldl -Wl,--no-as-needed")
endif()


add_executable(${EXECUTABLE}
  precision.f90
  GlobalData.f90
  DataInput.f90
  DataOutput.f90
  FaddeevaTable.f90
  ModelFunctions.f90
  PosteriorStatistics.f90
  main.f90
  )
