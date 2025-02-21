cmake_minimum_required(VERSION 3.16)

enable_language(Fortran)

project(bayesian-spectral-deconvolution Fortran)

set(EXECUTABLE mcmc.exe)

add_executable(${EXECUTABLE}
  precision.f90
  GlobalData.f90
  DataInput.f90
  test.f90
  )