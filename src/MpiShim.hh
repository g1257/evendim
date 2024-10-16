#ifndef EVENDIM_MPI_SHIM_HH
#define EVENDIM_MPI_SHIM_HH
#ifdef USE_MPI
#include "MpiShimActual.hh"
#else
#include "MpiShimSerial.hh"
#endif
#endif
