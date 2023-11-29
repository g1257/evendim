#ifndef LINEARTREEEXEC_HH
#define LINEARTREEEXEC_HH

#ifdef USE_XACC
#include "../xacc/LinearTreeExecXacc.hh"
#else
#include "LinearTreeExecCPU.hh"
#endif

#endif // LINEARTREEEXEC_HH
