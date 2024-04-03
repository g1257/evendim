# PsimagLite support is needed by Gep
LDFLAGS = -L../../PsimagLite/lib -lpsimaglite

CPPFLAGS += -I../../PsimagLite -I../../PsimagLite/src
CPPFLAGS += -I/usr/include/hdf5/serial

# Compiler to use. If using MPI then say mpicxx here (or mpic++)
# and also say -DUSE_MPI below
CXX = g++

# We're using ansi C++
CPPFLAGS += -pedantic -std=c++11

# Enable MPI (you must set the proper
# compiler wrapper under CXX above)
# CPPFLAGS += -DUSE_MPI

# Here add -lpthread if threading is needed and also
# set -DUSE_PTHREADS below
#LDFLAGS += -lpthread

# Enable pthreads
#CPPFLAGS += -DUSE_PTHREADS

# Enable warnings and treat warnings as errors
CPPFLAGS += -Wall -Werror -Wendif-labels

# This disables debugging
#CPPFLAGS += -DNDEBUG

# Optimization level here
CPPFLAGS += -O3

# This enables partial debugging (make sure to comment out previous line)
#CPPFLAGS +=   -g3

# This enables additional debugging
#CPPFLAGS += -D_GLIBCXX_DEBUG -D_GLIBCXX_PROFILE

# This makes the code use long instead of short integers
#CPPFLAGS +=-DUSE_LONG

# This makes the code use float instead of double
#CPPFLAGS += -DUSE_FLOAT

# This enables signals
#CPPFLAGS +=-DUSE_SIGNALS

# This enables gsl support
CPPFLAGS +=-DUSE_GSL
LDFLAGS += -lgsl -lgslcblas

LDFLAGS += -lblas -llapack

CPPFLAGS += -DUSE_BOOST

# This enables the custom allocator (use only for debugging)
#CPPFLAGS += -DUSE_CUSTOM_ALLOCATOR

# Specify the strip command to use or use true to disable
STRIP_COMMAND = true

