include xacc.make

CXX = g++ 
CPPFLAGS += -Wall -Werror -pedantic #-std=c++11

all: QuantumGEPXacc

QuantumGEPXacc.o: QuantumGEPXacc.cpp
	$(CXX) $(CPPFLAGS) -c QuantumGEPXacc.cpp

QuantumGEPXacc: QuantumGEPXacc.o
	$(CXX) $(LDFLAGS) -o QuantumGEPXacc QuantumGEPXacc.o
