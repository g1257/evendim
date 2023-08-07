include xacc.make

CXX = g++ 
CPPFLAGS += -Wall  -pedantic #-std=c++11

all: QuantumGEPXacc libQuantumGEP.so

QuantumGEPXacc.o: QuantumGEPXacc.cpp
	$(CXX) $(CPPFLAGS) -c QuantumGEPXacc.cpp

QuantumGEPXacc: QuantumGEPXacc.o
	$(CXX) $(LDFLAGS) -o QuantumGEPXacc QuantumGEPXacc.o

QuantumGEP_Plugin.o: QuantumGEP_Plugin.cpp
	$(CXX) $(CPPFLAGS) -c -fPIC -o QuantumGEP_Plugin.o QuantumGEP_Plugin.cpp

libQuantumGEP.so: QuantumGEP_Plugin.o
	$(CXX) -shared -o libQuantumGEP.so QuantumGEP_Plugin.o
clean:
	rm -f *.o QuantumGEPXacc *so
