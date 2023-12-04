include xacc.make

CXX = g++ 
CPPFLAGS += -Wall  -pedantic -I../../../PsimagLite/src -I../../../PsimagLite  #-std=c++11

all: testXacc 
#QuantumGEPXacc libQuantumGEP.so

testXacc.o: testXacc.cpp
	$(CXX) $(CPPFLAGS) -c testXacc.cpp

testXacc: testXacc.o
	$(CXX) $(LDFLAGS) testXacc.o -o testXacc

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
