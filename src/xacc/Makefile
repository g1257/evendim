include xacc.make

CXX = g++ -g3
CPPFLAGS += -Wall  -pedantic -I../../../PsimagLite/src -I../../../PsimagLite  #-std=c++11
CPPFLAGS += -I/usr/include/hdf5/serial
LDFLAGS += -L ../../../PsimagLite/lib -lpsimaglite

all: testXacc mweOptimizer mweNoAngles

#QuantumGEPXacc libQuantumGEP.so

mweNoAngles.o: mweNoAngles.cpp
	$(CXX) $(CPPFLAGS) -c mweNoAngles.cpp

mweNoAngles: mweNoAngles.o
	$(CXX) -o mweNoAngles mweNoAngles.o $(LDFLAGS)

mweOptimizer.o: mweOptimizer.cpp
	$(CXX) $(CPPFLAGS) -c mweOptimizer.cpp

mweOptimizer: mweOptimizer.o
	$(CXX)  -o mweOptimizer mweOptimizer.o $(LDFLAGS)

testXacc.o: testXacc.cpp
	$(CXX) $(CPPFLAGS) -c testXacc.cpp

testXacc: testXacc.o quantumStatics.o
	$(CXX) testXacc.o quantumStatics.o -o testXacc $(LDFLAGS)

quantumStatics.o: ../quantumStatics.cpp  Makefile
	        $(CXX) $(CPPFLAGS) -c ../quantumStatics.cpp

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
