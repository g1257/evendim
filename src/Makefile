include Config.make
CPPFLAGS += -I../../PsimagLite -I../../PsimagLite/src -IEngine
# Comment out below in order to run with XACC

all: gep2 quantumGep quantumIndividual graphGen isGraphConnected testPauliMatrixExpand testSchwinger test_redirect_output

quantumStatics.o: quantumStatics.cpp  Makefile   Config.make
	$(CXX) $(CPPFLAGS) -c quantumStatics.cpp

quantumGep.o: quantumGep.cpp  Makefile   Config.make
	$(CXX) $(CPPFLAGS) -c quantumGep.cpp

quantumGep: quantumGep.o quantumStatics.o
	$(CXX) -o  quantumGep quantumGep.o quantumStatics.o $(LDFLAGS) $(CPPFLAGS)
	$(STRIP_COMMAND) quantumGep

gep2.o: gep2.cpp  Makefile   Config.make
	$(CXX) $(CPPFLAGS) -c gep2.cpp

gep2: gep2.o
	$(CXX) -o  gep2 gep2.o $(LDFLAGS) $(CPPFLAGS)
	$(STRIP_COMMAND) gep2

testSchwinger.o: testSchwinger.cpp Makefile   Config.make
	$(CXX) $(CPPFLAGS) -c testSchwinger.cpp

testSchwinger: testSchwinger.o
	$(CXX) -o  testSchwinger testSchwinger.o $(LDFLAGS) 

graphGen: graphGen.o
	$(CXX) -o  graphGen graphGen.o $(LDFLAGS) $(CPPFLAGS)
	$(STRIP_COMMAND) graphGen

graphGen.o: graphGen.cpp  Makefile   Config.make
	$(CXX) $(CPPFLAGS) -c graphGen.cpp

isGraphConnected: isGraphConnected.o
	$(CXX) -o  isGraphConnected isGraphConnected.o $(LDFLAGS) $(CPPFLAGS)
	$(STRIP_COMMAND) isGraphConnected

isGraphConnected.o: isGraphConnected.cpp Makefile   Config.make
	$(CXX) $(CPPFLAGS) -c  isGraphConnected.cpp

quantumIndividual.o: quantumIndividual.cpp Makefile   Config.make
	$(CXX) $(CPPFLAGS) -c  quantumIndividual.cpp

quantumIndividual: quantumIndividual.o quantumStatics.o
	$(CXX) -o quantumIndividual quantumIndividual.o quantumStatics.o $(LDFLAGS) $(CPPFLAGS)
	$(STRIP_COMMAND) quantumIndividual

testPauliMatrixExpand: testPauliMatrixExpand.o
	$(CXX) -o testPauliMatrixExpand testPauliMatrixExpand.o $(LDFLAGS)

testPauliMatrixExpand.o: testPauliMatrixExpand.cpp xacc/ToPauliMatrices.hh
	$(CXX) $(CPPFLAGS) -c  testPauliMatrixExpand.cpp

../../PsimagLite/lib/libpsimaglite.a:
	$(MAKE) -f Makefile -C ../../PsimagLite/lib/

Makefile.dep: gep2.cpp
	$(CXX) $(CPPFLAGS) -MM gep2.cpp  quantumGep.cpp quantumIndividual.cpp graphGen.cpp isGraphConnected.cpp testSchwinger.cpp > Makefile.dep

clean: Makefile.dep
	rm -f core* gep2 *.o *.dep

include Makefile.dep
