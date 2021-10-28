include Config.make
CPPFLAGS += -I../../PsimagLite -I../../PsimagLite/src -IEngine
all: gep2 quantumGep

quantumGep.o: quantumGep.cpp  Makefile   Config.make
	$(CXX) $(CPPFLAGS) -c quantumGep.cpp

quantumGep: quantumGep.o
	$(CXX) -o  quantumGep quantumGep.o $(LDFLAGS) $(CPPFLAGS)
	$(STRIP_COMMAND) quantumGep

gep2.o: gep2.cpp  Makefile   Config.make
	$(CXX) $(CPPFLAGS) -c gep2.cpp

gep2: gep2.o
	$(CXX) -o  gep2 gep2.o $(LDFLAGS) $(CPPFLAGS)
	$(STRIP_COMMAND) gep2

../../PsimagLite/lib/libpsimaglite.a:
	$(MAKE) -f Makefile -C ../../PsimagLite/lib/

Makefile.dep: gep2.cpp
	$(CXX) $(CPPFLAGS) -MM gep2.cpp  quantumGep.cpp > Makefile.dep

clean: Makefile.dep
	rm -f core* gep2 *.o *.dep

include Makefile.dep
