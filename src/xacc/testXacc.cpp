#include "QuantumGEPXacc.hh"

int main(int argc, char** argv)
{
	Gep::QuantumGEPXacc quantumGEPXacc(argc, argv);

	Gep::QuantumGEPXacc::VecStringType mycircuit = { "Sx0" };

	quantumGEPXacc.testProgram(mycircuit);
}
