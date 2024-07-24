#include "../Engine/NodeFactory.h"
#include "HamiltonianXacc.hh"
#include "LinearTreeExecXacc.hh"
#include "XaccBackendActual.hh"
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

std::string implodeVecString(const std::vector<std::string>& vecStr)
{
	const char* const delim = ", ";

	std::ostringstream imploded;
	std::copy(vecStr.begin(), vecStr.end(), std::ostream_iterator<std::string>(imploded, delim));
	return imploded.str();
}

int main(int argc, char** argv)
{
	// This call xacc::init in its ctor and xacc:fin in its dtor
	Gep::XaccBackend xaccBackend(argc, argv);

	using DummyUnusedType = int;
	using HamiltonianType = Gep::Hamiltonian<std::complex<double>>;
	using LinearTreeExecType = Gep::LinearTreeExec<std::vector<std::complex<double>>, double, DummyUnusedType>;
	using HandleType = LinearTreeExecType::HandleType;

	// here is the circuit
	typename LinearTreeExecType::VecStringType mycircuit { "Sx0", "Ry1:0" };

	// here is the initial state
	std::vector<std::complex<double>> initVector(4);
	initVector[0] = 1;

	DummyUnusedType dummy = 0;
	LinearTreeExecType linearTreeExec(dummy);

	// create xacc program and store in handle
	constexpr int threadNum = 0; // no parallelization for now
	HandleType handle = linearTreeExec.getHandle(initVector, mycircuit, threadNum);

	// Does energy = <0|C H C |0>, with H = Z_0
	constexpr SizeType numberOfThreads = 1;
	constexpr SizeType sites = 2;
	HamiltonianType hamiltonian("Sx0 * Sx0", sites, numberOfThreads, "tnqvm");
	constexpr bool useXaccOptimizer = true;
	double energy = linearTreeExec.energy(handle, hamiltonian, useXaccOptimizer);
	std::cout << "Circuit is " << implodeVecString(mycircuit) << "\n";
	std::cout << "Energy is " << energy << "\n";
}
