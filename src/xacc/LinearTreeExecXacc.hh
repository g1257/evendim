#ifndef EVENDIM_QUANTUMGEPXACC_H_
#define EVENDIM_QUANTUMGEPXACC_H_

#include "../Engine/UnderlyingType.hh"
#include "AllocatorCpu.h"
#include "Complex.h"
#include "HamiltonianXacc.hh"
#include "QuantumGEPGate.hh"
#include "xacc.hpp"
#include <string>
#include <vector>

namespace Gep {

template <typename T1, typename T2>
struct TypesEqual {
	static const bool value = false;
};

template <typename T>
struct TypesEqual<T, T> {
	static const bool value = true;
};

template <bool b, typename T1, typename T2>
struct FirstOrSecondType {
	using type = T2;
};

template <typename T1, typename T2>
struct FirstOrSecondType<true, T1, T2> {
	using type = T1;
};

template <typename VecComplexType, typename AnglesType, typename CtorParamType>
class LinearTreeExec {

public:

	static constexpr bool HAS_XACC = true;

	using ComplexType = typename UnderlyingType<VecComplexType>::Type;
	using RealType = typename PsimagLite::Real<ComplexType>::Type;
	using VecStringType = std::vector<std::string>;
	using ProgramType = std::shared_ptr<xacc::CompositeInstruction>;
	using InstructionType = std::shared_ptr<xacc::Instruction>;
	using ProviderType = std::shared_ptr<xacc::IRProvider>;
	using BogusFirstType = typename FirstOrSecondType<TypesEqual<VecComplexType, std::vector<ComplexType>>::value, int*, double*>::type;
	using BogusSecondType = typename FirstOrSecondType<!TypesEqual<VecComplexType, std::vector<ComplexType>>::value, int*, double*>::type;
	using HamiltonianType = Hamiltonian<ComplexType>;

	struct HandleType {
		ProgramType program;
		std::vector<double> angles;
		SizeType numberOfBits = 0;
		SizeType threadNum = 0;
	};

	// Ctor not needed in the xacc version of LinearTreeExec
	LinearTreeExec(const CtorParamType& ctorParam)
	{
	}

	HandleType getHandle(const VecComplexType& initVector,
	                     const VecStringType& circuit,
	                     SizeType threadNum) const
	{
		VecStringType circuit2;
		pureVectorToXgates(circuit2, initVector);

		// Carefully insert incoming circuit but not its junk DNA (if any)
		unsigned int ncircuit = circuit.size();
		for (unsigned int i = 0; i < ncircuit; ++i) {
			if (circuit[i] == "0")
				break;
			circuit2.push_back(circuit[i]);
		}

		// Add identity gate if circuit is empty
		if (circuit2.size() == 0) {
			circuit2.push_back("I0");
		}

		SizeType numberOfBits = log2Exact(initVector.size());

		auto provider = xacc::getIRProvider("quantum");
		std::pair<ProgramType, std::vector<double>> programAndNparams = createProgram(circuit2, provider);
		// printCircuit(circuit2, std::cout);
		return HandleType { programAndNparams.first, programAndNparams.second, numberOfBits, threadNum };
	}

	RealType energy(HandleType& handle, const HamiltonianType& hamiltonian, bool useXaccOptimizer) const
	{
		if (handle.numberOfBits != hamiltonian.numberOfSites()) {
			throw std::runtime_error("Hamiltonian size incorrect\n");
		}

		return hamiltonian.energyXACC(handle.program, handle.angles, useXaccOptimizer);
	}

	void fillAngles(std::vector<double>& angles, const HandleType& handle) const
	{
		angles = handle.angles;
	}

private:

	static std::pair<ProgramType, std::vector<double>> createProgram(const VecStringType& circuit, const ProviderType& provider)
	{
		// Get the IRProvider and create an
		// empty CompositeInstruction
		std::vector<InstructionType> instructions;
		SizeType ngates = circuit.size();
		std::vector<double> angles;
		std::vector<std::string> total_params;
		for (SizeType i = 0; i < ngates; ++i) {
			if (circuit[i] == "0")
				break;
			QuantumGEPGate gate(circuit[i]);

			SizeType nparams = gate.numberOfParams();
			if (nparams == 0) {
				auto someGate = provider->createInstruction(gate.name(), gate.bits());
				instructions.push_back(someGate);
			}
			else {
				if (nparams != 1) {
					throw std::runtime_error(std::string(__FILE__) + " I can only deal with 1 param for now\n");
				}

				// std::cerr << "parametric " << gate.name() << "\n";
				std::string params = "t" + ttos(angles.size());
				total_params.push_back(params);
				auto someGate = provider->createInstruction(gate.name(), gate.bits(), { params });
				instructions.push_back(someGate);
				angles.push_back(gate.param(0));
			}
		}

		// Create X, Ry, CX, and Measure gates
		// auto x = provider->createInstruction("X", { 0 });
		// auto ry = provider->createInstruction("Ry", { 1 }, { "t" });
		// auto cx = provider->createInstruction("CNOT", { 1, 0 });
		// auto m0 = provider->createInstruction("Measure", { 0 });
		// instructions.push_back(m0);

		// create program
		auto program = provider->createComposite("foo", total_params);

		assert(total_params.size() == angles.size());
		// Add them to the CompositeInstruction
		program->addInstructions(instructions);
		return std::pair<ProgramType, std::vector<double>>(program, angles);
	}

	static void pureVectorToXgates(VecStringType& circuit,
	                               const VecComplexType& initVector)

	{
		for (SizeType j = 0; j < initVector.nonZeros(); ++j) {
			SizeType i = initVector.index(j);
			circuit.push_back("Sx" + ttos(i));
		}
	}

	// If n is 2^x, this function returns x
	// Else it throws
	static SizeType log2Exact(SizeType nn)
	{
		SizeType n = nn;
		SizeType x = 0;
		while (n > 0) {
			if (n & 1) {
				break;
			}

			n >>= 1;
			++x;
		}

		SizeType mustBeN = (1 << x);
		if (mustBeN != nn) {
			throw std::runtime_error("n is not a power of 2\n");
		}

		return x;
	}

	static void dieVectorNotPure(const std::vector<ComplexType>&, const std::string& msg)
	{
		throw std::runtime_error("initVector must be pure " + msg + "\n");
	}

	// unused for now
	static void printCircuit(const VecStringType& circuit, std::ostream& os)
	{
		os << "-----------------------\n";
		unsigned int n = circuit.size();
		os << n << "\n";
		for (unsigned int i = 0; i < n; ++i) {
			os << circuit[i] << " ";
		}

		os << "\n-----------------------\n\n";
	}
};
}
#endif // LINEARTREEEXEC_DUMMY_HH
