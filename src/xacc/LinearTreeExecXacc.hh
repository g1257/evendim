#ifndef EVENDIM_QUANTUMGEPXACC_H_
#define EVENDIM_QUANTUMGEPXACC_H_

#include "../Engine/UnderlyingType.hh"
#include "../Fitness/Hamiltonian.h"
#include "AllocatorCpu.h"
#include "Complex.h"
#include "QuantumGEPGate.hh"
#include "xacc.hpp"
#include <string>
#include <vector>

namespace Gep {

template <typename VecComplexType, typename AnglesType, typename CtorParamType, typename HamiltonianType>
class LinearTreeExec {

public:

	using ComplexType = typename UnderlyingType<VecComplexType>::Type;
	using RealType = typename PsimagLite::Real<ComplexType>::Type;
	using VecStringType = std::vector<std::string>;
	using ProgramType = std::shared_ptr<xacc::CompositeInstruction>;
	using InstructionType = std::shared_ptr<xacc::Instruction>;

	struct HandleType {
		ProgramType program;
		SizeType numberOfBits;
		SizeType threadNum;
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

		circuit2 += circuit;

		ProgramType program = createProgram(circuit);
		SizeType numberOfBits = log2Exact(initVector.size());
		return HandleType{program, numberOfBits, threadNum};
	}

	// Ignore hamiltonian for now and assume it's Z_0
	RealType energy(const HandleType& handle, const HamiltonianType& hamiltonian) const
	{
		if (handle.numberOfBits != hamiltonian.numberOfSites()) {
			throw std::runtime_error("Hamiltonian size incorrect\n");
		}

		auto buffer = xacc::qalloc(handle.numberOfBits);
		double angle = 0.;
		auto evaled = handle.program->operator()({ angle });
		auto accelerator = xacc::getAccelerator("tnqvm");
		accelerator->execute(buffer, evaled);
		return buffer->getExpectationValueZ();
	}

private:

	static ProgramType createProgram(const VecStringType& circuit)
	{		
		// Get the IRProvider and create an
		// empty CompositeInstruction
		auto provider = xacc::getIRProvider("quantum");
		auto program = provider->createComposite("foo", { "t" });
		std::vector<InstructionType> instructions;
		SizeType ngates = circuit.size();
		for (SizeType i = 0; i < ngates; ++i) {
			QuantumGEPGate gate(circuit[i]);
			if (!gate.isParametric()) {
				auto someGate = provider->createInstruction(gate.name(), gate.bits());
				instructions.push_back(someGate);
			}
			else {
				auto someGate = provider->createInstruction(gate.name(), gate.bits(), gate.params());
				instructions.push_back(someGate);
			}
		}

		// Create X, Ry, CX, and Measure gates
		// auto x = provider->createInstruction("X", { 0 });
		// auto ry = provider->createInstruction("Ry", { 1 }, { "t" });
		// auto cx = provider->createInstruction("CNOT", { 1, 0 });
		auto m0 = provider->createInstruction("Measure", { 0 });
		instructions.push_back(m0);

		// Add them to the CompositeInstruction
		program->addInstructions(instructions);
		return program;
	}

	static SizeType findPureState(const VecComplexType& initVector)
	{
		SizeType n = initVector.size();
		bool hasSeenNonZero = false;
		SizeType x = 0;
		for (SizeType i = 0; i < n; ++i) {
			if (std::norm(initVector[i]) > 0) {
				if (hasSeenNonZero) {
					dieVectorNotPure(initVector);
				}

				hasSeenNonZero = true;
				if (std::imag(initVector[i]) != 0) {
					dieVectorNotPure(initVector);
				}

				if (std::abs(std::real(initVector[i]) - 1) < 1e-4) {
					dieVectorNotPure(initVector);
				}

				x = i;
			}
		}

		if (!hasSeenNonZero) {
			throw std::runtime_error("initVector is zero\n");
		}

		return x;
	}

	void pureVectorToXgates(VecStringType& circuit,
	                     const VecComplexType& initVector)

	{
		SizeType x = findPureState(initVector);
		SizeType i = 0;
		while (x > 0) {
			if (x & 1) {
				circuit.push_back("Sx" + ttos(i));
			}

			x >>= 1;
			++i;
		}
	}

	// If n is 2^x, this function returns x
	// Else it throws
	static SizeType log2Exact(SizeType n)
	{
		SizeType x = 0;
		while (n > 0) {
			if (n & 1) {
				break;
			}

			n >>= 1;
			++x;
		}

		SizeType mustBeN = (1 << x);
		if (mustBeN != n) {
			throw std::runtime_error("n is not a power of 2\n");
		}


		return x;
	}

	static dieVectorNotPure(const VecComplexType&)
	{
		throw std::runtime_error("initVector must be pure\n");
	}
};
}
#endif // LINEARTREEEXEC_DUMMY_HH
