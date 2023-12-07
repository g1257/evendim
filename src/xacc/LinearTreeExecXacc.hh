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

	HandleType getHandle(const VecComplexType& initVector, // xacc init vector?
	                     const VecStringType& circuit,
	                     SizeType threadNum) const
	{
		ProgramType program = createProgram(circuit);
		SizeType numberOfBits = log2Exact(initVector.size());
		return HandleType{program, numberOfBits, threadNum};
	}

	// Ignore hamiltonian for now and assume it's Z_0
	RealType energy(const HandleType& handle, const HamiltonianType& hamiltonian) const
	{
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
};
}
#endif // LINEARTREEEXEC_DUMMY_HH
