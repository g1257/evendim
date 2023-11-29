#ifndef EVENDIM_QUANTUMGEPXACC_H_
#define EVENDIM_QUANTUMGEPXACC_H_

#include "AllocatorCpu.h"
#include "Complex.h"
#include "QuantumGEPGate.hh"
#include "xacc.hpp"
#include <string>
#include <vector>

namespace Gep {

template <typename T1, typename T2>
struct EnforceTypesEqual { };

template <typename T>
struct EnforceTypesEqual<T, T> {
	using Type = int;
};

template <typename VectorComplexType, typename RealType_>
class LinearTreeExec {

public:

	using ComplexType = typename VectorComplexType::value_type;
	using RealType = typename PsimagLite::Real<ComplexType>::Type;

	void bogus(EnforceTypesEqual<RealType, RealType_>::Type x) { }

	using VecStringType = std::vector<std::string>;
	using ProgramType = std::shared_ptr<xacc::CompositeInstruction>;
	using InstructionType = std::shared_ptr<xacc::Instruction>;

	LinearTreeExec(const VecStringType& vecStr, SizeType /* threadNum */)
	{
		program_ = createProgram(vecStr);
	}

	RealType energy() const
	{
		auto buffer = xacc::qalloc(2);
		double angle = 0.;
		auto evaled = program_->operator()({ angle });
		auto accelerator = xacc::getAccelerator("tnqvm");
		accelerator->execute(buffer, evaled);
		return buffer->getExpectationValueZ();
	}

private:

	ProgramType createProgram(const VecStringType& circuit)
	{
		// Get the IRProvider and create an
		// empty CompositeInstruction
		auto provider = xacc::getIRProvider("quantum");
		auto program = provider->createComposite("foo", { "t" });
		std::vector<InstructionType> instructions;
		SizeType ngates = circuit.size();
		for (SizeType i = 0; i < ngates; ++i) {
			QuantumGEPGate gate(circuit[i]);
			if (gate.isParametric()) {
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

	ProgramType program_;
};
}
#endif // LINEARTREEEXEC_DUMMY_HH
