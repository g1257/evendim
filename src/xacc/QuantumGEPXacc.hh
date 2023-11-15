#ifndef EVENDIM_QUANTUMGEPXACC_H_
#define EVENDIM_QUANTUMGEPXACC_H_

#include "AllocatorCpu.h"
#include "xacc.hpp"
#include <string>
#include <vector>
#include "QuantumGEPGate.hh"

namespace Gep
{

class QuantumGEPXacc
{

public:

	using VecStringType = std::vector<std::string>;
	using ProgramType = std::shared_ptr<xacc::CompositeInstruction>;
	using InstructionType = std::shared_ptr<xacc::Instruction>;

	QuantumGEPXacc(int argc, char* argv[])
	{
		xacc::Initialize(argc, argv);
		// Get reference to the Accelerator
	}

	~QuantumGEPXacc()
	{
		xacc::Finalize();
	}

	void testProgram(const VecStringType& circuit)
	{
		ProgramType program = createProgram(circuit);
		test(program);
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
			} else {
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

	void test(const ProgramType& program)
	{
		auto accelerator = xacc::getAccelerator("tnqvm");
		// Loop over [-pi, pi] and compute <Z0>
		auto angles = xacc::linspace(-xacc::constants::pi, xacc::constants::pi, 20);
		for (auto& a : angles) {
			auto buffer = xacc::qalloc(2);
			auto evaled = program->operator()({ a });
			accelerator->execute(buffer, evaled);
			std::cout << "<Z0>(" << a << ") = " << buffer->getExpectationValueZ()
				  << "\n";
		}
	}
};
}
#endif
