#include "Algorithm.hpp"
#include "Optimizer.hpp"
#include "PauliOperator.hpp"
#include "xacc.hpp"
#include "xacc_service.hpp"
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

int main(int argc, char** argv)
{

	xacc::set_verbose(true);
	using ProgramType = std::shared_ptr<xacc::CompositeInstruction>;
	using InstructionType = std::shared_ptr<xacc::Instruction>;
	using PauliOperatorType = xacc::quantum::PauliOperator;

	xacc::Initialize(argc, argv);
	auto provider = xacc::getIRProvider("quantum");
	std::vector<InstructionType> instructions;

	// Create X, Ry, CX, and Measure gates
	auto x = provider->createInstruction("X", { 0 });
	instructions.push_back(x);
	auto ry = provider->createInstruction("Ry", { 1 }, { "t0" });
	instructions.push_back(ry);
	auto cx = provider->createInstruction("CNOT", { 1, 0 });
	instructions.push_back(cx);
	auto m0 = provider->createInstruction("Measure", { 0 });
	instructions.push_back(m0);
	auto m1 = provider->createInstruction("Measure", { 1 });
	instructions.push_back(m1);

	// create program
	std::vector<std::string> total_params = { "t0" };
	ProgramType program = provider->createComposite("foo", total_params);

	// Add them to the CompositeInstruction
	program->addInstructions(instructions);

	auto buffer = xacc::qalloc(2);
	auto accelerator = xacc::getAccelerator("qsim");
	auto optimizer = xacc::getOptimizer("nlopt");

	xacc::Observable* pauliOperator = new PauliOperatorType("X0 X1");
	auto vqe = xacc::getService<xacc::Algorithm>("vqe");
	vqe->initialize({ { "ansatz", program },
	                  { "accelerator", accelerator },
	                  { "observable", pauliOperator },
	                  { "optimizer", optimizer } });
	vqe->execute(buffer);

	xacc::HeterogeneousMap extra_data;

	double energy = buffer->getInformation("opt-val").as<double>();

	std::cout << "Energy is " << energy << "\n";
	xacc::Finalize();
}
