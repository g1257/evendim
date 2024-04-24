#include "Algorithm.hpp"
#include "Optimizer.hpp"
#include "PauliOperator.hpp"
#include "xacc.hpp"
#include "xacc_service.hpp"
#include <iomanip>
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
	vqe->initialize({
	    { "ansatz", program },
	    { "accelerator", accelerator },
	    { "observable", pauliOperator },
	    //{ "optimizer", optimizer }
	});

	auto energy = vqe->execute(buffer, { 1.0 })[0];

	// double energy = buffer->getInformation("opt-val").as<double>();

	std::cout << std::setprecision(12) << "Energy is " << energy << "\n";
	xacc::Finalize();
}
