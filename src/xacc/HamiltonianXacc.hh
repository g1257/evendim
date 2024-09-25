#ifndef HAMILTONIAN_XACC_H
#define HAMILTONIAN_XACC_H

#include <type_traits>

#include "../Engine/InputCheck.h"
#include "../Fitness/IsingGraph.hh"
#include "Algorithm.hpp"
#include "InputNg.h"
#include "Optimizer.hpp"
#include "PauliOperator.hpp"
#include "PsimagLite.h"
#include "QuantumGEPGate.hh"
#include "Stream.hpp"
#include "ToPauliMatrices.hh"
#include "Vector.h"
#include "xacc.hpp"
#include "xacc_service.hpp"

namespace Gep {

template <typename ComplexType>
class Hamiltonian {

public:

	using VectorType = std::vector<ComplexType>;
	using RealType = typename PsimagLite::Real<ComplexType>::Type;
	using PauliOperatorType = xacc::quantum::PauliOperator;
	using BufferType = std::shared_ptr<xacc::AcceleratorBuffer>;
	using AcceleratorType = std::shared_ptr<xacc::Accelerator>;
	using InputNgType = PsimagLite::InputNg<InputCheck>;
	using VectorStringType = std::vector<std::string>;
	using ProgramType = std::shared_ptr<xacc::CompositeInstruction>;
	using IsingGraphType = IsingGraph<ComplexType>;

	Hamiltonian(typename InputNgType::Readable& io, SizeType /* numberOfThreads */)
	    : bits_(0)
	    , pauliOperator_(nullptr)
	    , verbose_("false")
	{
		io.readline(bits_, "NumberOfBits="); // == number of "sites"
		io.readline(ham_, "Hamiltonian=");
		io.readline(accel_, "Accelerator=");
		if (ham_.substr(0, 5) == "file:") {
			unimplemented("Hamiltonian=file:\n");
		}

		if (ham_ == "zxz" || ham_ == "xx") {
			throw std::runtime_error("Hamiltonian=zxz or xx unsupported in XACC mode\n");
		}

		if (ham_ == "IsingGraph" || ham_ == "zz") {
			PsimagLite::String graphFile = "zz";
			if (ham_ == "IsingGraph") {
				io.readline(graphFile, "GraphFile=");
				graphFile = "file:" + graphFile;
			}

			bool hasPeriodic = false;
			try {
				int tmp = 0;
				io.readline(tmp, "HamiltonianIsPeriodic=");
				hasPeriodic = true;
			}
			catch (std::exception&) {
			}
			if (hasPeriodic) {
				err("IsingGraph: HamiltonianIsPeriodic= line unsupported\n");
			}

			RealType coupling = 1;
			try {
				io.readline(coupling, "HamiltonianCoupling=");
			}
			catch (std::exception&) {
			}

			IsingGraphType ising_graph(bits_, coupling, false, graphFile);
			ham_ = ising_graph.buildExpression();
		}

		try {
			io.readline(verbose_, "XaccVerbosityLevel=");
		}
		catch (...) {
		}

		fromExpression(ham_);
	}

	Hamiltonian(const std::string& expression, SizeType bits, SizeType /* numberOfThreads */, const std::string& accel)
	    : bits_(bits)
	    , ham_(expression)
	    , accel_(accel)
	    , pauliOperator_(nullptr)
	{
		fromExpression(expression);
	}

	~Hamiltonian()
	{
		delete pauliOperator_;
		pauliOperator_ = nullptr;
	}

	std::vector<ProgramType> observe(ProgramType function) const
	{
		return pauliOperator_->observe(function);
	}

	double energyXACC(ProgramType program,
	                  std::vector<double>& angles,
	                  bool useXaccOptimizer) const
	{
		auto buffer = xacc::qalloc(bits_);
		auto accelerator = xacc::getAccelerator(accel_);
		if (useXaccOptimizer && angles.size() > 0) {
			return energyOptimizeAngles(angles, buffer, program, accelerator);
		}
		else {
			return energyFixedAngles(angles, buffer, program, accelerator);
		}
	}

	template <typename SomeType>
	double energy(const SomeType& y, SizeType threadNum) const
	{
		throw std::runtime_error("energy(y, thread) must not be called from XACC\n");
	}

	SizeType numberOfSites() const
	{
		return bits_;
	}

private:

	Hamiltonian(const Hamiltonian&) = delete;

	Hamiltonian& operator=(const Hamiltonian&) = delete;

	double energyFixedAngles(const std::vector<double>& angles,
	                         BufferType buffer,
	                         ProgramType program,
	                         AcceleratorType accelerator) const

	{
		auto vqe = xacc::getService<xacc::Algorithm>("vqe");
		vqe->initialize({ { "ansatz", program },
		                  { "accelerator", accelerator },
		                  { "observable", pauliOperator_ } });
		auto vec = vqe->execute(buffer, angles);
		// return buffer->getInformation("opt-val").as<double>();

		bool flag = (verbose_ == "full" || verbose_ == "debug");

		if (flag) {
			std::cerr << "Fixed angles for program\n";
		}

		assert(vec.size() > 0);
		return vec[0];
	}

	/* Here the initial angles for the optimizer are zero */
	double energyOptimizeAngles(std::vector<double>& angles,
	                            BufferType buffer,
	                            ProgramType program,
	                            AcceleratorType accelerator) const
	{
		auto optimizer = xacc::getOptimizer("nlopt");

		auto vqe = xacc::getService<xacc::Algorithm>("vqe");
		vqe->initialize({ { "ansatz", program }, { "accelerator", accelerator }, { "observable", pauliOperator_ }, { "cache-measurement-basis", true }, { "optimizer", optimizer } });

		bool flag = (verbose_ == "debug");
		if (verbose_ == "full") {
			xacc::set_verbose(true);
			flag = true;
		}

		vqe->execute(buffer);
		angles = buffer->getInformation("opt-params").as<std::vector<double>>();

		const std::string& program_as_string = program->toString();
		if (flag) {
			std::cerr << "Non fixed angles for the following program\n";
			std::cerr << program_as_string << "\n";
			std::cerr << "angles=";
			xacc::operator<<(std::cerr, angles);
			std::cerr << "\n";
		}

		return buffer->getInformation("opt-val").as<double>();
	}

	void fromExpression(const std::string& str)
	{
		std::string paulis = toPaulis(removeAllSpaces(str));

		std::cerr << "Asumming Hamiltonian Expression (XACC) " << paulis << "\n";
		pauliOperator_ = new PauliOperatorType(paulis);
	}

	static std::string removeAllSpaces(const std::string& str)
	{
		std::string str2;
		for (std::string::const_iterator it = str.begin(); it != str.end(); ++it) {
			if (*it == ' ' || *it == '\t')
				continue;
			str2 += *it;
		}

		return str2;
	}

	static void unimplemented(const std::string& msg)
	{
		throw std::runtime_error("XACC Backend: unimplemented: " + msg + "\n");
	}

	static std::string toPaulis(const std::string& str)
	{
		// split +
		VectorStringType terms;
		PsimagLite::split(terms, str, "+");

		std::string paulis;
		for (SizeType i = 0; i < terms.size(); ++i) {
			if (i > 0) {
				paulis += " + ";
			}
			std::string term = termToPauli(terms[i]);
			paulis += term;
		}

		return paulis;
	}

	static std::string termToPauli(const std::string& term)
	{
		// split *
		VectorStringType factors;
		PsimagLite::split(factors, term, "*");
		std::string paulis;
		for (SizeType i = 0; i < factors.size(); ++i) {
			std::string factor = factorToPauli(factors[i]);
			paulis += std::string(" ");
			paulis += factor;
		}

		return paulis;
	}

	static std::string factorToPauli(const std::string& factor)
	{
		if (isNumeric(factor))
			return factor;

		return pauliExpansion(factor);
	}

	// All characters are digits or .
	static bool isNumeric(const std::string& str)
	{
		for (std::string::const_iterator it = str.begin(); it != str.end(); ++it) {
			if (*it == '.' || *it == '+' || *it == '-')
				continue;
			if (std::isdigit(*it))
				continue;
			return false;
		}

		return true;
	}

	static std::string pauliExpansion(const std::string& str)
	{
		QuantumGEPGate gate(str);
		std::string name = gate.name();
		if (gate.bits().size() != 1) {
			err("pauliExpansion: Only one-bit gates supported here\n");
		}

		return name + ttos(gate.bits()[0]);
	}

	SizeType bits_;
	std::string ham_;
	std::string accel_;
	std::string verbose_;
	xacc::Observable* pauliOperator_;
};

}
#endif // HAMILTONIAN_XACC_H
