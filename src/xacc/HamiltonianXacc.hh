#ifndef HAMILTONIAN_XACC_H
#define HAMILTONIAN_XACC_H

#include "../Engine/InputCheck.h"
#include "InputNg.h"
#include "PauliOperator.hpp"
#include "QuantumGEPGate.hh"
#include "ToPauliMatrices.hh"
#include "xacc.hpp"

namespace Gep {

template <typename ComplexType>
class Hamiltonian {

public:

	using VectorType = std::vector<ComplexType>;
	using RealType = typename PsimagLite::Real<ComplexType>::Type;
	using PauliOperatorType = xacc::quantum::PauliOperator;
	using InputNgType = PsimagLite::InputNg<InputCheck>;
	using VectorStringType = std::vector<std::string>;
	using ProgramType = std::shared_ptr<xacc::CompositeInstruction>;

	Hamiltonian(typename InputNgType::Readable& io, SizeType /* numberOfThreads */)
	    : bits_(0)
	    , pauliOperator_(nullptr)
	{
		io.readline(bits_, "NumberOfBits="); // == number of "sites"
		io.readline(ham_, "Hamiltonian=");
		if (ham_.substr(0, 5) == "file:") {
			unimplemented("Hamiltonian=file:\n");
		}

		if (ham_ == "IsingGraph" || ham_ == "zz") {
			throw std::runtime_error("Hamiltonian=IsingGraph or zz\n");
		}

		if (ham_ == "zxz" || ham_ == "xx") {
			throw std::runtime_error("Hamiltonian=zxz or xx\n");
		}

		fromExpression(ham_);
	}

	Hamiltonian(const std::string& expression, SizeType bits, SizeType /* numberOfThreads */)
	    : bits_(bits)
	    , ham_(expression)
	    , pauliOperator_(nullptr)
	{
		fromExpression(expression);
	}

	std::vector<ProgramType> observe(ProgramType function) const
	{
		return pauliOperator_->observe(function);
	}

	double postProcess(std::shared_ptr<xacc::AcceleratorBuffer> buffer) const
	{
		xacc::HeterogeneousMap extra_data;
		return pauliOperator_->postProcess(buffer, xacc::Observable::PostProcessingTask::EXP_VAL_CALC, extra_data);
	}

	template <typename SomeType>
	RealType energy(const SomeType& y, SizeType threadNum) const
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

	void fromExpression(const std::string& str)
	{
		std::cerr << "Asumming Hamiltonian Expression (XACC) " << str << "\n";
		std::string paulis = toPaulis(removeAllSpaces(str));
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
		ToPauliMatrices toPauliMatrices(name);
		// name and bits <=== FIXME bits need adjustment
		return toPauliMatrices() + ttos(gate.bits());
	}

	SizeType bits_;
	std::string ham_;
	PauliOperatorType* pauliOperator_;
};

}
#endif // HAMILTONIAN_XACC_H
