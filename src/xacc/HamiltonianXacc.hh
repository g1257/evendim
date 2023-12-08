#ifndef HAMILTONIAN_XACC_H
#define HAMILTONIAN_XACC_H

#include "../Fitness/HamiltonianBase.hh"
#include "InputNg.h"

namespace Gep {

template <typename ComplexType>
class Hamiltonian {

public:

	using VectorType = std::vector<ComplexType>;
	using RealType = typename PsimagLite::Real<ComplexType>::Type;
	using PauliOperatorType = xacc::quantum::PauliOperator;

	Hamiltonian(typename InputNgType::Readable& io, SizeType /* numberOfThreads */ )
	    : pauliOperator_(nullptr)
	{
		io.readline(bits_, "NumberOfBits="); // == number of "sites"
		io.readline(ham, "Hamiltonian=");
		if (ham.substr(0, 5) == "file:") {
			unimplemented("Hamiltonian=file:\n");
		}

		if (ham == "IsingGraph" || ham == "zz") {
			throw std::runtime_error("Hamiltonian=IsingGraph or zz\n");
		}

		if (ham == "zxz" || ham == "xx") {
			throw std::runtime_error("Hamiltonian=zxz or xx\n");
		}

		fromExpression(str);;
	}

	HamiltonianXacc(const std::string& expression, SizeType /* numberOfThreads */ )
	: pauliOperator_(nullptr)
	{
		fromExpression(str);
	}

	RealType energy(const VectorType& y, SizeType threadNum) const
	{
		return 0;
	}

	SizeType numberOfSites() const
	{
		return bits_;
	}

private:

	HamiltonianXacc(const HamiltonianXacc&) = delete;

	HamiltonianXacc& operator=(const HamiltonianXacc&) = delete;

	void fromExpression(const std::string& str)
	{
		std::cerr << "Asumming Hamiltonian Expression (XACC)\n";
		pauliOperator_ = new PauliOperatorType(ham);
	}

	static void unimplemented(const std::string& msg)
	{
		throw std::runtime_error("XACC Backend: unimplemented: " + msg + "\n");
	}

	SizeType bits_;
	PauliOperatorType* pauliOperator_;
};

}
#endif // HAMILTONIAN_XACC_H
