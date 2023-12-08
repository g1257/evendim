#ifndef HAMILTONIAN_XACC_H
#define HAMILTONIAN_XACC_H

#include "../Engine/InputCheck.h"
#include "InputNg.h"
#include "PauliOperator.hpp"
#include "xacc.hpp"

namespace Gep {

template <typename ComplexType>
class Hamiltonian {

public:

	using VectorType = std::vector<ComplexType>;
	using RealType = typename PsimagLite::Real<ComplexType>::Type;
	using PauliOperatorType = xacc::quantum::PauliOperator;
	using InputNgType = PsimagLite::InputNg<InputCheck>;

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
		std::cerr << "Asumming Hamiltonian Expression (XACC)\n";
		pauliOperator_ = new PauliOperatorType(str);
	}

	static void unimplemented(const std::string& msg)
	{
		throw std::runtime_error("XACC Backend: unimplemented: " + msg + "\n");
	}

	SizeType bits_;
	std::string ham_;
	PauliOperatorType* pauliOperator_;
};

}
#endif // HAMILTONIAN_XACC_H
