#ifndef TOPAULIMATRICES_HH
#define TOPAULIMATRICES_HH
#include <stdexcept>
#include <string>

namespace Gep {

class ToPauliMatrices {
public:

	ToPauliMatrices(const std::string& name)
	    : name_(name)
	{
		expansion_ = expandIntoPaulis(name);
	}

	const std::string& operator()() const { return expansion_; }

private:

	static std::string expandIntoPaulis(const std::string& name)
	{
		if (name == "X" || name == "Y" || name == "Z") {
			return name;
		}

		throw std::runtime_error("I don't know how to exand " + name + " into Pauli matrices\n");
	}

	std::string name_;
	std::string expansion_;
};
}
#endif // TOPAULIMATRICES_HH
