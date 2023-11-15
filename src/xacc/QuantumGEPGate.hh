#ifndef QUANTUMGEP_GATE_HH_H
#define QUANTUMGEP_GATE_HH_H
#include "AllocatorCpu.h"
#include <string>
#include "xacc.hpp"

namespace Gep
{

class QuantumGEPGate
{
public:

	using VectorSizeType = std::vector<SizeType>;
	using ParamType = xacc::Variant<int, double, std::string>;
	using VectorParamType = std::vector<ParamType>;

	QuantumGEPGate(const std::string& s) : originalStr_(s), isParametric_(false)
	{
		// strip angle
		std::string str = stripPreviousAngleIfAny(originalStr_);
		throw std::runtime_error("Not fully implemented yet (sorry)\n");
	}

	const std::string& name() { return xaccName_; }

	const VectorSizeType& bits() { return bits_; }

	bool isParametric() const { return isParametric_; }

	const VectorParamType& params() { return params_; }

private:

	static std::string stripPreviousAngleIfAny(const std::string& str)
        {
                typename PsimagLite::String::const_iterator it = std::find(str.begin(),
                                                                           str.end(),
                                                                           ':');
                if (it == str.end()) return str; // no angle found

                return str.substr(0, it - str.begin());
        }

	std::string originalStr_;
	bool isParametric_;
	std::string xaccName_;
	VectorSizeType bits_;
	VectorParamType params_;
};
}

#endif
