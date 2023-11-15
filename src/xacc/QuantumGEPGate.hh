#ifndef QUANTUMGEP_GATE_HH_H
#define QUANTUMGEP_GATE_HH_H
#include "AllocatorCpu.h"
#include <string>
#include <cassert>
#include "xacc.hpp"

namespace Gep
{

class QuantumGEPGate
{
public:

	using VectorSizeType = std::vector<SizeType>;
	using ParamType = xacc::Variant<int, double, std::string>;
	using VectorParamType = std::vector<ParamType>;
	using PairStringSizeType = std::pair<std::string, SizeType>;

	// Non-parametric non-custom gates can be for example
	// Sx10 
	// This is the Pauli X gate acting on bit ten (10)
	//
	// Another example:
	// C0_11
	// This is the CNOT gate acting on bits 0 and 11, the 
	// underscore separates the two bits
	//
	// Parametric may or maynot come with an "angle":
	// Ry71:3.1415927
	// is the Ry gate action on bit 71, with angle 3.1415927
	// But this gate is also valid:
	// Ry71
	// which has no angle
	//
	// Custom gates (both parametric and non-parametric)
	// will be ignored for now (TODO)
	//
	QuantumGEPGate(const std::string& s) : originalStr_(s), isParametric_(false)
	{
		setMapOfGates();
		// strip angle if any (we'll ignore the angle here for now)
		std::string str = stripPreviousAngleIfAny(originalStr_);

		// set bits (at least one must be present)
		SizeType counter = setBits(str);	
		
		// get the GEP name
		std::string gepName = getGEPName(counter, str);

		// convert GEP name to XACC name
		xaccName_ = gepToXaccName(gepName);

	}

	const std::string& name() { return xaccName_; }

	const VectorSizeType& bits() { return bits_; }

	bool isParametric() const { return isParametric_; }

	const VectorParamType& params() { return params_; }

private:

	void setMapOfGates()
	{
		gepToXaccGates_["Sx"] = PairStringSizeType("X", 1);
	}

	static std::string stripPreviousAngleIfAny(const std::string& str)
        {
                typename PsimagLite::String::const_iterator it = std::find(str.begin(),
                                                                           str.end(),
                                                                           ':');
                if (it == str.end()) return str; // no angle found

                return str.substr(0, it - str.begin());
        }

	static SizeType readNumberFromTheEnd(SizeType& counter, const std::string& str)
        {
                // find first non digit starting from the end 
                const SizeType l = str.length();
		std::string buffer;
                for (SizeType i = counter; i < l; ++i) {
                        SizeType j = l - i - 1; // we go in reverse
                        unsigned int num = str[j];
                        if (num > 57 || num < 48) // if it's not a digit
                                break;
                        ++counter;

			buffer += str[j];
                }

		if (buffer.size() == 0) {
			 throw std::runtime_error("Gate " + str + " has two few bits!?\n");
		}

		std::reverse(buffer.begin(), buffer.end());
		std::cout<<buffer<<"\n";
		return std::stoi(buffer);
	}

	SizeType setBits(const std::string& str)
	{
		SizeType counter = 0;
		SizeType bit1 = readNumberFromTheEnd(counter, str);
		bits_.push_back(bit1);

		// is there an underscore
		if (str[counter] == '_') {
			++counter; // step over underscore
			SizeType bit2 = readNumberFromTheEnd(counter, str);
			bits_.push_back(bit2);
		}

		return counter;
	}

	std::string getGEPName(SizeType counter, const std::string& str) const
	{
		const SizeType l = str.length();
		// what remains is the gep name
                SizeType ll = l - counter;
                assert(ll > 0);
		std::string gepName = str.substr(0, ll);
		assert(bits_.size() > 0);
		checkNameAndArity(gepName, bits_.size());
		return gepName;
        }

	void checkNameAndArity(const std::string& name, SizeType nbits) const
	{
		if (gepToXaccGates_.count(name) == 0) {
			throw std::runtime_error("Gate " + name + " not found\n");
		}

		const PairStringSizeType& gate = gepToXaccGates_.at(name);
		if (gate.second != nbits) {
			throw std::runtime_error("Gate " + name + " found with wrong number of bits\n");
		}
	}

	std::string gepToXaccName(const std::string& name) const
	{
		if (gepToXaccGates_.count(name) == 0) {
                        throw std::runtime_error("Gate " + name + " not found\n");
                }

		return gepToXaccGates_.at(name).first;
	}

	std::map<std::string, PairStringSizeType> gepToXaccGates_;
	std::string originalStr_;
	bool isParametric_;
	std::string xaccName_;
	VectorSizeType bits_;
	VectorParamType params_;
};
}

#endif
