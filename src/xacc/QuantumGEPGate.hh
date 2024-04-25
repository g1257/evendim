#ifndef QUANTUMGEP_GATE_HH_H
#define QUANTUMGEP_GATE_HH_H
#include "AllocatorCpu.h"
#include "xacc.hpp"
#include <cassert>
#include <string>

namespace Gep {

class QuantumGEPGate {
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
	QuantumGEPGate(const std::string& s)
	    : originalStr_(s)
	{
		if (s == "0") {
			xaccName_ = "Measure";
			bits_ = { 0 };
			return;
		}

		std::string measure = "Measure";
		unsigned int lmeasure = measure.length();
		if (s.substr(0, lmeasure) == measure) {
			xaccName_ = "Measure";
			SizeType counter = 0;
			unsigned int bit = readNumberFromTheEnd(counter, s);
			bits_ = { bit };
			std::cerr << "QuantumGEPGate Measure with bit " << bit << "\n";
			return;
		}

		if (gepToXaccGates_.size() == 0) {
			setMapOfGates();
		}

		// strip angle if any (we'll ignore the angle here for now)
		std::pair<std::string, double> strPair = stripPreviousAngleIfAny(originalStr_);

		// set bits (at least one must be present)
		SizeType counter = setBits(strPair.first);

		// get the GEP name
		std::string gepName = getGEPName(counter, strPair.first);

		// convert GEP name to XACC name
		xaccName_ = gepToXaccName(gepName);

		bool isParametric = isGateParametric(gepName);
		// set number of params, for now all parametric
		// gates have all one param
		if (isParametric) {
			angles_.resize(1, strPair.second);
		}
	}

	const std::string& name() { return xaccName_; }

	const VectorSizeType& bits() { return bits_; }

	SizeType numberOfParams() { return angles_.size(); }

	double param(SizeType ind) const
	{
		assert(ind < angles_.size());
		return angles_[ind];
	}

private:

	// Taken from XACC's xacc/quantum/gate/ir/CommonGates.hpp
	static void setMapOfGates()
	{
		gepToXaccGates_["Sx"] = PairStringSizeType("X", 1);
		gepToXaccGates_["Sy"] = PairStringSizeType("Y", 1);
		gepToXaccGates_["Sz"] = PairStringSizeType("Z", 1);

		// gates with the same names and one bit
		std::vector<std::string> sameNames { "Rx", "Ry", "Rz", "T", "H", "I" };
		for (std::vector<std::string>::const_iterator it = sameNames.begin(); it != sameNames.end(); ++it) {
			gepToXaccGates_[*it] = PairStringSizeType(*it, 1);
		}

		// iSwap 2
		// fSim 2
		// I 1

		// H 1
		gepToXaccGates_["H"] = PairStringSizeType("H", 1);

		// CNOT 2
		gepToXaccGates_["C"] = PairStringSizeType("CNOT", 2);

		// AnnealingInstruction variable
		// U1 1
		// Swap 2
		// U 1
		// Rphi 1
		// XX 2

		// Measure 0 ==> dealt with elsewhere

		// CZ 2
		// Cphase 2
		// XY 2
		// S 1
		// Sdg 1
		// T 1
		// Tdg 1
		// CY 2
		// CH 2
		// CRZ 2
		// Reset 1
		// RZZ 2
	}

	static bool isGateParametric(const std::string& gep_name)
	{
		return (gep_name == "Rx" || gep_name == "Ry" || gep_name == "Rz");
	}

	static std::pair<std::string, double> stripPreviousAngleIfAny(const std::string& str)
	{
		typename PsimagLite::String::const_iterator it = std::find(str.begin(),
		                                                           str.end(),
		                                                           ':');
		if (it == str.end()) {
			return std::pair<std::string, double>(str, 0.); // no angle found
		}

		std::string first = str.substr(0, it - str.begin());
		std::string angle = str.substr(it + 1 - str.begin(), str.end() - it - 1);
		return std::pair<std::string, double>(first, std::stod(angle));
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
		// std::cout << " str = "<< str <<" and buffer = "<< buffer << "\n";
		return std::stoi(buffer);
	}

	SizeType setBits(const std::string& str)
	{
		SizeType counter = 0;
		SizeType bit1 = readNumberFromTheEnd(counter, str);
		bits_.push_back(bit1);

		// is there an underscore
		assert(counter < str.length());
		SizeType location = str.length() - counter - 1;
		if (str[location] == '_') {
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

	static std::map<std::string, PairStringSizeType> gepToXaccGates_;
	std::string originalStr_;
	std::string xaccName_;
	VectorSizeType bits_;
	std::vector<double> angles_;
};
}

#endif
