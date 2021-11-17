/*
Copyright (c) 2017-2021, UT-Battelle, LLC

evendim, Version 0.

This file is part of evendim.
evendim is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
evendim is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with evendim. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef EVENDIM_QUANTUM_CIRCUIT_H
#define EVENDIM_QUANTUM_CIRCUIT_H
#include "PsimagLite.h"
#include <cassert>
#include "QuantumOneBitGate.h"
#include "QuantumTwoBitGate.h"
#include "MersenneTwister.h"
#include "QuantumInput.h"
#include <numeric>

namespace Gep {

template<typename ValueType_>
class QuantumCircuit {

public:

	typedef typename PsimagLite::Vector<ValueType_>::Type VectorValueType;
	typedef typename ValueType_::value_type ComplexType;
	typedef typename PsimagLite::Real<ComplexType>::Type RealType;
	typedef typename PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef Node<VectorValueType, RealType> NodeType;
	typedef typename PsimagLite::Vector<NodeType*>::Type VectorNodeType;
	typedef NodeDc<VectorValueType> NodeDcType;
	typedef Plus<VectorValueType> PlusType;
	typedef Minus<VectorValueType> MinusType;
	typedef Times<VectorValueType> TimesType;
	typedef DividedBy<VectorValueType> DividedByType;
	typedef Input<VectorValueType> InputType;
	typedef NodeAdf<VectorValueType> NodeAdfType;
	typedef ValueType_ ValueType;
	typedef QuantumOneBitGate<VectorValueType> QuantumOneBitGateType;
	typedef QuantumTwoBitGate<VectorValueType> QuantumTwoBitGateType;
	typedef typename QuantumOneBitGateType::MatrixType MatrixType;
	typedef OneBitGateLibrary<typename ValueType::value_type> OneBitGateLibraryType;
	typedef TwoBitGateLibrary<typename ValueType::value_type> TwoBitGateLibraryType;

	static const bool hasAngles = true;

	QuantumCircuit(SizeType numberOfBits,
	               PsimagLite::String gates)
	    : maxArity_(0), rng_(1000), numberOfBits_(numberOfBits)
	{
		static const SizeType inputs = 1;

		PsimagLite::split(gates_, gates, ",");

		VectorStringType tmpGates = gates_;
		typename VectorStringType::const_iterator it = std::find(tmpGates.begin(),
		                                                         tmpGates.end(),
		                                                         "H");
		if (it != tmpGates.end()) {
			tmpGates.erase(it);
			// add Hadamard gates
			MatrixType hadamardGate;
			OneBitGateLibraryType::fillHadamard(hadamardGate);
			for (SizeType i = 0; i < numberOfBits; ++i) {
				NodeType* hadamard = new QuantumOneBitGateType("H", i, numberOfBits, hadamardGate);
				nodes_.push_back(hadamard);
			}
		}

		it = std::find(tmpGates.begin(), tmpGates.end(), "P");
		if (it != tmpGates.end()) {
			tmpGates.erase(it);
			// add PHASE gates
			MatrixType phaseGate;
			OneBitGateLibraryType::fillPhase(phaseGate);
			for (SizeType i = 0; i < numberOfBits; ++i) {
				NodeType* phase = new QuantumOneBitGateType("P", i, numberOfBits, phaseGate);
				nodes_.push_back(phase);
			}
		}

		it = std::find(tmpGates.begin(), tmpGates.end(), "R");
		if (it != tmpGates.end()) {
			tmpGates.erase(it);
			// add rotation gates
			for (SizeType dir = 0; dir < 3; ++dir) {
				MatrixType rotation;
				OneBitGateLibraryType::rotation(rotation, dir, 0); // 0 == angle
				PsimagLite::String rDir("R ");
				rDir[1] = OneBitGateLibraryType::directionIntegerToChar(dir);

				for (SizeType i = 0; i < numberOfBits; ++i) {
					NodeType* rot = new QuantumOneBitGateType(rDir, i, numberOfBits, rotation);
					nodes_.push_back(rot);
				}

				OneBitGateLibraryType::diffRotation(rotation, dir, 0); // 0 == angle
				rDir = "_R ";
				rDir[2] = OneBitGateLibraryType::directionIntegerToChar(dir);
				for (SizeType i = 0; i < numberOfBits; ++i) {
					NodeType* drot = new QuantumOneBitGateType(rDir, i, numberOfBits, rotation);
					nodes_.push_back(drot);
				}
			}
		}

		it = std::find(tmpGates.begin(), tmpGates.end(), "C");
		if (it != tmpGates.end()) {
			tmpGates.erase(it);
			// add CNOT gates
			MatrixType cnotGate;
			TwoBitGateLibraryType::fillCnot(cnotGate);
			for (SizeType i = 0; i < numberOfBits; ++i) {
				for (SizeType j = i + 1; j < numberOfBits; ++j) {
					NodeType* cnot = new QuantumTwoBitGateType("C", i, j, numberOfBits, cnotGate);
					nodes_.push_back(cnot);
				}
			}
		}

		if (tmpGates.size() > 0) {
			PsimagLite::String tmp = std::accumulate(tmpGates.begin(),
			                                         tmpGates.end(),
			                                         PsimagLite::String(" "));
			err("The following gates were not recognized: " + tmp + "\n");
		}

		for (SizeType i = 0; i < inputs; i++) {
			NodeType* input = new QuantumInput<VectorValueType>(numberOfBits);
			nodes_.push_back(input);
		}

		for (SizeType i=0;i<nodes_.size();i++) {
			if (nodes_[i]->isInput()) {
				terminals_.push_back(nodes_[i]->code());
			} else if (nodes_[i]->arity()>0 && nodes_[i]->code()[0] != '_') {
				nonTerminals_.push_back(nodes_[i]->code());
			}
		}

		for (SizeType i=0;i<nodes_.size();i++) {
			if (maxArity_ < nodes_[i]->arity())
				maxArity_ = nodes_[i]->arity();
		}
	}

	~QuantumCircuit()
	{
		for (SizeType i = 0; i < nodes_.size(); i++)
			delete nodes_[i];

		nodes_.clear();
	}

	const VectorNodeType& nodes() const { return nodes_; }

	const VectorStringType& nonTerminals() const
	{
		return nonTerminals_;
	}

	const VectorStringType& terminals() const
	{
		return terminals_;
	}

	SizeType arity() const { return maxArity_; }

	bool hasDc() const { return (dcValues_.size() > 0); }

	const VectorValueType& dcValues() const { return dcValues_; }

	const VectorStringType& dcArray() const { return dcArray_; }

	double rng() const { return rng_(); }

	SizeType numberOfBits() const { return numberOfBits_; }

private:

	void addConstants()
	{
		const SizeType n = dcValues_.size();
		if (n == 0) return;
		dcArray_.resize(n);
		for (SizeType i = 0; i < n; ++i) {
			dcValues_[i] = 10.0*rng_() - 10.0;
			dcArray_[i] = ttos(i);
		}

		NodeType* dc = new NodeDcType();
		nodes_.push_back(dc);
	}

	SizeType maxArity_;
	VectorValueType dcValues_;
	VectorStringType dcArray_;
	mutable PsimagLite::MersenneTwister rng_; //RandomForTests<double> rng_;
	const SizeType numberOfBits_;
	VectorNodeType nodes_;
	VectorStringType nonTerminals_;
	VectorStringType terminals_;
	VectorStringType gates_;
}; // class QuantumCircuit

} // namespace Gep

#endif // EVENDIM_QUANTUM_CIRCUIT_H
