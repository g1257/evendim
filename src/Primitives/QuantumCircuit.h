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
#include "Vector.h"
#include <cassert>
#include "QuantumOneBitGate.h"
#include "QuantumTwoBitGate.h"
#include "MersenneTwister.h"
#include "QuantumInput.h"

namespace Gep {

template<typename ValueType_>
class QuantumCircuit {

public:

	typedef typename PsimagLite::Vector<ValueType_>::Type VectorValueType;
	typedef typename PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef Node<VectorValueType> NodeType;
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

	QuantumCircuit(SizeType inputs,
	               SizeType genes,
	               SizeType numberOfBits)
	    : maxArity_(0), rng_(1000), numberOfBits_(numberOfBits)
	{

		// add Hadamard gates
		MatrixType hadamardGate;
		OneBitGateLibraryType::fillHadamard(hadamardGate);
		for (SizeType i = 0; i < numberOfBits; ++i) {
			NodeType* hadamard = new QuantumOneBitGateType('H', i, numberOfBits, hadamardGate);
			nodes_.push_back(hadamard);
		}

		// add PHASE gates
		MatrixType phaseGate;
		OneBitGateLibraryType::fillPhase(phaseGate);
		for (SizeType i = 0; i < numberOfBits; ++i) {
			NodeType* phase = new QuantumOneBitGateType('P', i, numberOfBits, phaseGate);
			nodes_.push_back(phase);
		}

		// add CNOT gates
		MatrixType cnotGate;
		TwoBitGateLibraryType::fillCnot(cnotGate);
		for (SizeType i = 0; i < numberOfBits; ++i) {
			for (SizeType j = i + 1; j < numberOfBits; ++j) {
				NodeType* cnot = new QuantumTwoBitGateType('C', i, j, numberOfBits, cnotGate);
				nodes_.push_back(cnot);
			}
		}

		for (SizeType i = 0; i < inputs; i++) {
			NodeType* input = new QuantumInput<VectorValueType>(numberOfBits);
			nodes_.push_back(input);
		}

		for (SizeType i=0;i<nodes_.size();i++) {
			if (nodes_[i]->isInput()) {
				terminals_.push_back(nodes_[i]->code());
			} else if (nodes_[i]->arity()>0) {
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
}; // class QuantumCircuit

} // namespace Gep

#endif // EVENDIM_QUANTUM_CIRCUIT_H
