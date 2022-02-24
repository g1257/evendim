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
#include "CanonicalFormQuantum.h"
#include "InputGatesUtil.h"

namespace Gep {

template<typename ValueType_>
class QuantumCircuit {

public:

	typedef QuantumCircuit<ValueType_> ThisType;
	typedef typename PsimagLite::Vector<ValueType_>::Type VectorValueType;
	typedef typename ValueType_::value_type ComplexType;
	typedef typename PsimagLite::Real<ComplexType>::Type RealType;
	typedef typename PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef Node<VectorValueType, RealType> NodeType;
	typedef typename PsimagLite::Vector<NodeType*>::Type VectorNodeType;
	typedef typename PsimagLite::Vector<VectorNodeType>::Type VectorVectorNodeType;
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
	typedef CanonicalFormQuantum<ValueType_, RealType> CanonicalFormType;
	typedef InputGatesUtil<ThisType> InputGatesUtilType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	QuantumCircuit(SizeType numberOfBits,
	               PsimagLite::String gates,
	               SizeType numberOfThreads)
	    : maxArity_(0), numberOfBits_(numberOfBits), rng_(1000)
	{
		PsimagLite::split(gates_, gates, ",");

		VectorNodeType nodes;
		makeNodes(nodes);

		for (SizeType i=0;i<nodes.size();i++) {
			if (nodes[i]->isInput()) {
				inputs_.push_back(i);
				terminals_.push_back(nodes[i]->code());
			} else if (nodes[i]->arity()>0 && nodes[i]->code()[0] != '_') {
				nonTerminals_.push_back(nodes[i]->code());
			}
		}

		for (SizeType i=0;i<nodes.size();i++) {
			if (maxArity_ < nodes[i]->arity())
				maxArity_ = nodes[i]->arity();
		}

		nodes_.push_back(nodes);

		for (SizeType i = 1; i < numberOfThreads; ++i) {
			VectorNodeType nodes;
			makeNodes(nodes);
			nodes_.push_back(nodes);
		}
	}

	~QuantumCircuit()
	{
		for (SizeType i = 0; i < nodes_.size(); i++) {
			for (SizeType j = 0; j < nodes_[i].size(); ++j) {
				delete nodes_[i][j];
				nodes_[i][j] = nullptr;
			}

			nodes_[i].clear();
		}

		nodes_.clear();
	}

	const VectorNodeType& nodes(SizeType threadNum) const
	{
		assert(threadNum < nodes_.size());
		return nodes_[threadNum];
	}

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

	SizeType numberOfInputs() const { return inputs_.size(); }

	void setInput(SizeType ind, ValueType x, SizeType threadId)
	{
		assert(ind < inputs_.size());
		assert(threadId < nodes_.size());
		assert(inputs_[ind] < nodes_[threadId].size());
		return nodes_[threadId][inputs_[ind]]->set(x);
	}

	void printInputs(std::ostream& os) const
	{
		assert(inputs_.size() > 0);
		assert(nodes_.size() > 0);
		assert(inputs_[inputs_.size() - 1] < nodes_[0].size());

		os<<"inputs= ";
		for (SizeType i = 0; i < inputs_.size(); i++)
			nodes_[0][inputs_[i]]->print(os);
		os<<"\n";
	}

	void sync()
	{

	}

private:

	void makeNodes(VectorNodeType& nodes)
	{
		static const SizeType inputs = 1;

		VectorStringType tmpGates = gates_;
		typename VectorStringType::const_iterator it = std::find(tmpGates.begin(),
		                                                         tmpGates.end(),
		                                                         "H");
		if (it != tmpGates.end()) {
			tmpGates.erase(it);
			// add Hadamard gates
			MatrixType hadamardGate;
			OneBitGateLibraryType::fillHadamard(hadamardGate);
			for (SizeType i = 0; i < numberOfBits_; ++i) {
				NodeType* hadamard = new QuantumOneBitGateType("H", i, numberOfBits_, hadamardGate);
				nodes.push_back(hadamard);
			}
		}

		it = std::find(tmpGates.begin(), tmpGates.end(), "P");
		if (it != tmpGates.end()) {
			tmpGates.erase(it);
			// add PHASE gates
			MatrixType phaseGate;
			OneBitGateLibraryType::fillPhase(phaseGate);
			for (SizeType i = 0; i < numberOfBits_; ++i) {
				NodeType* phase = new QuantumOneBitGateType("P", i, numberOfBits_, phaseGate);
				nodes.push_back(phase);
			}
		}

		// add Pauli matrices
		for (SizeType dir = 0; dir < 3; ++dir) {
			char charDir = OneBitGateLibraryType::directionIntegerToChar(dir);
			PsimagLite::String rDir("S ");
			rDir[1] = charDir;
			it = std::find(tmpGates.begin(), tmpGates.end(), rDir);
			if (it == tmpGates.end())
				continue;

			tmpGates.erase(it);
			MatrixType pauli;
			OneBitGateLibraryType::fillPauli(pauli, dir);
			for (SizeType i = 0; i < numberOfBits_; ++i) {
				NodeType* pauliNode = new QuantumOneBitGateType(rDir, i, numberOfBits_, pauli);
				nodes.push_back(pauliNode);
			}
		}

		// add rotation gates
		for (SizeType dir = 0; dir < 3; ++dir) {
			char charDir = OneBitGateLibraryType::directionIntegerToChar(dir);
			PsimagLite::String rDir("R ");
			rDir[1] = charDir;
			it = std::find(tmpGates.begin(), tmpGates.end(), rDir);
			if (it == tmpGates.end())
				continue;

			tmpGates.erase(it);
			MatrixType rotation;
			OneBitGateLibraryType::rotation(rotation, dir, 0); // 0 == angle

			for (SizeType i = 0; i < numberOfBits_; ++i) {
				NodeType* rot = new QuantumOneBitGateType(rDir, i, numberOfBits_, rotation);
				nodes.push_back(rot);
			}

			OneBitGateLibraryType::diffRotation(rotation, dir, 0); // 0 == angle
			rDir = "_R ";
			rDir[2] = charDir;
			for (SizeType i = 0; i < numberOfBits_; ++i) {
				NodeType* drot = new QuantumOneBitGateType(rDir, i, numberOfBits_, rotation);
				nodes.push_back(drot);
			}
		}

		it = std::find(tmpGates.begin(), tmpGates.end(), "C");
		if (it != tmpGates.end()) {
			tmpGates.erase(it);
			// add CNOT gates
			MatrixType cnotGate;
			TwoBitGateLibraryType::fillCnot(cnotGate);
			for (SizeType i = 0; i < numberOfBits_; ++i) {
				for (SizeType j = i + 1; j < numberOfBits_; ++j) {
					NodeType* cnot = new QuantumTwoBitGateType("C", i, j, numberOfBits_, cnotGate);
					nodes.push_back(cnot);
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
			NodeType* input = new QuantumInput<VectorValueType>(numberOfBits_);
			nodes.push_back(input);
		}
	}

	SizeType maxArity_;
	VectorValueType dcValues_;
	VectorStringType dcArray_;
	const SizeType numberOfBits_;
	VectorVectorNodeType nodes_;
	VectorStringType nonTerminals_;
	VectorStringType terminals_;
	VectorStringType gates_;
	VectorSizeType inputs_;
	mutable PsimagLite::MersenneTwister rng_; //RandomForTests<double> rng_;
}; // class QuantumCircuit

} // namespace Gep

#endif // EVENDIM_QUANTUM_CIRCUIT_H
