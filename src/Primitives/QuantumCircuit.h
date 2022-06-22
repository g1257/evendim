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
#include "InputNg.h"
#include "InputCheck.h"

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
	typedef PsimagLite::InputNg<InputCheck>::Readable InputNgReadableType;

	QuantumCircuit(SizeType numberOfBits,
	               PsimagLite::String gates,
	               InputNgReadableType& io)
	    : numberOfBits_(numberOfBits), io_(io)
	{
		PsimagLite::split(gates_, gates, ",");

		makeNodes(nodes_);
	}

	~QuantumCircuit()
	{
		for (SizeType i = 0; i < nodes_.size(); i++) {
			delete nodes_[i];
			nodes_[i] = nullptr;
		}
	}

	const VectorNodeType& nodesSerial() const
	{
		return nodes_;
	}

	const VectorValueType& dcValues() const { return dcValues_; }

	const VectorStringType& dcArray() const { return dcArray_; }

	SizeType numberOfBits() const { return numberOfBits_; }

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
			OneBitGateLibraryType::fillAnyGate(hadamardGate, "H");
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
			OneBitGateLibraryType::fillAnyGate(phaseGate, "P");
			for (SizeType i = 0; i < numberOfBits_; ++i) {
				NodeType* phase = new QuantumOneBitGateType("P", i, numberOfBits_, phaseGate);
				nodes.push_back(phase);
			}
		}

		it = std::find(tmpGates.begin(), tmpGates.end(), "T");
		if (it != tmpGates.end()) {
			tmpGates.erase(it);
			// add T gates
			MatrixType tGate;
			OneBitGateLibraryType::fillAnyGate(tGate, "T");
			for (SizeType i = 0; i < numberOfBits_; ++i) {
				NodeType* t = new QuantumOneBitGateType("T", i, numberOfBits_, tGate);
				nodes.push_back(t);
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

		customGates(nodes, tmpGates);

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

	void customGates(VectorNodeType& nodes, VectorStringType& gates)
	{
		for (SizeType i = 0; i < gates.size(); ++i) {
			const std::string& name = gates[i];

			MatrixType matrix;
			io_.read(matrix, name);

			std::cout<<"Trying to add custom gate named " + name + "\n";

			if (matrix.rows() != matrix.cols())
				err("Matrix named " + name + " must be square.\n");
			if (matrix.rows() == 2) {
				fillCustomOneBitGates(nodes, name, matrix);
			} else if (matrix.rows() == 4) {
				fillCustomTwoBitGates(nodes, name, matrix);
			} else {
				err("Matrix named " + name + " must have either two or four rows.\n");
			}
		}

		gates.clear();
	}

	void fillCustomTwoBitGates(VectorNodeType& nodes, const PsimagLite::String& name, const MatrixType& matrix)
	{
		for (SizeType i = 0; i < numberOfBits_; ++i) {
			for (SizeType j = i + 1; j < numberOfBits_; ++j) {
				NodeType* customGate = new QuantumTwoBitGateType(name, i, j, numberOfBits_, matrix);
				nodes.push_back(customGate);
			}
		}
	}

	void fillCustomOneBitGates(VectorNodeType& nodes, const PsimagLite::String& name, const MatrixType& matrix)
	{
		for (SizeType i = 0; i < numberOfBits_; ++i) {
			NodeType* customGate = new QuantumOneBitGateType(name, i, numberOfBits_, matrix);
			nodes.push_back(customGate);
		}
	}

	const SizeType numberOfBits_;
	InputNgReadableType& io_;
	VectorValueType dcValues_;
	VectorStringType dcArray_;
	VectorNodeType nodes_;
	VectorStringType gates_;
}; // class QuantumCircuit

} // namespace Gep

#endif // EVENDIM_QUANTUM_CIRCUIT_H
