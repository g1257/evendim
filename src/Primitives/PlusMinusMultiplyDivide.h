/*
Copyright (c) 2017, UT-Battelle, LLC

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
#ifndef PLUS_MINUS_MULT_DIV_H
#define PLUS_MINUS_MULT_DIV_H
#include "AST/Node.h"
#include "CanonicalFormEmpty.h"
#include "MersenneTwister.h"
#include "NodeAdf.h"
#include "NodeDc.h"
#include "Vector.h"
#include <cassert>

namespace Gep {

template <typename ValueType_>
class PlusMinusMultiplyDivide {

public:

	typedef typename PsimagLite::Vector<ValueType_>::Type VectorValueType;
	typedef PsimagLite::Node<VectorValueType> NodeType;
	typedef typename PsimagLite::Vector<NodeType*>::Type VectorNodeType;
	typedef NodeDc<VectorValueType> NodeDcType;
	typedef PsimagLite::Plus<VectorValueType> PlusType;
	typedef PsimagLite::Minus<VectorValueType> MinusType;
	typedef PsimagLite::Times<VectorValueType> TimesType;
	typedef PsimagLite::DividedBy<VectorValueType> DividedByType;
	typedef PsimagLite::Input<VectorValueType> InputType;
	typedef NodeAdf<VectorValueType> NodeAdfType;
	typedef ValueType_ ValueType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef CanonicalFormEmpty CanonicalFormType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	PlusMinusMultiplyDivide(SizeType inputs,
	                        SizeType genes,
	                        SizeType constants)
	    : dcValues_(constants)
	    , dcArray_(constants)
	    , rng_(1000)
	{
		addConstants();

		NodeType* plus = new PlusType();
		nodes_.push_back(plus);

		NodeType* minus = new MinusType();
		nodes_.push_back(minus);

		NodeType* times = new TimesType();
		nodes_.push_back(times);

		//		NodeType* dividedBy = new DividedByType();
		//		nodes_.push_back(dividedBy);

		for (SizeType i = 0; i < inputs; i++) {
			NodeType* input = new InputType(i);
			nodes_.push_back(input);
		}

		for (SizeType i = 0; i < genes; i++) {
			NodeType* adf = new NodeAdfType(i, 0);
			nodes_.push_back(adf);
		}
	}

	~PlusMinusMultiplyDivide()
	{
		for (SizeType i = 0; i < nodes_.size(); i++)
			delete nodes_[i];

		nodes_.clear();
	}

	const VectorValueType& dcValues() const { return dcValues_; }

	const VectorStringType& dcArray() const { return dcArray_; }

	const VectorNodeType& nodesSerial() const
	{
		return nodes_;
	}

private:

	void addConstants()
	{
		if (dcValues_.size() == 0)
			return;
		dcArray_.resize(dcValues().size());
		for (SizeType i = 0; i < dcValues_.size(); i++) {
			dcValues_[i] = 10.0 * rng_() - 10.0;
			dcArray_[i] = ttos(i);
		}

		NodeType* dc = new NodeDcType();
		nodes_.push_back(dc);
	}

	VectorValueType dcValues_;
	VectorStringType dcArray_;
	VectorNodeType nodes_;
	mutable PsimagLite::MersenneTwister rng_; // RandomForTests<double> rng_;
}; // class PlusMinusMultiplyDivide

} // namespace Gep

#endif // PLUS_MINUS_MULT_DIV_H
