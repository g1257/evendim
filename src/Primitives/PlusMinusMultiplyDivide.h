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
#include "Vector.h"
#include <cassert>
#include "Node.h"
#include "MersenneTwister.h"

namespace Gep {

template<typename ValueType_>
class PlusMinusMultiplyDivide {

public:

	typedef typename PsimagLite::Vector<ValueType_>::Type VectorValueType;
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
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	PlusMinusMultiplyDivide(SizeType inputs,
	                        SizeType genes,
	                        SizeType constants)
	    : maxArity_(0),dcValues_(constants),dcArray_(constants),rng_(1000)
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
			NodeType* input = new InputType(i,0);
			nodes_.push_back(input);
		}

		for (SizeType i = 0; i < genes; i++) {
			NodeType* adf = new NodeAdfType(i,0);
			nodes_.push_back(adf);
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

	~PlusMinusMultiplyDivide()
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

private:

	void addConstants()
	{
		if (dcValues_.size() == 0) return;
		dcArray_.resize(dcValues().size());
		for (SizeType i = 0; i < dcValues_.size(); i++) {
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
	VectorNodeType nodes_;
	VectorStringType nonTerminals_;
	VectorStringType terminals_;
}; // class PlusMinusMultiplyDivide

} // namespace Gep

#endif // PLUS_MINUS_MULT_DIV_H
