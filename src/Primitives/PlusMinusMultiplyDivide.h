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
#include "CanonicalFormEmpty.h"

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
	typedef CanonicalFormEmpty CanonicalFormType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	PlusMinusMultiplyDivide(SizeType inputs,
	                        SizeType genes,
	                        SizeType constants)
	    : dcValues_(constants), dcArray_(constants), rng_(1000)
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
			inputs_.push_back(nodes_.size());
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
	}

	~PlusMinusMultiplyDivide()
	{
		for (SizeType i = 0; i < nodes_.size(); i++)
			delete nodes_[i];

		nodes_.clear();
	}

	const VectorNodeType& nodes() const
	{
		return nodes_;
	}

	const VectorStringType& nonTerminals() const
	{
		return nonTerminals_;
	}

	const VectorStringType& terminals() const
	{
		return terminals_;
	}

	const VectorValueType& dcValues() const { return dcValues_; }

	const VectorStringType& dcArray() const { return dcArray_; }

	double rng() const { return rng_(); }

	SizeType numberOfInputs() const { return inputs_.size(); }

	void setInput(SizeType ind, SizeType x)
	{
		assert(ind < inputs_.size());
		assert(inputs_[ind] < nodes_.size());
		return nodes_[inputs_[ind]]->set(x);
	}

	void setInput(const VectorValueType& x) const
	{
		assert(x.size() == inputs_.size());
		SizeType n = std::min(x.size(), inputs_.size());
		assert(n > 0);
		assert(n < nodes_.size() + 1);
		for (SizeType i = 0; i < n; ++i) {
			nodes_[inputs_[i]]->set(x[i]);
		}
	}

	void printInputs(std::ostream& os) const
	{
		assert(inputs_.size() > 0);
		assert(inputs_[inputs_.size() - 1] < nodes_.size());

		os<<"inputs= ";
		for (SizeType i = 0; i < inputs_.size(); i++)
			nodes_[inputs_[i]]->print(os);
		os<<"\n";
	}

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

	VectorValueType dcValues_;
	VectorStringType dcArray_;
	VectorNodeType nodes_;
	VectorStringType nonTerminals_;
	VectorStringType terminals_;
	VectorSizeType inputs_;
	mutable PsimagLite::MersenneTwister rng_; //RandomForTests<double> rng_;
}; // class PlusMinusMultiplyDivide

} // namespace Gep

#endif // PLUS_MINUS_MULT_DIV_H
