/*
Copyright (c) 2017-2022, UT-Battelle, LLC

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
#ifndef NODEADF_H
#define NODEADF_H
#include "AST/Node.h"

namespace Gep {

template <typename VectorValueType, typename AnglesType_ = int>
class NodeAdf : public PsimagLite::Node<VectorValueType, AnglesType_> {

	typedef typename VectorValueType::value_type ValueType;

public:

	NodeAdf(SizeType i, ValueType input_)
	    : char_(i + 48)
	    , strOneChar_(" ")
	{
		strOneChar_[0] = char_;
	}

	NodeAdf* clone() const
	{
		return new NodeAdf(*this);
	}

	virtual PsimagLite::String code() const { return strOneChar_; }

	virtual SizeType arity() const { return 0; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 0);
		return value_;
	}

	virtual void set(const ValueType& x) const { value_ = x; }

private:

	char char_;
	PsimagLite::String strOneChar_;
	mutable ValueType value_;

}; // class NodeAdf
}
#endif // NODEADF_H
