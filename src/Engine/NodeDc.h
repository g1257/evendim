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
#ifndef NODEDC_H
#define NODEDC_H
#include "AST/Node.h"

namespace Gep {

template <typename VectorValueType>
class NodeDc : public PsimagLite::Node<VectorValueType> {

	typedef typename VectorValueType::value_type ValueType;

public:

	NodeDc* clone() const
	{
		return new NodeDc(*this);
	}

	virtual PsimagLite::String code() const { return "?"; }

	virtual SizeType arity() const { return 0; }

	virtual void setDcValue(const ValueType& value) const
	{
		value_ = value;
	}

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 0);
		return value_;
	}

private:

	mutable ValueType value_;

}; // class NodeDc

}
#endif // NODEDC_H
