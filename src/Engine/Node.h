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
#ifndef NODE_H
#define NODE_H
#include "Vector.h"

namespace Gep {

template<typename VectorValueType, typename AnglesType_ = int>
class Node {

public:

	typedef AnglesType_ AnglesType;
	typedef typename VectorValueType::value_type ValueType;
	typedef typename PsimagLite::Vector<AnglesType>::Type VectorAnglesType;

	virtual ~Node() {}

	virtual PsimagLite::String code() const = 0;

	virtual SizeType arity() const = 0;

	virtual ValueType exec(const VectorValueType&) const = 0;

	virtual ValueType exec(const VectorValueType&,
	                       const VectorAnglesType*,
	                       SizeType&) const
	{
		throw PsimagLite::RuntimeError("node::exec() long form\n");
	}

	virtual void set(const ValueType&) const
	{
		throw PsimagLite::RuntimeError("node::set\n");
	}

	virtual void print(std::ostream&) {}

	virtual void setDcValue(const ValueType&) const
	{
		throw PsimagLite::RuntimeError("node::setDcValue\n");
	}

	virtual bool isInput() const  { return false; }

}; // class Node

template<typename VectorValueType>
class NodeDc : public Node<VectorValueType> {

	typedef typename VectorValueType::value_type ValueType;

public:

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

template<typename VectorValueType>
class Plus : public Node<VectorValueType> {

	typedef typename VectorValueType::value_type ValueType;

public:

	virtual PsimagLite::String code() const { return "+"; }

	virtual SizeType arity() const { return 2; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 2);
		return v[0] + v[1];
	}

}; // class Plus

template<typename VectorValueType>
class Minus : public Node<VectorValueType> {

	typedef typename VectorValueType::value_type ValueType;

public:

	virtual PsimagLite::String code() const { return "-"; }

	virtual SizeType arity() const { return 2; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 2);
		return v[0] - v[1];
	}

}; // class Minus

template<typename VectorValueType>
class Times : public Node<VectorValueType> {

	typedef typename VectorValueType::value_type ValueType;

public:

	virtual PsimagLite::String code() const { return "*"; }

	virtual SizeType arity() const { return 2; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 2);
		return v[0] * v[1];
	}

}; // class Times

template<typename VectorValueType>
class DividedBy : public Node<VectorValueType> {

	typedef typename VectorValueType::value_type ValueType;

public:

	virtual PsimagLite::String code() const { return "/"; }

	virtual SizeType arity() const { return 2; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 2);
		if (fabs(v[1]) < 1e-6) return v[0];

		return v[0] / v[1];
	}

}; // class DividedBy

template<typename VectorValueType>
class IfGtZero : public Node<VectorValueType> {

	typedef typename VectorValueType::value_type ValueType;

public:

	virtual char code() const { return 'g'; }

	virtual SizeType arity() const { return 1; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 1);

		return (v[0] > 0) ? 1 : 0;
	}

}; // class IfGtZero

template<typename VectorValueType>
class Int : public Node<VectorValueType> {

	typedef typename VectorValueType::value_type ValueType;

public:

	virtual char code() const { return 'i'; }

	virtual SizeType arity() const { return 1; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 1);

		SizeType x = static_cast<SizeType>(v[0]);
		return x;
	}

}; // class Int

template<typename VectorValueType>
class Input : public Node<VectorValueType> {

	typedef typename VectorValueType::value_type ValueType;

public:

	Input(SizeType i, ValueType)
	    : char_(i+48),strOneChar_(" ")
	{
		strOneChar_[0] = char_;
	}

	virtual PsimagLite::String code() const { return strOneChar_; }

	virtual SizeType arity() const { return 0; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 0);
		return input_;
	}

	virtual void set(const ValueType& x) const { input_ = x; }

	virtual bool isInput() const  { return true; }

private:

	char char_;
	PsimagLite::String strOneChar_;
	mutable ValueType input_;

}; // class Input

template<typename VectorValueType>
class NodeAdf : public Node<VectorValueType> {

	typedef typename VectorValueType::value_type ValueType;

public:

	NodeAdf(SizeType i,ValueType input_)
	    : char_(i+48),strOneChar_(" ")
	{
		strOneChar_[0] = char_;
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

} // namespace Gep
#endif // NODE_H
