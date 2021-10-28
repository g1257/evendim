#ifndef QUANTUMGATES_H
#define QUANTUMGATES_H
#include "Node.h"

namespace Gep {

template<typename VectorValueType>
class Hadamard : public Node<VectorValueType> {

	typedef typename VectorValueType::value_type ValueType;

public:

	Hadamard(SizeType bitNumber) : bitNumber_(bitNumber)
	{}

	virtual PsimagLite::String code() const { return "H"; }

	virtual SizeType arity() const { return 1; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 2);
		return v[0] + v[1];
	}

private:

	SizeType bitNumber_;

}; // class Plus

}

#endif // QUANTUMGATES_H
