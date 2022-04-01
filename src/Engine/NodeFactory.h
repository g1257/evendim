#ifndef NODEFACTORY_H
#define NODEFACTORY_H
#include "Vector.h"

namespace Gep {

/* Thread safe: Try the following
 *
 * if threadId > 0 copy the node and delete old angle/add new angle
 *
 * mark node
 *
 * return node
 *
 *
 * PROBLEMS WITH ABOVE APPROACH
 *
 * Return type is reference
 *
 * Too many copies?
 */
template<typename NodeType>
class NodeFactory {

public:

	typedef typename PsimagLite::Vector<NodeType*>::Type VectorNodeType;

	NodeFactory(const VectorNodeType& nodes) : nodes_(nodes)
	{}

	const NodeType& findNodeFromCode(PsimagLite::String codeStr,
	                                 const typename NodeType::ValueType& value,
	                                 bool isCell) const
	{
		PsimagLite::String codeStripped = stripPreviousAngleIfAny(codeStr);

		for (SizeType i = 0; i < nodes_.size(); ++i) {
			if (isCell && nodes_[i]->isInput()) continue;
			PsimagLite::String ncode = stripPreviousAngleIfAny(nodes_[i]->code());
			if (ncode == codeStripped) {
				if (codeStr == "?") nodes_[i]->setDcValue(value);
				nodes_[i]->setAngle(codeStr);
				return *nodes_[i];
			}
		}

		throw PsimagLite::RuntimeError("findNodeWithCode\n");
	}

	static PsimagLite::String stripPreviousAngleIfAny(PsimagLite::String str)
	{
		typename PsimagLite::String::const_iterator it = std::find(str.begin(),
		                                                           str.end(),
		                                                           ':');
		if (it == str.end()) return str; // no angle found

		return str.substr(0, it - str.begin());
	}

private:

	const VectorNodeType& nodes_;
};
}
#endif // NODEFACTORY_H
