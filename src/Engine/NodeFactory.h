#ifndef NODEFACTORY_H
#define NODEFACTORY_H
#include "Concurrency.h"
#include "Matrix.h"
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
template <typename NodeType>
class NodeFactory {

public:

	typedef typename PsimagLite::Vector<NodeType*>::Type VectorNodeType;
	typedef PsimagLite::Matrix<int> MatrixIntType;
	typedef PsimagLite::Vector<PsimagLite::Concurrency::PthreadtType>::Type VectorThreadIdType;

	NodeFactory(const VectorNodeType& nodes)
	    : nodes_(nodes)
	    , nthreads_(PsimagLite::Concurrency::codeSectionParams.npthreads)
	    , newNodes_(nodes.size() * nthreads_)
	{
	}

	~NodeFactory()
	{
		clearNewNodes();
	}

	const NodeType& findNodeFromCode(PsimagLite::String codeStr,
	                                 const typename NodeType::ValueType& value,
	                                 bool isCell,
	                                 SizeType threadNum) const
	{
		PsimagLite::String codeStripped = stripPreviousAngleIfAny(codeStr);

		for (SizeType i = 0; i < nodes_.size(); ++i) {
			if (isCell && nodes_[i]->isInput())
				continue;
			PsimagLite::String ncode = stripPreviousAngleIfAny(nodes_[i]->code());
			if (ncode == codeStripped) {
				NodeType* newNode = findOrCreateCombo(i, threadNum);
				if (codeStr == "?")
					newNode->setDcValue(value);
				newNode->setAngle(codeStr);
				return *newNode;
			}
		}

		throw PsimagLite::RuntimeError("findNodeWithCode\n");
	}

	void sync()
	{
		// std::cerr << "FINAL-->" << nodes_.size() << " vs. " << newNodes_.size() << "\n";
		clearNewNodes();
	}

	static PsimagLite::String stripPreviousAngleIfAny(PsimagLite::String str)
	{
		typename PsimagLite::String::const_iterator it = std::find(str.begin(),
		                                                           str.end(),
		                                                           ':');
		if (it == str.end())
			return str; // no angle found

		return str.substr(0, it - str.begin());
	}

	NodeType& node(SizeType ind, SizeType threadNum)
	{
		NodeType* newNode = findOrCreateCombo(ind, threadNum);
		return *newNode;
	}

	const NodeType& node(SizeType ind, SizeType threadNum) const
	{
		NodeType* newNode = findOrCreateCombo(ind, threadNum);
		return *newNode;
	}

	SizeType numberOfNodes() const { return nodes_.size(); }

private:

	NodeType* findOrCreateCombo(SizeType ind, SizeType threadNum) const
	{
		int tId = threadNum + ind * nthreads_;
		assert(static_cast<SizeType>(tId) < newNodes_.size());
		if (!newNodes_[tId]) {
			newNodes_[tId] = nodes_[ind]->clone();
		}

		assert(static_cast<SizeType>(tId) < newNodes_.size());
		assert(newNodes_[tId]);
		return newNodes_[tId];
	}

	void clearNewNodes()
	{
		SizeType n = newNodes_.size();
		for (SizeType i = 0; i < n; ++i) {
			delete newNodes_[i];
			newNodes_[i] = nullptr;
		}
	}

	const VectorNodeType& nodes_;
	SizeType nthreads_;
	mutable VectorNodeType newNodes_;
};
}
#endif // NODEFACTORY_H
