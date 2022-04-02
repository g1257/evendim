#ifndef NODEFACTORY_H
#define NODEFACTORY_H
#include "Vector.h"
#include "Concurrency.h"
#include "Matrix.h"

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
	typedef PsimagLite::Matrix<int> MatrixIntType;
	typedef PsimagLite::Vector<PsimagLite::Concurrency::PthreadtType>::Type VectorThreadIdType;

	NodeFactory(const VectorNodeType& nodes)
	    : nodes_(nodes),
	      nthreads_(PsimagLite::Concurrency::codeSectionParams.npthreads),
	      threadNum_(0),
	      newNodes_(nodes.size()*nthreads_),
	      table_(nthreads_, nodes.size())
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
				NodeType* newNode = findOrCreateCombo(i);
				if (codeStr == "?") newNode->setDcValue(value);
				newNode->setAngle(codeStr);
				return *newNode;
			}
		}

		throw PsimagLite::RuntimeError("findNodeWithCode\n");
	}

	void setThreadId(SizeType threadNum)
	{
		threadNum_ = threadNum;
	}

	void sync()
	{
		std::cerr<<"FINAL-->"<<nodes_.size()<<" vs. "<<newNodes_.size()<<"\n";
		clearNewNodes();
		//throw PsimagLite::RuntimeError("testing sync\n");
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

	NodeType* findOrCreateCombo(SizeType ind) const
	{
		if (threadNum_ == 0) return nodes_[ind];

		int tId = table_(threadNum_, ind);
		if (tId < 0) {
			auto newNode = nodes_[ind]->clone();
			SizeType loc = threadNum_ + ind*nthreads_;
			assert(loc < newNodes_.size());
			newNodes_[loc] = newNode;
			table_(threadNum_, ind) = loc;
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
	SizeType threadNum_;
	mutable VectorNodeType newNodes_;
	mutable MatrixIntType table_;
};
}
#endif // NODEFACTORY_H
