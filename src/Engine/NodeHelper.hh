#ifndef NODEHELPER_HH
#define NODEHELPER_HH
#include "../Fitness/LinearTreeExec.hh"
#include "NodeFactory.h"
#include "UnderlyingType.hh"

namespace Gep {

template <typename NodeType>
class NodeHelper {

public:

	using ValueType = typename NodeType::ValueType;
	using AnglesType = typename NodeType::AnglesType;
	using NodeFactoryType = NodeFactory<NodeType>;
	using VectorNodeType = std::vector<NodeType*>;
	using VectorSizeType = std::vector<SizeType>;
	using VectorStringType = std::vector<std::string>;
	using VectorValueType = std::vector<ValueType>;
	using ComplexType = typename UnderlyingType<ValueType>::Type;
	using RealType = typename PsimagLite::Real<ComplexType>::Type;
	using LinearTreeExecType = LinearTreeExec<ValueType, AnglesType, NodeFactoryType>;

	NodeHelper(const VectorNodeType& nodes)
	    : nodeFactory_(nodes)
	    , linearTreeExec_(nodeFactory_)
	{
	}

	const LinearTreeExecType& linearTreeExec() const { return linearTreeExec_; }

	void setInput(SizeType ind, ValueType x, SizeType threadNum)
	{
		assert(ind < inputs_.size());
		assert(inputs_[ind] < nodeFactory_.numberOfNodes());
		return nodeFactory_.node(inputs_[ind], threadNum).set(x);
	}

	SizeType numberOfInputs() const { return inputs_.size(); }

	void setInput(const VectorValueType& x) const
	{
		assert(x.size() == inputs_.size());
		SizeType n = std::min(x.size(), inputs_.size());
		assert(n > 0);
		assert(n < nodeFactory_.numberOfNodes() + 1);
		const SizeType threadNum = 0;
		for (SizeType i = 0; i < n; ++i) {
			nodeFactory_.node(inputs_[i], threadNum).set(x[i]);
		}
	}

	void printInputs(std::ostream& os) const
	{
		assert(nodeFactory_.numberOfNodes() > 0);

		const SizeType threadNum = 0;
		os << "inputs= ";
		for (SizeType i = 0; i < inputs_.size(); i++) {
			SizeType j = inputs_[i];
			nodeFactory_.node(j, threadNum).print(os);
		}

		os << "\n";
	}

	SizeType maxArity() const
	{
		SizeType threadNum = 0;
		SizeType maxArity = 0;
		for (SizeType i = 0; i < nodeFactory_.numberOfNodes(); ++i) {
			if (maxArity < nodeFactory_.node(i, threadNum).arity())
				maxArity = nodeFactory_.node(i, threadNum).arity();
		}

		return maxArity;
	}

	void setInputsTerminalsAndNonTerminals(VectorStringType& terminals,
	                                       VectorStringType& nonTerminals)
	{
		SizeType threadNum = 0;
		for (SizeType i = 0; i < nodeFactory_.numberOfNodes(); ++i) {
			if (nodeFactory_.node(i, threadNum).isInput()) {
				inputs_.push_back(i);
				terminals.push_back(nodeFactory_.node(i, threadNum).code());
			}
			else if (nodeFactory_.node(i, threadNum).arity() > 0 && nodeFactory_.node(i, threadNum).code()[0] != '_') {
				nonTerminals.push_back(nodeFactory_.node(i, threadNum).code());
			}
		}
	}

	const NodeFactoryType& nodeFactory() const { return nodeFactory_; }

	// Try to remove?
	NodeFactoryType& nodeFactory() { return nodeFactory_; }

private:

	NodeFactoryType nodeFactory_;
	VectorSizeType inputs_;
	LinearTreeExecType linearTreeExec_;
};

}
#endif // NODEHELPER_HH
