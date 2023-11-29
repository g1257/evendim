#ifndef LINEARTREEEXEC_CPU_HH
#define LINEARTREEEXEC_CPU_HH
#include <vector>
#include <string>
#include <complex>
#include "AST/Node.h"
#include "NodeFactory.h"
#include "Hamiltonian.h"

namespace Gep {

template<typename RealType>
class LinearTreeExecCPU {

public:

	using VectorStringType = std::vector<std::string>;
	using ComplexType = std::vector<RealType>;
	using VectorComplexType = std::vector<ComplexType>;
	using VectorVectorComplexType = std::vector<VectorComplexType>;
	using AnglesType = RealType;
	using NodeType = PsimagLite::Node<VectorVectorComplexType, AnglesType>;
	using NodeFactorType = NodeFactory<NodeType>;
	using HamiltonianType = Hamiltonian<ComplexType>;
	using HandleType = std::pair<VectorComplexType, SizeType>;

	LinearTreeExecCPU(const NodeFactorType& nodeFactory)
	    : nodeFactory_(nodeFactory)
	{}

	HandleType getHandle(const VectorComplexType& initVector,
	                     const VectorStringType& circuit,
	                     SizeType threadNum)
	{
		static const VectorComplexType value;
		constexpr bool isCell = false;
		SizeType ngates = circuit.size();
		VectorComplexType v = initVector;
		VectorComplexType w;
		// here we could use commutation relations, order by site, etc TODO FIXME
		for (SizeType i = 0; i < ngates; ++i) {
			const NodeType& node = nodeFactory_.findNodeFromCode(circuit[i],
			                                                     value,
			                                                     isCell,
			                                                     threadNum);
			w = node.exec(v);
			v.swap(w);
		}

		return HandleType(v, threadNum);
	}

	RealType energy(const HandleType& handle,
	                const HamiltonianType& hamiltonian)
	{
		return hamiltonian.energy(handle.first, handle.second);
	}

private:

	const NodeFactorType& nodeFactory_;
};
}
#endif // LINEARTREEEXEC_CPU_HH
