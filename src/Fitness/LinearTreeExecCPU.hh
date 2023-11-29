#ifndef LINEARTREEEXEC_CPU_HH
#define LINEARTREEEXEC_CPU_HH
#include <vector>
#include <string>
#include <complex>
#include "AST/Node.h"
#include "NodeFactory.h"
#include "Hamiltonian.h"

namespace Gep {

template<typename T>
struct UnderlyingType {
	using Type = T;
};

template<typename T>
struct UnderlyingType<QuasiVector<T> > {
	using Type = T;
};

template<typename NodeType>
class LinearTreeExecCPU {

public:

	using VectorStringType = std::vector<std::string>;
	using AnglesType = typename NodeType::AnglesType;
	using NodeFactoryType = NodeFactory<NodeType>;
	using ValueType = typename NodeType::ValueType;
	using ComplexType = typename UnderlyingType<ValueType>::Type;
	using RealType = typename PsimagLite::Real<ComplexType>::Type;
	using HamiltonianType = Hamiltonian<ComplexType>;
	using HandleType = std::pair<ValueType, SizeType>;

	explicit LinearTreeExecCPU(const NodeFactoryType& nodeFactory)
	    : nodeFactory_(nodeFactory)
	{}

	HandleType getHandle(const ValueType& initVector,
	                     const VectorStringType& circuit,
	                     SizeType threadNum) const
	{
		static const ValueType value;
		constexpr bool isCell = false;
		SizeType ngates = circuit.size();
		ValueType v = initVector;
		ValueType w;
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
	                const HamiltonianType& hamiltonian) const
	{
		return hamiltonian.energy(handle.first, handle.second);
	}

private:

	LinearTreeExecCPU(const LinearTreeExecCPU&) = delete;

	LinearTreeExecCPU& operator=(const LinearTreeExecCPU&) = delete;

	const NodeFactoryType& nodeFactory_;
};
}
#endif // LINEARTREEEXEC_CPU_HH
