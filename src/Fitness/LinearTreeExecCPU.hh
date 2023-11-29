#ifndef LINEARTREEEXEC_CPU_HH
#define LINEARTREEEXEC_CPU_HH
#include "AST/Node.h"
#include "Hamiltonian.h"
#include "NodeFactory.h"
#include <complex>
#include <string>
#include <vector>

namespace Gep {

template <typename T>
struct UnderlyingType {
	using Type = T;
};

template <typename T>
struct UnderlyingType<QuasiVector<T>> {
	using Type = T;
};

template <typename ValueType, typename AnglesType_>
class LinearTreeExecCPU {

public:

	using VectorValueType = typename std::vector<ValueType>;
	using VectorStringType = std::vector<std::string>;
	using AnglesType = AnglesType_;
	using NodeType = PsimagLite::Node<VectorValueType, AnglesType>;
	using NodeFactoryType = NodeFactory<NodeType>;
	using ComplexType = typename UnderlyingType<ValueType>::Type;
	using RealType = typename PsimagLite::Real<ComplexType>::Type;
	using HamiltonianType = Hamiltonian<ComplexType>;
	using HandleType = std::pair<ValueType, SizeType>;

	explicit LinearTreeExecCPU(const NodeFactoryType& nodeFactory)
	    : nodeFactory_(nodeFactory)
	{
	}

	HandleType getHandle(const ValueType& initVector, // input has already been set
	                     const VectorStringType& circuit,
	                     SizeType threadNum) const
	{
		static const ValueType value;
		constexpr bool isCell = false;
		SizeType ngates = circuit.size();
		VectorValueType v(1, initVector);
		ValueType w;
		// here we could use commutation relations, order by site, etc TODO FIXME
		for (SizeType i = 0; i < ngates; ++i) {
			if (circuit[i] == "0")
				break;
			const NodeType& node = nodeFactory_.findNodeFromCode(circuit[i],
			                                                     value,
			                                                     isCell,
			                                                     threadNum);
			w = node.exec(v);
			v[0].swap(w);
		}

		return HandleType(v[0], threadNum);
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
