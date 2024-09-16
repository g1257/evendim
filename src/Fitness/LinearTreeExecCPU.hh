#ifndef LINEARTREEEXEC_CPU_HH
#define LINEARTREEEXEC_CPU_HH
#include "AST/Node.h"
#include "HamiltonianCPU.hh"
#include "NodeFactory.h"
#include "UnderlyingType.hh"
#include <complex>
#include <string>
#include <vector>

template <typename T1, typename T2>
struct TypesMustBeEqual { };

template <typename T>
struct TypesMustBeEqual<T, T> {
	using Type = int;
};

namespace Gep {

template <typename ValueType, typename AnglesType_, typename NodeFactoryType_>
class LinearTreeExec {

public:

	static constexpr bool HAS_XACC = false;

	using VectorValueType = typename std::vector<ValueType>;
	using VectorStringType = std::vector<std::string>;
	using AnglesType = AnglesType_;
	using NodeType = PsimagLite::Node<VectorValueType, AnglesType>;
	using NodeFactoryType = NodeFactory<NodeType>;
	using ComplexType = typename UnderlyingType<ValueType>::Type;
	using RealType = typename PsimagLite::Real<ComplexType>::Type;
	using HamiltonianType = Hamiltonian<ComplexType>;
	using HandleType = std::pair<ValueType, SizeType>;

	explicit LinearTreeExec(const NodeFactoryType& nodeFactory,
	                        const typename TypesMustBeEqual<NodeFactoryType_, NodeFactoryType>::Type = 0)
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
	                const HamiltonianType& hamiltonian,
	                bool useXaccOptimizer) const
	{
		assert(!useXaccOptimizer);
		return hamiltonian.energy(handle.first, handle.second);
	}

	void fillAngles(std::vector<double>&, const HandleType&) const
	{
		err("fillAngles should not be called unless XACC is used\n");
	}

private:

	LinearTreeExec(const LinearTreeExec&) = delete;

	LinearTreeExec& operator=(const LinearTreeExec&) = delete;

	const NodeFactoryType& nodeFactory_;
};
}
#endif // LINEARTREEEXEC_CPU_HH
