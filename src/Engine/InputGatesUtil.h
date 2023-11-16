#ifndef INPUTGATESUTIL_H
#define INPUTGATESUTIL_H
#include "PsimagLite.h"

namespace Gep {

template <typename PrimitivesType>
class InputGatesUtil {

public:

	typedef typename PrimitivesType::VectorValueType VectorValueType;
	typedef typename VectorValueType::value_type ValueType;

	InputGatesUtil(const PrimitivesType& primitives,
	               SizeType numberOfThreads)
	    : primitives_(primitives)
	    , numberOfThreads_(numberOfThreads)
	    , inputs_(0)
	{
		const auto& nodes = primitives_.nodes(0);

		for (SizeType i = 0; i < nodes.size(); ++i) {
			if (nodes[i]->isInput()) {
				++inputs_;
			}
		}
	}

	void setInput(SizeType kk, ValueType x)
	{
		SizeType k = 0;
		for (SizeType j = 0; j < numberOfThreads_; ++j) {
			auto& nodes = primitives_.nodes(j);
			for (SizeType i = 0; i < nodes.size(); ++i) {
				if (!nodes[i]->isInput())
					continue;
				if (k != kk) {
					++k;
					continue;
				}

				nodes[i]->set(x);
				++k;
			}
		}
	}

	void setInput(const VectorValueType& x) const
	{
		for (SizeType j = 0; j < numberOfThreads_; ++j) {
			auto& nodes = primitives_.nodes(j);
			SizeType k = 0;
			for (SizeType i = 0; i < nodes.size(); ++i) {
				if (!nodes[i]->isInput())
					continue;

				assert(k < x.size());
				nodes[i]->set(x[k++]);
			}
		}
	}

	void printInputs(std::ostream& os) const
	{
		os << "inputs= ";
		auto& nodes = primitives_.nodes(0);
		for (SizeType i = 0; i < nodes.size(); ++i) {
			if (!nodes[i]->isInput())
				continue;
			nodes[i]->print(os);
		}

		os << "\n";
	}

private:

	const PrimitivesType& primitives_;
	SizeType numberOfThreads_;
	SizeType inputs_;
};
}
#endif // INPUTGATESUTIL_H
