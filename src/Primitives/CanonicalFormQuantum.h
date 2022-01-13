#ifndef CANONICALFORMQUANTUM_H
#define CANONICALFORMQUANTUM_H
#include "Vector.h"
#include "Node.h"
#include "ProgramGlobals.h"
#include "Sort.h"
#include <cassert>
#include "PsimagLite.h"
#include <queue>

namespace Gep {

template<typename ValueType_, typename AnglesType>
class CanonicalFormQuantum {

public:

	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef typename PsimagLite::Vector<ValueType_>::Type VectorValueType;
	typedef Node<VectorValueType, AnglesType> NodeType;
	typedef typename PsimagLite::Vector<NodeType*>::Type VectorNodeType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef std::queue<PsimagLite::String> QueueStringType;
	typedef PsimagLite::Vector<QueueStringType>::Type VectorQueueStringType;
	typedef PsimagLite::Vector<bool>::Type VectorBoolType;

	CanonicalFormQuantum(const VectorStringType& data,
	                     const VectorNodeType& nodes)
	    : data_(data), needsChange_(false)
	{
		getEffectiveAndJunk();
		needsChange_ |= orderGatesByBit(effective_, nodes);
	}

	void changeIfNeeded(VectorStringType& vstr) const
	{
		if (!needsChange_) return;
		vstr.resize(effective_.size() + junkDna_.size());
		for (SizeType i = 0; i < effective_.size(); ++i)
			vstr[i] = effective_[i];

		for (SizeType i = 0; i < junkDna_.size(); ++i)
			vstr[i + effective_.size()] = junkDna_[i];
	}

private:

	void getEffectiveAndJunk()
	{
		bool flag = true;
		for (auto it = data_.begin(); it != data_.end(); ++it) {
			if (flag)
				effective_.push_back(*it);
			else
				junkDna_.push_back(*it);

			if (*it == "0") flag = false;
		}
	}

	static bool orderGatesByBit(VectorStringType& effective, const VectorNodeType& nodes)
	{
		constexpr SizeType maxConstant = 12345;
		static const ValueType_ value;
		constexpr bool isCell = false;
		VectorQueueStringType v(100);
		VectorSizeType order(effective.size(), maxConstant);
		bool flag = false;

		for (auto it = effective.begin(); it != effective.end(); ++it) {
			const NodeType& node = ProgramGlobals::findNodeFromCode<NodeType>(*it,
			                                                                  nodes,
			                                                                  value,
			                                                                  isCell);

			if (node.isInput()) continue;

			if (node.arity() > 1) continue; // FIXME TODO more than one bit gates

			SizeType bit = getBit(node.code());
			assert(bit < v.size());
			v[bit].push(*it);
			order[it - effective.begin()] = bit;
			flag = true;
		}

		if (!flag) return false;

		VectorSizeType iperm(order.size());
		PsimagLite::Sort<VectorSizeType> sort;
		sort.sort(order, iperm);

		VectorStringType newData(effective.size());

		const SizeType n = effective.size();
		for (SizeType i = 0; i < n; ++i) {
			const SizeType bit = order[i]; // order is ordered
			if (bit >= maxConstant) {
				newData[i] = effective[i];
				continue;
			}

			assert(bit < v.size());
			assert(v[bit].size() > 0);
			newData[i] = v[bit].front();

			v[bit].pop();
		}

		const bool b = isEqual(newData, effective);
		if (b) return false;

		effective = newData;
		return true;
	}

	static SizeType getBit(PsimagLite::String code)
	{
		VectorStringType tokens;
		PsimagLite::split(tokens, code, ":");
		if (tokens.size() > 2)
			err("getBit: code " + code + " with two or more colons\n");

		if (tokens.size() == 2) code = tokens[0]; // ignore angles

		PsimagLite::String buffer;
		const SizeType n = code.size();
		for (SizeType i = 0; i < n; ++i) {
			const SizeType j = n - i - 1;
			unsigned char c = code[j];
			if (!std::isdigit(c)) break;
			buffer += c;
		}

		std::reverse(buffer.begin(), buffer.end());
		return PsimagLite::atoi(buffer);
	}

	static bool isEqual(const VectorStringType& v1, const VectorStringType& v2)
	{
		const SizeType n = v1.size();
		if (n != v2.size()) return false;
		for (SizeType i = 0; i < n; ++i)
			if (v1[i] != v2[i]) return false;
		return true;
	}

	const VectorStringType& data_;
	bool needsChange_;
	VectorStringType effective_;
	VectorStringType junkDna_;
};
}
#endif // CANONICALFORMQUANTUM_H
