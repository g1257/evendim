/*
Copyright (c) 2017, UT-Battelle, LLC

evendim, Version 0.

This file is part of evendim.
evendim is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
evendim is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with evendim. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef GENE_H
#define GENE_H

#include "Vector.h"
#include "TypeToString.h"
#include "PsimagLite.h"
#include "ProgramGlobals.h"

namespace Gep {

template<typename TreeType,typename EvolutionType>
class Gene {

public:

	typedef typename EvolutionType::PrimitivesType::NodeType NodeType;
	typedef typename TreeType::VectorValueType VectorValueType;
	typedef typename NodeType::ValueType ValueType;
	typedef typename PsimagLite::Vector<TreeType*>::Type VectorTreeType;
	typedef typename PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef Gene<TreeType,EvolutionType> GeneType;

	Gene(const Gene& other)
	    : head_ (other.head_),
	      tail_(other.tail_),
	      vecStr_(other.vecStr_),
	      vt_(other.vt_.size(), nullptr)
	{
		const SizeType  n = vt_.size();
		for (SizeType i = 0; i < n; ++i) {
			vt_[i] = new TreeType(*other.vt_[i]);
		}
	}

	Gene(SizeType head,
	     bool isCell,
	     const EvolutionType& evolution,
	     const VectorStringType& vecStr,
	     SizeType threadNum = 0)
	    : head_(head),
	      tail_(evolution.tail(head)),
	      vecStr_(vecStr)
	{
		if (!isCell) evolution.checkStringNonCell(vecStr_, head);

		SizeType headPlusTail = head_ + tail_;

		fromString(vt_, evolution, vecStr, headPlusTail, isCell, threadNum);
	}

	~Gene()
	{
		deleteAll();
	}

	static void fromString(VectorTreeType& vt,
	                       const EvolutionType& evolution,
	                       const VectorStringType& vecStr,
	                       SizeType effectiveSize,
	                       bool isCell,
	                       SizeType threadNum)
	{
		PsimagLite::Vector<SizeType>::Type va;
		const SizeType dcLength = vecStr.size() - effectiveSize;
		SizeType sumOfA = 1;
		SizeType dcIndex = 0;
		PsimagLite::String dcStr = (dcLength > 0) ? vecStr[dcIndex + effectiveSize] : "0";
		assert(dcStr.length() == 1 and dcStr[0] >= 48);
		SizeType dcNumber = dcStr[0] - 48;
		const VectorValueType& dcArray = evolution.primitives().dcValues();
		assert(dcLength == 0 || dcNumber < dcArray.size());
		ValueType dcValue = (dcLength > 0) ? dcArray[dcNumber] : ValueType(0);
		const auto& nodes = evolution.primitives().nodes(threadNum);
		for (SizeType i = 0; i < effectiveSize; i++) {
			PsimagLite::String cStr = vecStr[i];
			const NodeType& node = ProgramGlobals::findNodeFromCode<NodeType>(cStr,
			                                                                  nodes,
			                                                                  dcValue,
			                                                                  isCell);
			if (cStr == "?") {
				assert(dcLength > 0);
				dcIndex++;
				assert(dcIndex < dcLength);
				dcStr = vecStr[dcIndex + effectiveSize];
				dcNumber = dcStr[0] - 48;
				assert(dcNumber < dcArray.size());
				dcValue = dcArray[dcNumber];
			}

			SizeType a = node.arity();
			sumOfA += (a - 1);
			TreeType* tree = new TreeType(evolution.primitives(),
			                              node,
			                              evolution.verbose());

			va.push_back(a);
			vt.push_back(tree);
			if (sumOfA == 0) break;
		}

		SizeType k = 0;
		for (SizeType i = 0; i < vt.size(); i++) {
			SizeType a = va[i];
			if (a == 0 || !vt[i]) continue;
			for (SizeType j = k+1; j < k+a+1; j++) {
				if (j>=vt.size()) continue;
				vt[i]->setDescendants(*vt[j]);
			}

			k += a;
		}
	}

	const VectorStringType& vecString() const
	{
		return vecStr_;
	}

	const TreeType& getExpression() const
	{
		return *vt_[0];
	}

	const SizeType head() const { return head_; }

	SizeType effectiveSize() const
	{
		return vt_.size();
	}

private:

	Gene& operator=(const Gene& other) = delete;

	void deleteAll()
	{
		for (SizeType i = 0; i < vt_.size(); i++) {
			if (vt_[i]) delete vt_[i];
			vt_[i] = 0;
		}
	}

	SizeType head_;
	SizeType tail_;
	VectorStringType vecStr_;
	VectorTreeType vt_;
}; // class Gene

} // namespace Gep

#endif // GENE_H
