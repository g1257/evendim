#ifndef GENE_H
#define GENE_H

#include "Vector.h"
#include "TypeToString.h"

namespace Gep {

template<typename TreeType,typename EvolutionType>
class Gene {

	typedef typename EvolutionType::PrimitivesType::NodeType NodeType;
	typedef typename TreeType::VectorValueType VectorValueType;
	typedef typename NodeType::ValueType ValueType;
	typedef typename PsimagLite::Vector<TreeType*>::Type VectorTreeType;
	typedef Gene<TreeType,EvolutionType> GeneType;

public:

	Gene(SizeType head,
	     bool isCell,
	     const EvolutionType& evolution,
	     const PsimagLite::String& str)
	    : head_(head),
	      tail_(evolution.tail(head)),
	      str_(str)
	{
		if (!isCell) evolution.checkStringNonCell(str_,head);

		SizeType headPlusTail = head_ + tail_;

		fromString(vt_,evolution,str,headPlusTail,isCell);
	}

	~Gene()
	{
		deleteAll();
	}

	static void fromString(VectorTreeType& vt,
	                       const EvolutionType& evolution,
	                       const PsimagLite::String& str,
	                       SizeType effectiveSize,
	                       bool isCell)
	{
		PsimagLite::Vector<SizeType>::Type va;
		PsimagLite::String dc = str.substr(effectiveSize);

		SizeType sumOfA = 1;
		SizeType dcIndex = 0;
		char dcChar = (dc.length() > 0) ? dc[dcIndex] : '0';
		assert(dcChar >= 48);
		SizeType dcNumber = dcChar - 48;
		const VectorValueType& dcArray = evolution.primitives().dcValues();
		assert(dc.length() ==0 || dcNumber < dcArray.size());
		ValueType dcValue = (dc.length() > 0) ? dcArray[dcNumber] : 0;
		for (SizeType i = 0; i < effectiveSize; i++) {
			char c = str[i];
			const NodeType& node = evolution.findNodeWithCode(c,dcValue,isCell);
			if (c == '?') {
				assert(dc.length() > 0);
				dcIndex++;
				assert(dcIndex < dc.length());
				dcChar = dc[dcIndex];
				dcNumber = dcChar - 48;
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

	const PsimagLite::String& string() const
	{
		return str_;
	}

	const TreeType& getExpression() const
	{
		return *vt_[0];
	}

	const SizeType head() const { return head_; }

	PsimagLite::String effectiveString() const
	{
		SizeType index = vt_.size();
		return str_.substr(0,index);
	}

private:

	void deleteAll()
	{
		for (SizeType i = 0; i < vt_.size(); i++) {
			if (vt_[i]) delete vt_[i];
			vt_[i] = 0;
		}
	}

	SizeType head_;
	SizeType tail_;
	PsimagLite::String str_;
	VectorTreeType vt_;
}; // class Gene

} // namespace Gep

#endif // GENE_H
