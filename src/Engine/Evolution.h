#ifndef EVOLUTION_H
#define EVOLUTION_H

#include "Vector.h"
#include <cassert>
#include <iostream>
#include "TypeToString.h"

namespace Gep {

template<typename PrimitivesType_>
class Evolution {

	typedef typename PrimitivesType_::VectorNodeType VectorNodeType;
	typedef typename PrimitivesType_::NodeType NodeType;
	typedef typename PrimitivesType_::VectorValueType VectorValueType;
	typedef typename NodeType::ValueType ValueType;

public:

	typedef PrimitivesType_ PrimitivesType;

	Evolution(const PrimitivesType& primitives,
	          SizeType r,
	          bool verbose)
	    : primitives_(primitives),
	      verbose_(verbose)
	{
		const VectorNodeType& nodes = primitives_.nodes();

		for (SizeType i = 0; i < nodes.size(); i++) {
			if (nodes[i]->isInput()) {
				assert(inputs_.size() == static_cast<SizeType>(nodes[i]->code() - 48));
				inputs_.push_back(nodes[i]);
			}
		}
	}

	bool verbose() const { return verbose_; }

	const NodeType& findNodeWithCode(char c,
	                                 const ValueType& value = 0,
	                                 bool isCell = false) const
	{
		const VectorNodeType& nodes = primitives_.nodes();

		for (SizeType i = 0; i < nodes.size(); i++) {
			if (isCell && nodes[i]->isInput()) continue;
			if (nodes[i]->code() == c) {
				if (c == '?') nodes[i]->setDcValue(value);
				return *nodes[i];
			}
		}

		throw PsimagLite::RuntimeError("findNodeWithCode\n");
	}

	void setInput(SizeType i,ValueType x) const
	{
		assert(i < inputs_.size());
		inputs_[i]->set(x);
	}

	void setInput(const VectorValueType& x) const
	{
		assert(x.size() == inputs_.size());
		SizeType n = std::min(x.size(),inputs_.size());
		for (SizeType i = 0; i < n; ++i)
			inputs_[i]->set(x[i]);
	}

	void printInputs(std::ostream& os) const
	{
		os<<"inputs= ";
		for (SizeType i = 0; i < inputs_.size(); i++)
			inputs_[i]->print(os);
		os<<"\n";
	}

	SizeType inputs() const { return inputs_.size(); }

	SizeType tail(SizeType head) const
	{
		return head*(primitives_.arity() - 1) + 1;
	}

	PsimagLite::String randomGene(SizeType head) const
	{
		SizeType tail1 = tail(head);
		PsimagLite::String str = primitives_.nonTerminals() + primitives_.terminals();
		PsimagLite::String str1 = selectRandomFrom(head,str);
		PsimagLite::String str2 = selectRandomFrom(tail1,primitives_.terminals());
		SizeType dc = (primitives_.hasDc()) ? tail1 : 0;
		PsimagLite::String str3 = selectRandomFrom(dc,primitives_.dcArray());
		return str1 + str2 + str3;
	}

	PsimagLite::String randomAdf(SizeType chead,SizeType genes) const
	{
		PsimagLite::String terminals = "";
		for (SizeType i = 0; i < genes; i++)
			terminals += ttos(i);
		SizeType ctail = tail(chead);
		PsimagLite::String str = primitives_.nonTerminals() + terminals;
		PsimagLite::String str1 = selectRandomFrom(chead,str);
		PsimagLite::String str2 = selectRandomFrom(ctail,terminals);
		return str1 + str2;
	}

	const PrimitivesType& primitives() const { return primitives_; }

	PsimagLite::String mutate(const PsimagLite::String& str,
	                          SizeType head,
	                          SizeType genes,
	                          bool isCell) const
	{
		SizeType len = str.length();

		SizeType index = static_cast<SizeType>(primitives_.rng() * len);

		PsimagLite::String something = getStringForRegion(index,head,genes,isCell);

		return str.substr(0,index) + selectRandomFrom(1,something) + str.substr(index+1);

	}

	PsimagLite::String getStringForRegion(SizeType index,
	                                      SizeType head,
	                                      SizeType genes,
	                                      bool isCell) const
	{
		if (isCell)
			return getStringForRegionCell(index,head,genes);
		else
			return getStringForRegionNonCell(index,head);
	}

	PsimagLite::String getStringForRegionNonCell(SizeType index,
	                                             SizeType head) const
	{
		SizeType tail1 = tail(head);
		if (index < head)
			return primitives_.terminals() + primitives_.nonTerminals();

		if (index < head + tail1)
			return primitives_.terminals();

		return primitives_.dcArray();
	}

	PsimagLite::String getStringForRegionCell(SizeType index,
	                                             SizeType head,
	                                             SizeType genes) const
	{

		PsimagLite::String terminals = "";
		for (SizeType i = 0; i < genes; i++)
			terminals += ttos(i);

		if (index < head)
			return terminals + primitives_.nonTerminals();

		return terminals;
	}

	PsimagLite::String invert(const PsimagLite::String& str,SizeType head) const
	{
		PsimagLite::String ret = str;
		for (SizeType i = 0; i < head; i++) {
			ret[i] = str[head-i-1];
		}
		return ret;
	}

	void checkStringNonCell(const PsimagLite::String& str,
	                        SizeType head) const
	{
		SizeType tail1 = tail(head);
		SizeType len = str.length();
		SizeType dc = (primitives_.hasDc())? tail(head) : 0;

		if (len != head + tail1 + dc) {
			PsimagLite::String errorMessage(__FILE__);
			errorMessage += " " + ttos(__LINE__) + "\n";
			errorMessage += "Length= " + ttos(len);
			errorMessage += " expected " + ttos(head + tail1 + dc) + "\n";
			throw PsimagLite::RuntimeError(errorMessage);
		}

		PsimagLite::String terminals = primitives_.terminals();

		for (SizeType i = head; i < len -dc; i++) {
			if (terminals.find(str[i]) != PsimagLite::String::npos) continue;
			PsimagLite::String errorMessage(__FILE__);
			errorMessage += " " + ttos(__LINE__) + "\n";
			errorMessage += "head= " + ttos(head);
			errorMessage += " string " + str + "\n";
			throw PsimagLite::RuntimeError(errorMessage);
		}

		for (SizeType i = head+tail1; i < len; i++) {
			char c = str[i];
			SizeType x = c - 48;
			if (x >= primitives_.dcArray().length()) {
				PsimagLite::String errorMessage(__FILE__);
				errorMessage += " " + ttos(__LINE__) + "\n";
				throw PsimagLite::RuntimeError(errorMessage);
			}
		}
	}

	void checkStringCell(const PsimagLite::String& str,
	                     SizeType head,
	                     SizeType genes) const
	{
		SizeType tail1 = tail(head);
		SizeType len = str.length();

		if (len != head + tail1) {
			PsimagLite::String errorMessage(__FILE__);
			errorMessage += " " + ttos(__LINE__) + "\n";
			errorMessage += "Length= " + ttos(len);
			errorMessage += " expected " + ttos(head + tail1) + "\n";
			throw PsimagLite::RuntimeError(errorMessage);
		}

		PsimagLite::String terminals = "";
		for (SizeType i = 0; i < genes; i++)
			terminals += ttos(i);

		for (SizeType i = head; i < len; i++) {
			if (terminals.find(str[i]) != PsimagLite::String::npos) continue;
			PsimagLite::String errorMessage(__FILE__);
			errorMessage += " " + ttos(__LINE__) + "\n";
			errorMessage += "head= " + ttos(head);
			errorMessage += " string " + str + "\n";
			throw PsimagLite::RuntimeError(errorMessage);
		}
	}

private:

	PsimagLite::String selectRandomFrom(SizeType head,const PsimagLite::String& str) const
	{
		PsimagLite::String ret = "";

		for (SizeType i = 0; i < head; i++) {
			SizeType index = static_cast<SizeType>(primitives_.rng()*str.length());
			ret += str[index];
		}

		return ret;
	}

	const PrimitivesType& primitives_;
	bool verbose_;
	VectorNodeType inputs_;

}; // class Evolution

} // namespace Gep
#endif // EVOLUTION_H
