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
#ifndef EVOLUTION_H
#define EVOLUTION_H

#include "Vector.h"
#include <cassert>
#include <iostream>
#include "TypeToString.h"
#include "PsimagLite.h"
#include "ProgramGlobals.h"

namespace Gep {

template<typename PrimitivesType_>
class Evolution {

	typedef typename PrimitivesType_::VectorNodeType VectorNodeType;
	typedef typename PrimitivesType_::NodeType NodeType;
	typedef typename PrimitivesType_::VectorValueType VectorValueType;
	typedef typename NodeType::ValueType ValueType;
	typedef typename PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

public:

	typedef PrimitivesType_ PrimitivesType;

	Evolution(const PrimitivesType& primitives,
	          SizeType r,
	          bool verbose)
	    : primitives_(primitives),
	      verbose_(verbose)
	{
		const VectorNodeType& nodes = primitives_.nodes(0);

		for (SizeType i = 0; i < nodes.size(); i++) {
			if (nodes[i]->isInput()) {
				inputs_.push_back(nodes[i]);
			}
		}
	}

	bool verbose() const { return verbose_; }

	void setInput(SizeType i, ValueType x) const
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

	VectorStringType randomGene(SizeType head) const
	{
		const SizeType tail1 = tail(head);
		VectorStringType str = primitives_.nonTerminals();
		ProgramGlobals::pushVector(str, primitives_.terminals());
		VectorStringType str1 = selectRandomFrom(head,str);
		const VectorStringType str2 = selectRandomFrom(tail1, primitives_.terminals());
		const SizeType dc = (primitives_.hasDc()) ? tail1 : 0;
		const VectorStringType str3 = selectRandomFrom(dc, primitives_.dcArray());
		ProgramGlobals::pushVector(str1, str2);
		ProgramGlobals::pushVector(str1, str3);
		return str1;
	}

	VectorStringType randomAdf(SizeType chead,SizeType genes) const
	{
		VectorStringType terminals(genes);
		for (SizeType i = 0; i < genes; i++)
			terminals[i] = ttos(i);
		SizeType ctail = tail(chead);
		VectorStringType str = primitives_.nonTerminals();
		ProgramGlobals::pushVector(str, terminals);
		VectorStringType str1 = selectRandomFrom(chead, str);
		const VectorStringType str2 = selectRandomFrom(ctail, terminals);
		ProgramGlobals::pushVector(str1, str2);
		return str1;
	}

	const PrimitivesType& primitives() const { return primitives_; }

	VectorStringType mutate(const VectorStringType& str,
	                        SizeType head,
	                        SizeType genes,
	                        bool isCell) const
	{
		SizeType len = str.size();

		SizeType index = static_cast<SizeType>(primitives_.rng() * len);

		VectorStringType something = getStringForRegion(index, head, genes, isCell);

		VectorStringType ret(index);
		for (SizeType i = 0; i < index; ++i)
			ret[i] = str[i];

		VectorStringType vRandom = selectRandomFrom(1, something);
		ProgramGlobals::pushVector(ret, vRandom);

		for (SizeType i = index + 1; i < str.size(); ++i)
			ret.push_back(str[i]);

		return ret;
	}

	VectorStringType getStringForRegion(SizeType index,
	                                    SizeType head,
	                                    SizeType genes,
	                                    bool isCell) const
	{
		if (isCell)
			return getStringForRegionCell(index,head,genes);
		else
			return getStringForRegionNonCell(index,head);
	}

	VectorStringType getStringForRegionNonCell(SizeType index,
	                                           SizeType head) const
	{
		SizeType tail1 = tail(head);
		if (index < head) {
			VectorStringType ret = primitives_.terminals();
			ProgramGlobals::pushVector(ret, primitives_.nonTerminals());
			return ret;
		}

		if (index < head + tail1)
			return primitives_.terminals();

		return primitives_.dcArray();
	}

	VectorStringType getStringForRegionCell(SizeType index,
	                                        SizeType head,
	                                        SizeType genes) const
	{

		VectorStringType terminals(genes);
		for (SizeType i = 0; i < genes; ++i)
			terminals[i] = ttos(i);

		if (index < head)
			ProgramGlobals::pushVector(terminals, primitives_.nonTerminals());

		return terminals;
	}

	VectorStringType invert(const VectorStringType& str,SizeType head) const
	{
		VectorStringType ret = str;
		for (SizeType i = 0; i < head; ++i)
			ret[i] = str[head - i - 1];

		return ret;
	}

	void checkStringNonCell(const VectorStringType& vecStr,
	                        SizeType head) const
	{
		SizeType tail1 = tail(head);
		SizeType len = vecStr.size();
		SizeType dc = (primitives_.hasDc())? tail(head) : 0;

		if (len != head + tail1 + dc) {
			PsimagLite::String errorMessage(__FILE__);
			errorMessage += " " + ttos(__LINE__) + "\n";
			errorMessage += "Length= " + ttos(len);
			errorMessage += " expected " + ttos(head + tail1 + dc) + "\n";
			err(errorMessage);
		}

		VectorStringType terminals = primitives_.terminals();

		for (SizeType i = head; i < len -dc; i++) {
			if (std::find(terminals.begin(),terminals.end(), vecStr[i]) != terminals.end())
				continue;
			PsimagLite::String errorMessage(__FILE__);
			errorMessage += " " + ttos(__LINE__) + "\n";
			errorMessage += "head= " + ttos(head);
			errorMessage += " string " + ProgramGlobals::vecStrToStr(vecStr, "") + "\n";
			throw PsimagLite::RuntimeError(errorMessage);
		}

		for (SizeType i = head + tail1; i < len; ++i) {
			PsimagLite::String cStr = vecStr[i];
			if (cStr.size() != 1) {
				PsimagLite::String errorMessage(__FILE__);
				errorMessage += " " + ttos(__LINE__) + "\n";
				err(errorMessage);
			}

			const char c = cStr[0];
			SizeType x = c - 48;
			if (x >= primitives_.dcArray().size()) {
				PsimagLite::String errorMessage(__FILE__);
				errorMessage += " " + ttos(__LINE__) + "\n";
				err(errorMessage);
			}
		}
	}

	void checkStringCell(const VectorStringType& str,
	                     SizeType head,
	                     SizeType genes) const
	{
		SizeType tail1 = tail(head);
		SizeType len = str.size();

		if (len != head + tail1) {
			PsimagLite::String errorMessage(__FILE__);
			errorMessage += " " + ttos(__LINE__) + "\n";
			errorMessage += "Length= " + ttos(len);
			errorMessage += " expected " + ttos(head + tail1) + "\n";
			throw PsimagLite::RuntimeError(errorMessage);
		}

		VectorStringType terminals(genes);
		for (SizeType i = 0; i < genes; ++i)
			terminals[i] = ttos(i);

		for (SizeType i = head; i < len; i++) {
			if (std::find(terminals.begin(), terminals.end(), str[i]) != terminals.end())
				continue;
			PsimagLite::String errorMessage(__FILE__);
			errorMessage += " " + ttos(__LINE__) + "\n";
			errorMessage += "head= " + ttos(head);
			errorMessage += " string " + ProgramGlobals::vecStrToStr(str, "") + "\n";
			err(errorMessage);
		}
	}

private:

	VectorStringType selectRandomFrom(SizeType head, const VectorStringType& str) const
	{
		VectorStringType ret(head);

		for (SizeType i = 0; i < head; i++) {
			SizeType index = static_cast<SizeType>(primitives_.rng()*str.size());
			ret[i] = str[index];
		}

		return ret;
	}

	const PrimitivesType& primitives_;
	bool verbose_;
	VectorNodeType inputs_;

}; // class Evolution

} // namespace Gep
#endif // EVOLUTION_H
