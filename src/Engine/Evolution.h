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

#include "MersenneTwister.h"
#include "NodeHelper.hh"
#include "ProgramGlobals.h"
#include "PsimagLite.h"
#include "TypeToString.h"
#include "Vector.h"
#include <cassert>
#include <iostream>

namespace Gep {

template <typename PrimitivesType_>
class Evolution {

public:

	typedef typename PrimitivesType_::VectorNodeType VectorNodeType;
	typedef typename PrimitivesType_::NodeType NodeType;
	typedef typename PrimitivesType_::VectorValueType VectorValueType;
	typedef typename NodeType::ValueType ValueType;
	typedef typename PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef PrimitivesType_ PrimitivesType;
	using NodeHelperType = NodeHelper<NodeType>;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	Evolution(PrimitivesType& primitives,
	          SizeType r,
	          bool verbose)
	    : primitives_(primitives)
	    , verbose_(verbose)
	    , maxArity_(0)
	    , nodeHelper_(primitives.nodesSerial())
	    , rng_(r)
	{
		maxArity_ = nodeHelper_.maxArity();

		nodeHelper_.setInputsTerminalsAndNonTerminals(terminals_, nonTerminals_);
	}

	SizeType geneLength(SizeType head) const
	{
		bool hasDc = (primitives_.dcValues().size() > 0);
		SizeType dc = (hasDc) ? this->tail(head) : 0;

		return head + this->tail(head) + dc;
	}

	bool verbose() const { return verbose_; }

	SizeType tail(SizeType head) const
	{
		return head * (maxArity_ - 1) + 1;
	}

	VectorStringType randomGene(SizeType head) const
	{
		const SizeType tail1 = tail(head);
		VectorStringType str = nonTerminals_;
		ProgramGlobals::pushVector(str, terminals_);
		VectorStringType str1 = selectRandomFrom(head, str);
		const VectorStringType str2 = selectRandomFrom(tail1, terminals_);
		bool hasDc = (primitives_.dcValues().size() > 0);
		const SizeType dc = (hasDc) ? tail1 : 0;
		const VectorStringType str3 = selectRandomFrom(dc, primitives_.dcArray());
		ProgramGlobals::pushVector(str1, str2);
		ProgramGlobals::pushVector(str1, str3);
		return str1;
	}

	VectorStringType randomAdf(SizeType chead, SizeType genes) const
	{
		VectorStringType terminals(genes);
		for (SizeType i = 0; i < genes; i++)
			terminals[i] = ttos(i);
		SizeType ctail = tail(chead);
		VectorStringType str = nonTerminals_;
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

		SizeType index = static_cast<SizeType>(rng_() * len);

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
			return getStringForRegionCell(index, head, genes);
		else
			return getStringForRegionNonCell(index, head);
	}

	VectorStringType getStringForRegionNonCell(SizeType index,
	                                           SizeType head) const
	{
		SizeType tail1 = tail(head);
		if (index < head) {
			VectorStringType ret = terminals_;
			ProgramGlobals::pushVector(ret, nonTerminals_);
			return ret;
		}

		if (index < head + tail1)
			return terminals_;

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
			ProgramGlobals::pushVector(terminals, nonTerminals_);

		return terminals;
	}

	VectorStringType invert(const VectorStringType& str, SizeType head) const
	{
		VectorStringType ret = str;
		for (SizeType i = 0; i < head; ++i)
			ret[i] = str[head - i - 1];

		return ret;
	}

	void checkStringNonCell(const VectorStringType& vecStr,
	                        SizeType head,
	                        bool isCell) const
	{
		SizeType tail1 = tail(head);
		SizeType len = vecStr.size();
		bool hasDc = (primitives_.dcValues().size() > 0);
		SizeType dc = (hasDc) ? tail(head) : 0;

		if (len != head + tail1 + dc) {
			PsimagLite::String errorMessage(__FILE__);
			errorMessage += " " + ttos(__LINE__) + "\n";
			errorMessage += "Length= " + ttos(len);
			errorMessage += " expected " + ttos(head + tail1 + dc) + "\n";
			err(errorMessage);
		}

		VectorStringType terminals = terminals_;

		for (SizeType i = head; i < len - dc; i++) {
			if (std::find(terminals.begin(), terminals.end(), vecStr[i]) != terminals.end())
				continue;
			if (isCell && isAnInteger(vecStr[i]))
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

	NodeHelperType& nodeHelper() { return nodeHelper_; }

	const NodeHelperType& nodeHelper() const { return nodeHelper_; }

	double rng() const { return rng_(); }

	static bool isAnInteger(PsimagLite::String str)
	{
		SizeType len = str.length();
		for (SizeType i = 0; i < len; ++i) {
			if (str[i] < '0' || str[i] > '9')
				return false;
		}

		return true;
	}

private:

	VectorStringType selectRandomFrom(SizeType head, const VectorStringType& str) const
	{
		VectorStringType ret(head);

		for (SizeType i = 0; i < head; i++) {
			SizeType index = static_cast<SizeType>(rng_() * str.size());
			ret[i] = str[index];
		}

		return ret;
	}

	PrimitivesType& primitives_;
	bool verbose_;
	SizeType maxArity_;
	NodeHelperType nodeHelper_;
	VectorStringType nonTerminals_;
	VectorStringType terminals_;
	mutable PsimagLite::MersenneTwister rng_; // RandomForTests<double> rng_;
}; // class Evolution

} // namespace Gep
#endif // EVOLUTION_H
