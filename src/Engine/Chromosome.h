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
#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include "Gene.h"
#include "ProgramGlobals.h"
#include "Vector.h"

namespace Gep {

template <typename TreeType, typename EvolutionType_, typename ParametersType>
class Chromosome {

public:

	typedef EvolutionType_ EvolutionType;
	typedef Gene<TreeType, EvolutionType> GeneType;
	typedef typename EvolutionType::PrimitivesType PrimitivesType;
	typedef typename PrimitivesType::ValueType ValueType;
	typedef typename TreeType::NodeType NodeType;
	typedef typename NodeType::AnglesType AnglesType;
	typedef typename PsimagLite::Vector<AnglesType>::Type VectorAnglesType;
	typedef typename PsimagLite::Vector<ValueType>::Type VectorValueType;
	typedef typename PsimagLite::Vector<GeneType*>::Type VectorGeneType;
	typedef typename GeneType::VectorStringType VectorStringType;
	typedef Chromosome<TreeType, EvolutionType, ParametersType> ChromosomeType;
	typedef std::pair<VectorStringType, VectorStringType> PairVectorStringType;
	typedef typename PsimagLite::Vector<VectorStringType>::Type VectorVectorStringType;
	typedef std::pair<VectorVectorStringType, VectorVectorStringType> PairVectorVectorStringType;

	Chromosome(const ParametersType& params,
	           const EvolutionType& evolution,
	           const VectorStringType& vecStr,
	           SizeType threadNum)
	    : evolution_(evolution)
	    , params_(params)
	    , geneLength_(evolution.geneLength(params.head))
	    , isLinearTree_(true)
	{
		SizeType len = vecStr.size();
		if (len == 0)
			throw PsimagLite::RuntimeError("Chromosome::ctor()\n");

		SizeType index = 0;
		VectorStringType buffer;
		for (SizeType i = 0; i < geneLength_; i++)
			buffer.push_back(" ");

		for (SizeType i = 0; i < len; i++) {
			buffer[index] = vecStr[i];
			index++;
			if (index == geneLength_) {
				index = 0;
				GeneType* gene = new GeneType(params.head,
				                              false,
				                              evolution,
				                              buffer,
				                              threadNum);

				// All genes must be linear for the chromosome to be linear
				if (isLinearTree_ && !gene->isLinearTree())
					isLinearTree_ = false;

				genes_.push_back(gene);

				if (genes_.size() == params.genes)
					break;
			}
		}

		assert(genes_.size() == params.genes);

		index = 0;

		SizeType start = geneLength_ * genes_.size();
		SizeType cgeneLength = params.chead + evolution.tail(params.chead);

		buffer.clear();
		for (SizeType i = 0; i < cgeneLength; i++)
			buffer.push_back(" ");

		for (SizeType i = start; i < len; i++) {
			buffer[index] = vecStr[i];
			index++;
			if (index == cgeneLength) {
				index = 0;
				GeneType* gene = new GeneType(params.chead,
				                              true,
				                              evolution,
				                              buffer,
				                              threadNum);
				adfs_.push_back(gene);
				if (adfs_.size() == params.adfs)
					break;
			}
		}

		for (SizeType i = 0; i < genes_.size(); i++) {
			const SizeType effectiveSize = genes_[i]->effectiveSize();
			ProgramGlobals::pushVector(effectiveVecStr_, genes_[i]->vecString(), effectiveSize);
		}

		assert(adfs_.size() == params.adfs);
		if (adfs_.size() == 0)
			return;

		isLinearTree_ = false; // disable linear tree when there are ADFs
		assert(adfs_.size() == 1);

		adfsVecStr_ = adfs_[0]->vecString();
		const SizeType adfsEffective = adfs_[0]->effectiveSize();
		const VectorStringType& adfsVec = adfs_[0]->vecString();

		ProgramGlobals::pushVector(effectiveVecStr_, adfsVec, adfsEffective);
	}

	~Chromosome()
	{
		for (SizeType i = 0; i < genes_.size(); i++)
			delete genes_[i];
		genes_.clear();

		for (SizeType i = 0; i < adfs_.size(); i++)
			delete adfs_[i];
		adfs_.clear();
	}

	Chromosome& operator=(const Chromosome& other)
	{
		//		evolution_ = other.evolution_;

		//		params_ = other.params_;

		effectiveVecStr_ = other.effectiveVecStr_;

		adfsVecStr_ = other.adfsVecStr_;

		isLinearTree_ = other.isLinearTree_;

		const SizeType n = genes_.size();
		for (SizeType i = 0; i < n; ++i) {
			delete genes_[i];
			genes_[i] = new GeneType(*other.genes_[i]);
		}

		const SizeType m = adfs_.size();
		for (SizeType i = 0; i < m; ++i) {
			delete adfs_[i];
			adfs_[i] = new GeneType(*other.adfs_[i]);
		}

		return *this;
	}

	bool isLinearTree() const { return isLinearTree_; }

	SizeType geneLength() const { return geneLength_; }

	VectorStringType vecString() const
	{
		VectorStringType ret;

		for (SizeType i = 0; i < genes_.size(); i++) {
			ProgramGlobals::pushVector(ret, genes_[i]->vecString());
		}

		ProgramGlobals::pushVector(ret, adfsVecStr_);

		return ret;
	}

	VectorStringType vecString(SizeType i) const
	{
		if (i < genes_.size())
			return genes_[i]->vecString();

		SizeType index = i - genes_.size();
		assert(index < adfs_.size());
		return adfs_[index]->vecString();
	}

	const VectorStringType& effectiveVecString() const
	{
		return effectiveVecStr_;
	}

	ValueType exec(SizeType outputIndex) const
	{
		assert(genes_.size() > 0);
		assert(outputIndex < genes_.size());

		VectorValueType values(genes_.size());
		for (SizeType i = 0; i < values.size(); i++) {
			values[i] = genes_[i]->exec();
		}

		if (adfs_.size() == 0)
			return values[outputIndex];

		PsimagLite::String msg("Chromosome::exec(): ");
		if (outputIndex > 0)
			throw PsimagLite::RuntimeError(msg + "outputIndex>0 only with adfs==0\n");
		// reinterpret inputs
		for (SizeType i = 0; i < adfs_.size(); i++) {
			adfs_[i]->set(values);
		}

		if (adfs_.size() != 1)
			throw PsimagLite::RuntimeError(msg + "adfs must be 1\n");

		ValueType tmp = adfs_[0]->exec();

		return tmp;
	}

	SizeType effectiveSize() const { return effectiveVecStr_.size(); }

	SizeType size() const { return genes_.size(); }

	PairVectorStringType recombine(const ChromosomeType& other,
	                               SizeType points) const
	{
		SizeType genes = genes_.size();
		SizeType index = static_cast<SizeType>(evolution_.rng() * (genes + adfs_.size()));
		GeneType* gene = (index >= genes) ? adfs_[index - genes] : genes_[index];

		bool isCell = (index >= genes);
		SizeType indexCorrected = (isCell) ? genes_.size() : index;

		VectorStringType firstVec;
		for (SizeType i = 0; i < indexCorrected; i++) {
			ProgramGlobals::pushVector(firstVec, genes_[i]->vecString());
		}

		VectorStringType lastVec;
		if (!isCell) {
			for (SizeType i = index + 1; i < genes_.size(); i++) {
				ProgramGlobals::pushVector(lastVec, genes_[i]->vecString());
			}

			ProgramGlobals::pushVector(lastVec, adfsVecStr_);
		}

		PairVectorStringType p = (points == 1) ? recombine1(gene->vecString(), other.vecString(index)) : recombine2(gene->vecString(), other.vecString(index));

		VectorStringType vecStr1 = firstVec;
		ProgramGlobals::pushVector(vecStr1, p.first);
		ProgramGlobals::pushVector(vecStr1, lastVec);

		VectorStringType vecStr2 = firstVec;
		ProgramGlobals::pushVector(vecStr2, p.second);
		ProgramGlobals::pushVector(vecStr2, lastVec);

		return PairVectorStringType(vecStr1, vecStr2);
	}

	VectorStringType evolve(const PsimagLite::String& action) const
	{
		SizeType genes = genes_.size();
		SizeType index = static_cast<SizeType>(evolution_.rng() * (genes + adfs_.size()));

		GeneType* gene = (index >= genes) ? adfs_[index - genes] : genes_[index];
		bool isCell = (index >= genes);

		VectorStringType firstVec;
		for (SizeType i = 0; i < index; i++) {
			ProgramGlobals::pushVector(firstVec, genes_[i]->vecString());
		}

		VectorStringType lastVec;

		if (!isCell) {
			for (SizeType i = index + 1; i < genes_.size(); ++i) {
				ProgramGlobals::pushVector(lastVec, genes_[i]->vecString());
			}

			ProgramGlobals::pushVector(lastVec, adfsVecStr_);
		}

		VectorStringType ret = firstVec;
		if (action == "mutate") {
			ProgramGlobals::pushVector(ret, evolution_.mutate(gene->vecString(), gene->head(), genes, isCell));
			ProgramGlobals::pushVector(ret, lastVec);
			return ret;
		}
		else if (action == "invert") {
			ProgramGlobals::pushVector(ret, evolution_.invert(gene->vecString(), gene->head()));
			ProgramGlobals::pushVector(ret, lastVec);
			return ret;
		}
		else if (action == "swap") {
			ProgramGlobals::pushVector(ret, swap(gene->vecString(), gene->head(), isCell));
			ProgramGlobals::pushVector(ret, lastVec);
			return ret;
		}

		throw PsimagLite::RuntimeError("Chromosome::evolve()\n");
	}

	const ParametersType& params() const { return params_; }

private:

	Chromosome(const Chromosome&) = delete;

	VectorStringType swap(const VectorStringType& str,
	                      SizeType head,
	                      bool isCell) const
	{
		assert(head > 0);
		VectorStringType ret = str;
		SizeType index = head;
		while (index + 1 >= head) {
			index = static_cast<SizeType>(evolution_.rng() * str.size());
		}

		ret[index] = str[index + 1];
		ret[index + 1] = str[index];

		if (isCell)
			evolution_.checkStringCell(ret, head, genes_.size());
		else
			evolution_.checkStringNonCell(ret, head, false);
		return ret;
	}

	PairVectorStringType recombine1(const VectorStringType& str1,
	                                const VectorStringType& str2) const
	{
		assert(str1.size() == str2.size());

		SizeType len = str1.size();
		SizeType index = static_cast<SizeType>(evolution_.rng() * len);
		PairVectorStringType newVecStrings;
		newVecStrings.first = recombine(str1, str2, index);
		newVecStrings.second = recombine(str2, str1, index);
		return newVecStrings;
	}

	VectorStringType recombine(const VectorStringType& str1,
	                           const VectorStringType& str2,
	                           SizeType index) const
	{
		VectorStringType ret(str2.size());
		for (SizeType i = 0; i < index; ++i)
			ret[i] = str1[i];

		for (SizeType i = index; i < str2.size(); ++i)
			ret[i] = str2[i];

		return ret;
	}

	VectorStringType recombine(const VectorStringType& str1,
	                           const VectorStringType& str2,
	                           SizeType index1,
	                           SizeType index2) const
	{
		VectorStringType tmp1(str1.size());
		for (SizeType i = 0; i < index1; ++i)
			tmp1[i] = str1[i];

		for (SizeType i = index1; i < index2; ++i)
			tmp1[i] = str2[i];

		for (SizeType i = index2; i < str1.size(); ++i)
			tmp1[i] = str1[i];

		return tmp1;
	}

	PairVectorStringType recombine2(const VectorStringType& str1,
	                                const VectorStringType& str2) const
	{
		assert(str1.size() == str2.size());

		SizeType len = str1.size();

		SizeType i1 = static_cast<SizeType>(evolution_.rng() * len);
		SizeType i2 = static_cast<SizeType>(evolution_.rng() * len);
		SizeType index1 = (i1 < i2) ? i1 : i2;
		SizeType index2 = (i1 < i2) ? i2 : i1;

		PairVectorStringType newVecStrings;
		newVecStrings.first = recombine(str1, str2, index1, index2);
		assert(newVecStrings.first.size() == len);
		newVecStrings.second = recombine(str2, str1, index1, index2);
		assert(newVecStrings.second.size() == len);
		return newVecStrings;
	}

	const EvolutionType& evolution_;
	const ParametersType& params_;
	SizeType geneLength_;
	bool isLinearTree_;
	VectorStringType effectiveVecStr_;
	VectorStringType adfsVecStr_;
	VectorGeneType genes_;
	VectorGeneType adfs_;
}; // class Fitness

} // namespace Gep

#endif // CHROMOSOME_H
