#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include "Vector.h"
#include "Gene.h"

namespace Gep {

template<typename TreeType,typename EvolutionType,typename ParametersType>
class Chromosome {

	typedef Gene<TreeType,EvolutionType> GeneType;
	typedef typename EvolutionType::PrimitivesType PrimitivesType;
	typedef typename PrimitivesType::ValueType ValueType;
	typedef typename PsimagLite::Vector<ValueType>::Type VectorValueType;
	typedef typename PsimagLite::Vector<GeneType*>::Type VectorGeneType;
	typedef Chromosome<TreeType,EvolutionType,ParametersType> ChromosomeType;
	typedef std::pair<PsimagLite::String,PsimagLite::String> PairOfStringsType;

public:

	Chromosome(const ParametersType& params,
	           const EvolutionType& evolution,
	           const PsimagLite::String& str)
	    : evolution_(evolution),
	      params_(params),
	      effectiveString_(""),
	      adfsString_("")
	{
		SizeType len = str.length();
		SizeType dc = (evolution_.primitives().hasDc())? evolution.tail(params.head) : 0;

		SizeType geneLength = params.head + evolution.tail(params.head) + dc;

		if (len == 0)
			throw PsimagLite::RuntimeError("Chromosome::ctor()\n");

		SizeType index = 0;
		PsimagLite::String buffer = "";
		for (SizeType i = 0; i < geneLength; i++)
			buffer += " ";

		for (SizeType i = 0; i < len; i++) {
			buffer[index] = str[i];
			index++;
			if (index == geneLength) {
				index = 0;
				GeneType* gene = new GeneType(params.head,false,evolution,buffer);
				genes_.push_back(gene);
				if (genes_.size() == params.genes) break;
			}
		}

		assert(genes_.size() == params.genes);

		index = 0;

		SizeType start = geneLength*genes_.size();
		SizeType cgeneLength = params.chead + evolution.tail(params.chead);

		buffer = "";
		for (SizeType i = 0; i < cgeneLength; i++)
			buffer += " ";

		for (SizeType i = start; i < len; i++) {
			buffer[index] = str[i];
			index++;
			if (index == cgeneLength) {
				index = 0;
				GeneType* gene = new GeneType(params.chead,true,evolution,buffer);
				adfs_.push_back(gene);
				if (adfs_.size() == params.adfs) break;
			}
		}

		for (SizeType i = 0; i < genes_.size(); i++) {
			effectiveString_ += genes_[i]->effectiveString();
		}

		assert(adfs_.size() == params.adfs);
		if (adfs_.size() == 0) return;
		assert(adfs_.size() == 1);

		adfsString_ = adfs_[0]->string();
		PsimagLite::String adfsEffective = adfs_[0]->effectiveString();

		effectiveString_ += adfsEffective;

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

	PsimagLite::String string(const PsimagLite::String& sep = "") const
	{
		PsimagLite::String ret = "";

		for (SizeType i = 0; i < genes_.size(); i++) {
			ret += genes_[i]->string() + sep;
		}

		return ret + adfsString_;
	}

	PsimagLite::String string(SizeType i) const
	{
		if (i < genes_.size())
			return genes_[i]->string();

		SizeType index = i - genes_.size();
		assert(index < adfs_.size());
		return adfs_[index]->string();
	}

	const PsimagLite::String& effectiveString() const
	{
		return effectiveString_;
	}

	ValueType exec(SizeType outputIndex) const
	{
		assert(genes_.size() > 0);
		assert(outputIndex < genes_.size());
		VectorValueType values(genes_.size());
		for (SizeType i = 0; i < values.size(); i++) {
			values[i] = genes_[i]->getExpression().exec();
		}

		if (adfs_.size() == 0) return values[outputIndex];

		PsimagLite::String msg("Chromosome::exec(): ");
		if (outputIndex > 0)
			throw PsimagLite::RuntimeError(msg + "outputIndex>0 only with adfs==0\n");
		// reinterpret inputs
		for (SizeType i = 0; i < adfs_.size(); i++) {
			adfs_[i]->getExpression().set(values);
		}

		if (adfs_.size() != 1)
			throw PsimagLite::RuntimeError(msg + "adfs must be 1\n");

		ValueType tmp = adfs_[0]->getExpression().exec();

		return tmp;
	}

	SizeType effectiveSize() const { return effectiveString_.size(); }

	SizeType size() const { return genes_.size(); }

	PairOfStringsType recombine(const ChromosomeType& other,
	                            SizeType points) const
	{
		const PrimitivesType& primitives = evolution_.primitives();
		SizeType genes = genes_.size();
		SizeType index = static_cast<SizeType>(primitives.rng() * (genes+adfs_.size()));
		GeneType *gene = (index >= genes) ? adfs_[index - genes] : genes_[index];

		bool isCell = (index >= genes);
		SizeType indexCorrected = (isCell) ? genes_.size() : index;

		PsimagLite::String first = "";
		for (SizeType i = 0; i < indexCorrected; i++) {
			first += genes_[i]->string();
		}

		PsimagLite::String last = "";
		if (!isCell) {
			for (SizeType i = index+1; i < genes_.size(); i++) {
				last += genes_[i]->string();
			}
			last += adfsString_;
		}

		PairOfStringsType p = (points == 1) ?
		            recombine1(gene->string(),other.string(index)) :
		            recombine2(gene->string(),other.string(index));

		PsimagLite::String str1 = first + p.first + last;
		PsimagLite::String str2 = first + p.second + last;

		return PairOfStringsType(str1,str2);
	}

	PsimagLite::String evolve(const PsimagLite::String& action) const
	{
		const PrimitivesType& primitives = evolution_.primitives();

		SizeType genes = genes_.size();
		SizeType index = static_cast<SizeType>(primitives.rng() * (genes+adfs_.size()));

		GeneType *gene = (index >= genes) ? adfs_[index - genes] : genes_[index];
		bool isCell = (index >= genes);

		PsimagLite::String first = "";
		for (SizeType i = 0; i < index; i++) {
			first += genes_[i]->string();
		}

		PsimagLite::String last = "";

		if (!isCell) {
			for (SizeType i = index+1; i < genes_.size(); i++) {
				last += genes_[i]->string();
			}
			last +=  adfsString_;
		}

		if (action == "mutate") {
			return first + evolution_.mutate(gene->string(),
			                                 gene->head(),
			                                 genes,
			                                 isCell) + last;
		} else if (action == "invert") {
			return first + evolution_.invert(gene->string(),
			                                 gene->head()) + last;
		} else if (action == "swap") {
			return first + swap(gene->string(),
			                    gene->head(),isCell) + last;
		}
		throw PsimagLite::RuntimeError("Chromosome::evolve()\n");
	}

private:

	PsimagLite::String swap(const PsimagLite::String& str,
	                        SizeType head,
	                        bool isCell) const
	{
		PsimagLite::String ret = str;
		const PrimitivesType& primitives = evolution_.primitives();
		SizeType index = static_cast<SizeType>(primitives.rng() * str.length());
		if (index+1 == head) index = 0;
		if (index+1 == str.length())
			index--;
		if (index+1 == head + evolution_.tail(head)) index = 0;
		ret[index] = str[index+1];
		ret[index+1] = str[index];

		if (isCell) evolution_.checkStringCell(ret,head,genes_.size());
		else evolution_.checkStringNonCell(ret,head);
		return ret;
	}

	PairOfStringsType recombine1(const PsimagLite::String& str1,
	                            const PsimagLite::String& str2) const
	{
		assert(str1.length() == str2.length());

		SizeType len = str1.length();
		const PrimitivesType& primitives = evolution_.primitives();
		SizeType index = static_cast<SizeType>(primitives.rng() * len);
		PairOfStringsType newStrings;
		newStrings.first = recombine(str1,str2,index);
		newStrings.second = recombine(str2,str1,index);
		return newStrings;
	}

	PsimagLite::String recombine(const PsimagLite::String& str1,
	                            const PsimagLite::String& str2,
	                            SizeType index) const
	{
		return str1.substr(0,index) + str2.substr(index);
	}

	PsimagLite::String recombine(const PsimagLite::String& str1,
	                             const PsimagLite::String& str2,
	                             SizeType index1,
	                             SizeType index2) const
	{
		PsimagLite::String tmp1 = str1.substr(0,index1);
		PsimagLite::String tmp2 = str2.substr(index1,index2-index1);
		return tmp1 + tmp2 + str1.substr(index2);
	}

	PairOfStringsType recombine2(const PsimagLite::String& str1,
	                            const PsimagLite::String& str2) const
	{
		assert(str1.length() == str2.length());

		SizeType len = str1.length();

		const PrimitivesType& primitives = evolution_.primitives();
		SizeType i1 = static_cast<SizeType>(primitives.rng() * len);
		SizeType i2 = static_cast<SizeType>(primitives.rng() * len);
		SizeType index1 = (i1 < i2) ? i1 : i2;
		SizeType index2 = (i1 < i2) ? i2 : i1;

		PairOfStringsType newStrings;
		newStrings.first = recombine(str1,str2,index1,index2);
		assert(newStrings.first.length() == len);
		newStrings.second = recombine(str2,str1,index1,index2);
		assert(newStrings.second.length() == len);
		return newStrings;
	}

	const EvolutionType& evolution_;
	const ParametersType& params_;
	PsimagLite::String effectiveString_;
	PsimagLite::String adfsString_;
	VectorGeneType genes_;
	VectorGeneType adfs_;

}; // class Fitness

} // namespace Gep

#endif // CHROMOSOME_H
