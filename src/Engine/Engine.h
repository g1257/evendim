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
#ifndef ENGINE_H
#define ENGINE_H

#include "Tree.h"
#include "Chromosome.h"
#include "ParametersEngine.h"
#include "Sort.h"
#include "Parallelizer2.h"

namespace Gep {

template<template<typename> class FitnessTemplate, typename EvolutionType>
class Engine {

public:

	typedef typename EvolutionType::PrimitivesType PrimitivesType;
	typedef typename PrimitivesType::ValueType ValueType;
	typedef typename PrimitivesType::CanonicalFormType CanonicalFormType;
	typedef Tree<PrimitivesType> TreeType;
	typedef double RealType;
	typedef ParametersEngine<RealType> ParametersEngineType_;
	typedef Chromosome<TreeType,EvolutionType,ParametersEngineType_> ChromosomeType;
	typedef FitnessTemplate<ChromosomeType> FitnessType;
	typedef typename FitnessType::FitnessParamsType FitnessParamsType;
	typedef typename PsimagLite::Vector<ChromosomeType*>::Type VectorChromosomeType;
	typedef typename PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename ChromosomeType::PairVectorStringType PairVectorStringType;
	typedef typename ChromosomeType::PairVectorVectorStringType PairVectorVectorStringType;
	typedef typename ChromosomeType::VectorVectorStringType VectorVectorStringType;
	typedef typename ChromosomeType::VectorAnglesType VectorAnglesType;
	typedef PsimagLite::Vector<long unsigned int>::Type VectorLongUnsignedType;
	typedef ParametersEngineType_ ParametersEngineType;

	/* PSIDOC Engine::ctor
The engine constructor creates the initial individuals randomly.
	 */
	Engine(const ParametersEngineType& params,
	       EvolutionType& evolution,
	       FitnessParamsType* fitnessParams = nullptr)
	    : params_(params),
	      evolution_(evolution),
	      fitness_(params.samples, evolution, fitnessParams)
	{
		constexpr SizeType threadNum = 0;
		for (SizeType i = 0; i< params_.population; ++i) {
			VectorStringType vecStr;
			for (SizeType j = 0; j < params_.genes; ++j)
				ProgramGlobals::pushVector(vecStr, evolution_.randomGene(params_.head));
			for (SizeType j = 0; j < params_.adfs; ++j)
				ProgramGlobals::pushVector(vecStr, evolution_.randomAdf(params_.chead,
				                                                        params_.genes));
			ChromosomeType* chromosome = new ChromosomeType(params_,
			                                                evolution_,
			                                                vecStr,
			                                                threadNum);
			chromosomes_.push_back(chromosome);
		}
	}

	~Engine()
	{
		deleteAll();
	}

	/* PSIDOC Engine::evolve
Engine::evolve() function starts by considering all parent chromosomes.
It then computes the fitness of these parent chromosomes.
It then applies one-point recombination, two-point recombination, mutation, inversion,
and swap algorithms to all parent chromosomes to generate the descendants for this generation.
It then canonicalizes them and selects the best p chromosomes and discards the ones with lowest
fitness, where p is the population number set from the input file or the command line.
	 */
	bool evolve(SizeType ind)
	{
		PairVectorVectorStringType newChromosomes;
		VectorRealType parentFitness(chromosomes_.size());
		SizeType totalChromosomes = chromosomes_.size();
		for (SizeType i = 0; i < totalChromosomes; i++) {
			const VectorStringType vecStr = chromosomes_[i]->vecString();

			const VectorStringType& effectiveVec = chromosomes_[i]->effectiveVecString();
			if (notAdded(newChromosomes.second, effectiveVec)) {
				newChromosomes.first.push_back(vecStr);
				newChromosomes.second.push_back(effectiveVec);
			}
		}

		PsimagLite::CodeSectionParams codeParams = PsimagLite::Concurrency::codeSectionParams;
		codeParams.npthreads = std::min(totalChromosomes,
		                                PsimagLite::Concurrency::codeSectionParams.npthreads);

		VectorLongUnsignedType seeds = fitness_.createSeeds(totalChromosomes);
		PsimagLite::Parallelizer2<> parallelizer2(codeParams);
		parallelizer2.parallelFor(0,
		                          totalChromosomes,
		                          [&parentFitness, &seeds, this](SizeType ind, SizeType threadNum) {
			parentFitness[ind] = -fitness_.getFitness(*chromosomes_[ind], seeds[ind], threadNum);
		});

		evolution_.nodeFactory().sync();

		recombination(newChromosomes, parentFitness, 1);

		recombination(newChromosomes, parentFitness, 2);

		evolve(newChromosomes, "mutate");

		evolve(newChromosomes,"invert");

		evolve(newChromosomes,"swap");

		if (ind > 0 && !params_.options.isSet("noncanonical"))
			canonicalizeAll(newChromosomes.first);

		return selectBest(newChromosomes.first);
	}

private:

	void canonicalizeAll(VectorVectorStringType& newChromosomes)
	{
		for (SizeType i = 0; i < newChromosomes.size(); i++) {
			CanonicalFormType canonicalForm(newChromosomes[i],
			                                evolution_.nodeFactory());
			canonicalForm.changeIfNeeded(newChromosomes[i]);
		}
	}

	void deleteAll()
	{
		for (SizeType i = 0; i < chromosomes_.size(); i++)
			delete chromosomes_[i];
		chromosomes_.clear();
	}

	void addWithCare(PairVectorVectorStringType& newChromosomes,
	                 const VectorStringType& vecStr) const
	{
		constexpr SizeType threadNum = 0;
		ChromosomeType chromosome(params_, evolution_, vecStr, threadNum);
		VectorStringType vEff(chromosome.effectiveSize());
		for (SizeType i = 0; i < chromosome.effectiveSize(); ++i)
			vEff[i] = chromosome.vecString()[i];
		if (notAdded(newChromosomes.second, vEff))
			newChromosomes.first.push_back(vecStr);
	}

	void recombination(PairVectorVectorStringType& newChromosomes,
	                   const VectorRealType& parentFitness,
	                   SizeType points) const
	{
		for (SizeType i = 0; i < params_.descendants; i++) {
			SizeType index1 = selectAccordingToFitness(parentFitness);
			SizeType index2 = selectAccordingToFitness(parentFitness);
			PairVectorStringType newStrings = chromosomes_[index1]->
			        recombine(*chromosomes_[index2],
			                  points);

			addWithCare(newChromosomes, newStrings.first);

			addWithCare(newChromosomes, newStrings.second);
		}
	}

	SizeType selectAccordingToFitness(const VectorRealType& parentFitness) const
	{
		RealType minFitness = -1e50;
		for (SizeType i = 0; i < parentFitness.size(); i++)
			if (minFitness < parentFitness[i])
				minFitness = parentFitness[i];

		minFitness = -minFitness;

		RealType totalFitness = 0;
		for (SizeType i = 0; i < parentFitness.size(); i++)
			totalFitness += (-parentFitness[i] - minFitness);

		if (totalFitness < 0)
			throw PsimagLite::RuntimeError("totalFitness<=0");

		if (totalFitness == 0)
			return static_cast<SizeType>(evolution_.rng() * parentFitness.size());

		RealType r = evolution_.rng()*totalFitness;
		RealType min = 0;
		RealType max = 0;
		for (SizeType i = 0; i < parentFitness.size(); i++) {
			max += (-parentFitness[i] - minFitness);
			if (r <= max && r >= min) return i;
			min = max;
		}

		throw PsimagLite::RuntimeError("selectAccordingToFitness\n");
	}

	void evolve(PairVectorVectorStringType& newChromosomes,
	            const PsimagLite::String& action) const
	{
		SizeType population = chromosomes_.size();
		for (SizeType i = 0; i < params_.mutation; i++) {
			SizeType index = static_cast<SizeType>(fitness_.rng() * population);
			VectorStringType newVecStr = chromosomes_[index]->evolve(action);

			addWithCare(newChromosomes, newVecStr);
		}
	}

	bool selectBest(VectorVectorStringType& newChromosomes)
	{
		assert(chromosomes_.size() > 0);
		typename PsimagLite::Vector<RealType>::Type fitness(newChromosomes.size());

		if (newChromosomes.size() < chromosomes_.size()) {
			PsimagLite::String errorMessage("selectBest, newChromosomes= ");
			errorMessage += ttos(newChromosomes.size()) + " oldChromosomes= ";
			errorMessage += ttos(chromosomes_.size());
			throw PsimagLite::RuntimeError(errorMessage);
		}

		const SizeType totalChromosomes = newChromosomes.size();
		PsimagLite::CodeSectionParams codeParams = PsimagLite::Concurrency::codeSectionParams;
		codeParams.npthreads = std::min(totalChromosomes,
		                       PsimagLite::Concurrency::codeSectionParams.npthreads);

		assert(codeParams.npthreads > 0);
		bool withProgressBar = (codeParams.npthreads > 1) ? false
		                                                  : params_.options.isSet("progressBar");

		bool isVerbose = (evolution_.verbose() && codeParams.npthreads == 1);
		VectorLongUnsignedType seeds = fitness_.createSeeds(totalChromosomes);

		PsimagLite::Parallelizer2<> parallelizer2(codeParams);
		parallelizer2.parallelFor(0,
		                          totalChromosomes,
		                          [&newChromosomes,
		                          &fitness,
		                          &seeds,
		                          isVerbose,
		                          withProgressBar,
		                          this](SizeType ind, SizeType threadNum) {
			ChromosomeType chromosome(params_, evolution_, newChromosomes[ind], threadNum);
			if (isVerbose)
				std::cout<<"About to exec chromosome= "<<newChromosomes[ind]<<"\n";
			fitness[ind] = -fitness_.getFitness(chromosome, seeds[ind], threadNum);
			newChromosomes[ind] = chromosome.vecString();
			const int status = fitness_.status();
			const PsimagLite::String symbol = (status == 0) ? "." : "*";
			if (withProgressBar) std::cerr<<symbol;
		});

		evolution_.nodeFactory().sync();

		if (withProgressBar) std::cerr<<"\n";

		PsimagLite::Sort<typename PsimagLite::Vector<RealType>::Type> sort;
		PsimagLite::Vector<SizeType>::Type iperm(fitness.size());
		sort.sort(fitness,iperm);

		VectorVectorStringType newChromosomes2 = newChromosomes;
		for (SizeType i = 0; i < newChromosomes.size(); i++)
			newChromosomes[i] = newChromosomes2[iperm[i]];

		orderBySize(newChromosomes, fitness);

		SizeType population = chromosomes_.size();
		RealType fraction = 0.8;
		SizeType point = static_cast<SizeType>(population*fraction);

		deleteAll();

		RealType maxFitness = params_.samples;
		for (SizeType i = 0; i < point; i++) {
			RealType f = -fitness[i];
			addChromosome(newChromosomes[i], f);
			if (i==0 && f == maxFitness) return true;
		}

		for (SizeType i = point; i < population; i++) {
			SizeType index = point +
			        static_cast<SizeType>(fitness_.rng() * population * (1.0-fraction));
			assert(index >= point);
			addChromosome(newChromosomes[index],-fitness[index]);
		}

		std::cout<<"----------------\n";
		return false;
	}

	void addChromosome(const VectorStringType& str, const RealType& f)
	{
		constexpr SizeType threadNum = 0;
		ChromosomeType* chromosome = new ChromosomeType(params_,
		                                                evolution_,
		                                                str,
		                                                threadNum);

		assert(chromosome);

		chromosomes_.push_back(chromosome);

		std::cout<<ProgramGlobals::vecStrToStr(chromosome->vecString(), " ");
		const auto fit = (params_.options.isSet("printcompact")) ? " fit " : " fitness ";
		std::cout<<fit<<f<<" "<<fitness_.info(*chromosome);
		const auto esize = (params_.options.isSet("printcompact")) ? " #= " : " effective size= ";
		std::cout<<esize<<chromosome->effectiveSize()<<"\n";
	}

	bool notAdded(const VectorVectorStringType& newChromosomes,
	              const VectorStringType& newStr) const
	{
		return (find(newChromosomes.begin(),
		             newChromosomes.end(),
		             newStr) == newChromosomes.end());
	}

	SizeType getWeightedIndex(VectorSizeType& added,
	                          const VectorRealType& fitness,
	                          SizeType populationOver2) const
	{
		ValueType r = evolution_.rng();
		ValueType f = (1.0-exp(-r))/(1.0-exp(-1.0));
		SizeType index = 0;
		do {
			index = findIndexWithFitness(f,fitness,populationOver2,added);
		} while (find(added.begin(),added.end(),index) != added.end());

		added.push_back(index);
		return index;
	}

	SizeType findIndexWithFitness(const ValueType& f,
	                              const VectorRealType& fitness,
	                              SizeType populationOver2,
	                              const VectorSizeType& added) const
	{
		SizeType population = 2 * populationOver2;
		SizeType index = 0;
		ValueType value = f*fitness_.maxFitness();
		ValueType min = getMax(fitness,value);

		for (SizeType i = populationOver2; i < population; i++) {
			if (find(added.begin(),added.end(),i) != added.end())
				continue;
			if (fabs(-fitness[i]-value) < min) {
				min = fabs(-fitness[i]-value);
				index = i;
			}
		}
		assert(index >= populationOver2);
		return index;
	}

	ValueType getMax(const VectorRealType& fitness,const ValueType& value) const
	{
		ValueType max = 0;
		for (SizeType i = 0; i < fitness.size(); i++) {
			if (fabs(-fitness[i]-value) > max) {
				max = fabs(-fitness[i]-value);
			}
		}
		return max;
	}

	void orderBySize(VectorVectorStringType& newChromosomes,const VectorRealType& fitness) const
	{
		constexpr SizeType threadNum = 0;
		RealType value = -fitness_.maxFitness();
		VectorRealType bestSize;

		for (SizeType i = 0; i < fitness.size(); i++) {
			if (fitness[i] != value) break;
			ChromosomeType chromosome(params_, evolution_, newChromosomes[i], threadNum);
			bestSize.push_back(chromosome.effectiveSize());

		}
		if (bestSize.size() == 0) return;

		PsimagLite::Vector<SizeType>::Type iperm(bestSize.size());
		PsimagLite::Sort<typename PsimagLite::Vector<RealType>::Type> sort;
		sort.sort(bestSize,iperm);

		VectorVectorStringType oldChromosomes = newChromosomes;
		for (SizeType i = 0; i < fitness.size(); i++) {
			if (fitness[i] != value) break;
			newChromosomes[i] = oldChromosomes[iperm[i]];
		}

	}

	const ParametersEngineType& params_;
	EvolutionType& evolution_;
	FitnessType fitness_;
	VectorChromosomeType chromosomes_;
}; // class Engine

} // namespace Gep
#endif // ENGINE_H
