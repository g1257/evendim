/*
Copyright (c) 2017-2021, UT-Battelle, LLC

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
#ifndef EVEDIM_QUANTUM_FITNESS_H
#define EVEDIM_QUANTUM_FITNESS_H
#include "BaseFitness.h"
#include "MersenneTwister.h"
#include "Minimizer.h"
#include "MinimizerParams.h"
#include "ProgramGlobals.h"
#include "PsimagLite.h"

namespace Gep {

template <typename ChromosomeType, typename EvolutionType, typename ComplexType>
class FunctionToMinimize {
public:

	typedef typename EvolutionType::ValueType QuasiVectorType;
	typedef typename PsimagLite::Real<ComplexType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef RealType FieldType;
	typedef typename ChromosomeType::VectorStringType VectorStringType;
	typedef typename PsimagLite::Vector<QuasiVectorType>::Type VectorQuasiVectorType;
	typedef typename EvolutionType::NodeHelperType::NodeFactoryType NodeFactoryType;

	enum class FunctionEnum { FITNESS,
		                  DIFFERENCE };

	FunctionToMinimize(EvolutionType& evolution,
	                   const ChromosomeType& chromosome,
	                   SizeType samples,
	                   SizeType threadNum)
	    : evolution_(evolution)
	    , chromosome_(chromosome)
	    , threadNum_(threadNum)
	    , vecVec_(samples)
	    , inVector_((1 << evolution.primitives().numberOfBits()))
	    , outVector_(inVector_.size())
	{
		// constexpr ComplexType two = 2;
		// constexpr ComplexType minusOne = -1;
		numberOfAngles_ = findNumberOfAngles(chromosome.effectiveVecString());
		PsimagLite::MersenneTwister rng(static_cast<int>(evolution_.rng() * 1000));
		for (SizeType i = 0; i < samples; ++i) {
			vecVec_[i].randomize(inVector_.size(), rng, 2., -1.);
		}
	}

	SizeType size() const { return numberOfAngles_; }

	RealType operator()(const VectorRealType& angles)
	{
		// the case without angles is handled elsewhere
		// if it reaches here, it's an internal error
		assert(numberOfAngles_ > 0);

		assert(angles.size() == numberOfAngles_);

		// false means don't be verbose here
		return fitness(&angles, FunctionEnum::DIFFERENCE, false);
	}

	void df(VectorRealType& dest, const VectorRealType& angles)
	{
		VectorStringType vecStr = chromosome_.vecString();
		encodeAngles(vecStr, angles);
		const ChromosomeType* chromosome = new ChromosomeType(chromosome_.params(),
		                                                      evolution_,
		                                                      vecStr,
		                                                      threadNum_);

		dest.resize(angles.size());

		const SizeType samples = vecVec_.size();
		for (SizeType angleIndex = 0; angleIndex < numberOfAngles_; ++angleIndex) {
			for (SizeType i = 0; i < samples; ++i) {
				setInVector(i);
				evolution_.nodeHelper().setInput(0, inVector_, threadNum_);
				functionF(outVector_, inVector_);

				computeDifferentialVector(differential_, angles, angleIndex);

				const RealType tmp = diffVectorDiff2(chromosome->exec(0),
				                                     outVector_,
				                                     differential_);
				dest[angleIndex] += tmp;
			}

			dest[angleIndex] /= samples;
		}
	}

	RealType fitness(const VectorRealType* angles, FunctionEnum functionEnum, bool verbose)
	{
		const ChromosomeType* chromosome = nullptr;

		VectorStringType vecStr = chromosome_.vecString();

		if (angles) {
			encodeAngles(vecStr, *angles);
			chromosome = new ChromosomeType(chromosome_.params(), evolution_, vecStr, threadNum_);
		}
		else {
			chromosome = &chromosome_;
		}

		const SizeType samples = vecVec_.size();
		RealType sum = 0;
		for (SizeType i = 0; i < samples; ++i) {
			setInVector(i);
			evolution_.nodeHelper().setInput(0, inVector_, threadNum_);
			if (verbose)
				evolution_.nodeHelper().printInputs(std::cout);

			functionF(outVector_, inVector_);

			const RealType tmp = vectorDiff2(chromosome->exec(0),
			                                 outVector_);
			sum += fabs(tmp);
		}

		if (angles) {
			delete chromosome;
			chromosome = nullptr;
		}

		sum /= samples;
		return (functionEnum == FunctionEnum::DIFFERENCE) ? sum : 1 - sum;
	}

	static void encodeAngles(VectorStringType& vecStr, const VectorRealType& angles)
	{
		const SizeType n = vecStr.size();
		SizeType currentIndex = 0;
		bool flag = true;
		for (SizeType i = 0; i < n; ++i) {

			if (isInputGate(vecStr[i])) {
				flag = false;
				continue;
			}

			if (numberOfAnglesOneGate(vecStr[i]) == 0 || !flag)
				continue;
			PsimagLite::String str = vecStr[i];
			str = NodeFactoryType::stripPreviousAngleIfAny(str);
			if (angles.size() < currentIndex)
				err("encodeAngles: too many angles for rotations in this individual!?\n");
			str += ":" + ttos(angles[currentIndex++]);
			vecStr[i] = str;
		}

		if (currentIndex != angles.size())
			err("encodeAngles: too few angles for rotations in this individual!?\n");
	}

	template <typename SomeRngType>
	static void initAngles(VectorRealType& angles, const VectorStringType& vStr, SomeRngType& rng)
	{
		const SizeType n = vStr.size();
		SizeType currentIndex = 0;
		for (SizeType i = 0; i < n; ++i) {
			if (numberOfAnglesOneGate(vStr[i]) == 0)
				continue;
			if (angles.size() < currentIndex)
				err("initAngles: too few angles for individual\n");

			initAngle(angles[currentIndex++], vStr[i], rng);
		}

		if (angles.size() != currentIndex)
			err("initAngles: too many angles for individual\n");
	}

	template <typename SomeRngType>
	static void initAngle(RealType& angle, PsimagLite::String str, SomeRngType& rng)
	{
		typename PsimagLite::String::const_iterator it = std::find(str.begin(),
		                                                           str.end(),
		                                                           ':');
		if (it == str.end()) {
			angle = 2 * M_PI * rng();
			return;
		}

		PsimagLite::String angleStr = str.substr(it - str.begin() + 1, str.end() - it - 1);
		angle = std::stod(angleStr);
	}

private:

	static SizeType findNumberOfAngles(const VectorStringType& vstr)
	{
		SizeType n = vstr.size();

		SizeType count = 0;
		for (SizeType i = 0; i < n; ++i) {
			count += numberOfAnglesOneGate(vstr[i]);
		}

		return count;
	}

	static SizeType numberOfAnglesOneGate(PsimagLite::String str)
	{
		if (str.length() == 0)
			return 0;
		return (str[0] == 'R' || str.substr(0, 2) == "PG") ? 1 : 0;
	}

	static bool isInputGate(PsimagLite::String str)
	{
		if (str.length() == 0)
			return false;
		return (str[0] == '0') ? true : false;
	}

	// This is the function to fit to
	static void functionF(QuasiVectorType& dest, const QuasiVectorType& src)
	{
		// flip the first bit
		QuasiVectorType dest2 = src;
		flipABit(dest2, src, 0);

		// flip the second bit
		flipABit(dest, dest2, 1);
	}

	// Flip the m-th bit
	// Example if m=0
	// 0.1*|0000> -0.2|1110>
	// src[0] = 0.1;   src[14] = -0.2 src[..] = 0
	// 0.1*|0001> - 0.2|1111>
	// dest[1] = 0.1;  dest[15] = -0.2 dest[...] = 0
	// Flip the first bit
	// 1111 <--- 15 --> i
	// 0001 <--- 1
	// 1110 <--- 14 --> j
	static void flipABit(QuasiVectorType& dest, const QuasiVectorType& src, SizeType bit)
	{
		dest.flipABit(src, bit);
	}

	void computeDifferentialVector(QuasiVectorType& differential,
	                               const VectorRealType& angles,
	                               SizeType angleIndex)
	{
		assert(angles.size() == numberOfAngles_);

		VectorStringType cString = chromosome_.effectiveVecString();

		SizeType geneLength = chromosome_.geneLength();

		VectorStringType tmpString = replaceOneR(cString, angleIndex, geneLength);

		assert(tmpString.size() == geneLength);

		// create derivative individual angle-th
		ChromosomeType newChromosome(chromosome_.params(), evolution_, tmpString, threadNum_);

		// apply to inVector
		differential = newChromosome.exec(0);
	}

	// adds padding as well
	VectorStringType replaceOneR(VectorStringType& v, SizeType angleIndex, SizeType geneLength)
	{
		VectorStringType w(geneLength);

		const SizeType n = v.size();
		assert(n <= geneLength);
		SizeType count = 0;
		for (SizeType i = 0; i < n; ++i) {
			w[i] = v[i];
			SizeType n = numberOfAnglesOneGate(v[i]);
			if (n == 0)
				continue;
			if (n == 1) {
				if (count++ == angleIndex)
					w[i] = "_" + v[i];
			}
		}

		assert(count == numberOfAngles_);

		for (SizeType i = n; i < geneLength; ++i)
			w[i] = "0";

		return w;
	}

	void setInVector(SizeType jnd)
	{
		const SizeType n = vecVec_[jnd].size();
		assert(inVector_.size() == n);
		inVector_ = vecVec_[jnd];
	}

	EvolutionType& evolution_;
	const ChromosomeType& chromosome_;
	SizeType threadNum_;
	SizeType numberOfAngles_;
	VectorQuasiVectorType vecVec_;
	QuasiVectorType inVector_;
	QuasiVectorType outVector_;
	QuasiVectorType differential_;
};

template <typename ChromosomeType>
class QuantumFitness : public BaseFitness<ChromosomeType> {

public:

	typedef BaseFitness<ChromosomeType> BaseType;
	typedef typename ChromosomeType::EvolutionType EvolutionType;
	typedef typename EvolutionType::PrimitivesType PrimitivesType;
	typedef typename PrimitivesType::ValueType VectorType;
	typedef typename VectorType::value_type ComplexType;
	typedef typename PsimagLite::Real<ComplexType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef MinimizerParams<RealType> MinimizerParamsType;
	typedef MinimizerParamsType FitnessParamsType;
	typedef FunctionToMinimize<ChromosomeType, EvolutionType, ComplexType>
	    FunctionToMinimizeType;
	typedef typename PsimagLite::Minimizer<RealType, FunctionToMinimizeType> MinimizerType;
	typedef typename ChromosomeType::VectorStringType VectorStringType;

	QuantumFitness(SizeType samples, EvolutionType& evolution, MinimizerParamsType* minParams)
	    : samples_(samples)
	    , evolution_(evolution)
	    , minParams_(*minParams)
	    , status_(0)
	{
		if (evolution.nodeHelper().numberOfInputs() != 1)
			err("QuantumFitness::ctor(): 1 input expected\n");
	}

	RealType getFitness(const ChromosomeType& chromosome,
	                    long unsigned int seed,
	                    SizeType threadNum)
	{
		if (threadNum > 0)
			err("QuantumFitness: Threading not supported yet (sorry)\n");

		FunctionToMinimizeType f(evolution_, chromosome, samples_, threadNum);

		if (f.size() == 0) {
			return f.fitness(nullptr,
			                 FunctionToMinimizeType::FunctionEnum::FITNESS,
			                 evolution_.verbose());
		}

		MinimizerType min(f, minParams_.maxIter, minParams_.verbose);

		int used = 0;
		VectorRealType angles(f.size());
		PsimagLite::MersenneTwister rng(seed);
		FunctionToMinimizeType::initAngles(angles, chromosome.effectiveVecString(), rng);
		if (minParams_.algo == MinimizerParamsType::SIMPLEX) {
			used = min.simplex(angles,
			                   minParams_.delta,
			                   minParams_.tol);
		}
		else if (minParams_.algo == MinimizerParamsType::NONE) {
			used = 1;
		}
		else {
			used = min.conjugateGradient(angles,
			                             minParams_.delta,
			                             minParams_.delta2,
			                             minParams_.tol,
			                             minParams_.saveEvery);
		}

		status_ = (min.status() == MinimizerType::GSL_SUCCESS) ? 0 : 1;

		if (status_ == 0) {
			VectorStringType vecStr = chromosome.vecString();
			FunctionToMinimizeType::encodeAngles(vecStr, angles);
			const ChromosomeType* chromosome2 = new ChromosomeType(chromosome.params(),
			                                                       evolution_,
			                                                       vecStr,
			                                                       threadNum);

			ChromosomeType& chromosomeNonconst = const_cast<ChromosomeType&>(chromosome);
			chromosomeNonconst = *chromosome2;

			delete chromosome2;
			chromosome2 = nullptr;
		}

		const bool printFooter = minParams_.verbose;
		// const int returnStatus = (used > 0) ? 0 : 1;

		RealType value = f.fitness(&angles,
		                           FunctionToMinimizeType::FunctionEnum::FITNESS,
		                           evolution_.verbose());

		if (!printFooter)
			return value; // <--- EARLY EXIT HERE

		std::cerr << "QuantumFitness::minimize(): ";
		if (min.status() == MinimizerType::GSL_SUCCESS) {
			std::cerr << " converged after ";
		}
		else {
			std::cerr << "NOT CONVERGED after ";
		}

		++used;
		std::cerr << used << " iterations. " << toString(angles) << "\n";

		return value;
	}

	RealType maxFitness() const { return samples_; }

private:

	static PsimagLite::String toString(const VectorRealType& angles)
	{
		const SizeType n = angles.size();
		PsimagLite::String str;
		for (SizeType i = 0; i < n; ++i)
			str += ttos(angles[i]) + " ";
		return str;
	}

	SizeType samples_;
	EvolutionType& evolution_;
	const MinimizerParamsType minParams_;
	int status_;
}; // class QuantumFitness
} // namespace Gep

#endif // EVEDIM_QUANTUM_FITNESS_H
