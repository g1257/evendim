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
#ifndef EVENDIM_GROUND_STATE_FITNESS_H
#define EVENDIM_GROUND_STATE_FITNESS_H
#include "PsimagLite.h"
#include "Minimizer.h"
#include "MinimizerParams.h"
#include "BaseFitness.h"
#include "MersenneTwister.h"
#include "GroundStateParams.h"
#include "Hamiltonian.h"

namespace Gep {

template<typename ChromosomeType, typename EvolutionType, typename GroundStateParamsType>
class FunctionToMinimize2 {
public:

	typedef typename GroundStateParamsType::ComplexType ComplexType;
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorType;
	typedef typename PsimagLite::Real<ComplexType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef RealType FieldType;
	typedef typename ChromosomeType::VectorStringType VectorStringType;
	typedef PsimagLite::Matrix<ComplexType> MatrixType;
	typedef typename EvolutionType::NodeFactoryType NodeFactoryType;

	enum class FunctionEnum {FITNESS, DIFFERENCE};

	FunctionToMinimize2(EvolutionType& evolution,
	                    const ChromosomeType& chromosome,
	                    const GroundStateParamsType& groundStateParams,
	                    SizeType thread)
	    : evolution_(evolution),
	      chromosome_(chromosome),
	      groundStateParams_(groundStateParams),
	      outVector_(groundStateParams_.inVector.size()),
	      threadNum_(thread)
	{
		numberOfAngles_ = findNumberOfAngles(chromosome.effectiveVecString());
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
		evolution_.setInput(0, groundStateParams_.inVector, threadNum_);

		const VectorType& inVector = groundStateParams_.inVector;
		for (SizeType angleIndex = 0; angleIndex < numberOfAngles_; ++angleIndex) {
			evolution_.setInput(0, inVector, threadNum_);

			computeDifferentialVector(differential_, angles, angleIndex);

			const RealType tmp = diffVectorDiff2(chromosome->exec(0),
			                                     outVector_,
			                                     differential_);
			dest[angleIndex] += tmp;
		}
	}

	RealType fitness(const VectorRealType* angles, FunctionEnum functionEnum, bool verbose)
	{
		const ChromosomeType* chromosome = nullptr;

		VectorStringType vecStr = chromosome_.vecString();

		if (angles) {
			encodeAngles(vecStr, *angles);
			chromosome = new ChromosomeType(chromosome_.params(),
			                                evolution_,
			                                vecStr,
			                                threadNum_);
		} else {
			chromosome = &chromosome_;
		}

		evolution_.setInput(0, groundStateParams_.inVector, threadNum_);
		if (verbose) evolution_.printInputs(std::cout);

		// oracle goes here
		RealType e = groundStateParams_.hamiltonian.energy(chromosome->exec(0), threadNum_);

		if (angles) {
			delete chromosome;
			chromosome = nullptr;
		}

		return (functionEnum == FunctionEnum::DIFFERENCE) ? e : -e;
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

			if (numberOfAnglesOneGate(vecStr[i]) == 0 || !flag) continue;
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

	static void initAngles(VectorRealType& angles,
	                       const VectorStringType& vStr,
	                       long unsigned int seed)
	{
		const SizeType n = vStr.size();
		SizeType currentIndex = 0;
		PsimagLite::MersenneTwister rng(seed);
		for (SizeType i = 0; i < n; ++i) {
			if (numberOfAnglesOneGate(vStr[i]) == 0) continue;
			if (angles.size() < currentIndex)
				err("initAngles: too few angles for individual\n");

			initAngle(angles[currentIndex++], vStr[i], rng);
		}

		if (angles.size() != currentIndex)
			err("initAngles: too many angles for individual\n");
	}

private:

	template<typename SomeRngType>
	static void initAngle(RealType& angle,
	                      PsimagLite::String str,
	                      SomeRngType& rng)
	{
		typename PsimagLite::String::const_iterator it = std::find(str.begin(),
		                                                           str.end(),
		                                                           ':');
		if (it == str.end()) {
			angle = 2*M_PI*rng();
			return;
		}

		PsimagLite::String angleStr = str.substr(it - str.begin() + 1, str.end() - it - 1);
		angle = std::stod(angleStr);
	}

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
		if (str.length() == 0) return 0;
		return (str[0] == 'R' || str.substr(0, 2) == "PG") ? 1 : 0;
	}

	static bool isInputGate(PsimagLite::String str)
	{
		if (str.length() == 0) return false;
		return (str[0] == '0') ? true : false;
	}

	static RealType vectorDiff2(const VectorType& v1, const VectorType& v2)
	{
		const SizeType n = v1.size();
		assert(n == v2.size());
		RealType sum = 0;
		for (SizeType i = 0; i < n; ++i)
			sum += std::abs(v2[i] - v1[i]);

		return sum/n;
	}

	static RealType diffVectorDiff2(const VectorType& v1,
	                                const VectorType& v2,
	                                const VectorType& v3)
	{
		const SizeType n = v1.size();
		assert(n == v2.size());
		assert(n == v3.size());
		RealType sum = 0;
		for (SizeType i = 0; i < n; ++i) {
			RealType denom = std::abs(v2[i] - v1[i]);
			if (denom == 0) denom = 1;
			const RealType re1 = PsimagLite::real(v2[i] - v1[i]);
			const RealType im1 = PsimagLite::imag(v2[i] - v1[i]);
			sum += (PsimagLite::real(v3[i])*re1 + PsimagLite::imag(v3[i])*im1)/denom;
		}

		return sum/n;
	}

	void computeDifferentialVector(VectorType& differential,
	                               const VectorRealType& angles,
	                               SizeType angleIndex)
	{
		differential.resize(angles.size());
		std::fill(differential.begin(), differential.end(), 0);

		assert(angles.size() == numberOfAngles_);

		VectorStringType cString = chromosome_.effectiveVecString();

		VectorStringType tmpString = replaceOneR(cString, angleIndex, chromosome_.vecString(0).size());

		assert(tmpString.size() == chromosome_.vecString(0).size());

		// create derivative individual angle-th
		ChromosomeType newChromosome(chromosome_.params(), evolution_, tmpString, threadNum_);

		// apply to inVector
		evolution_.setInput(0, groundStateParams_.inVector, threadNum_);
		differential = newChromosome.exec(0);
	}

	// adds padding as well
	VectorStringType replaceOneR(VectorStringType& v, SizeType angleIndex, SizeType fullLength)
	{
		VectorStringType w(fullLength);

		const SizeType n = v.size();
		assert(n <= fullLength);
		SizeType count = 0;
		for (SizeType i = 0; i < n; ++i) {
			w[i] = v[i];
			SizeType m = numberOfAnglesOneGate(v[i]);
			if (m == 0) continue;
			if (m == 1) {
				if (count++ == angleIndex)
					w[i] = "_" + v[i];
			}
		}

		assert(count == numberOfAngles_);

		for (SizeType i = n; i < fullLength; ++i)
			w[i] = "0";

		return w;
	}

	EvolutionType& evolution_;
	const ChromosomeType& chromosome_;
	const GroundStateParamsType& groundStateParams_;
	SizeType numberOfAngles_;
	VectorType outVector_;
	VectorType differential_;
	SizeType threadNum_;
};

/* PSIDOC GroundStateFitnessClass
When Runtype=``GroundState'' in the input file, quantumGEP finds
circuits that, when applied to the initial state, yield the ground
state of a chosen Hamiltonian. Class GroundStateFitness implements
the fitness for this case. The fitness of an individual equals minus the value of
$\langle v|H|v\rangle,$ where $|v\rangle$ is the vector produced
by the individual (that is, the quantum circuit) when applied to the
initial state, and $H$ is the Hamiltonian.
*/
template<typename ChromosomeType>
class GroundStateFitness : public BaseFitness<ChromosomeType> {

public:

	typedef typename ChromosomeType::EvolutionType EvolutionType;
	typedef BaseFitness<ChromosomeType> BaseType;
	typedef typename EvolutionType::PrimitivesType PrimitivesType;
	typedef typename PrimitivesType::ValueType VectorType;
	typedef typename VectorType::value_type ComplexType;
	typedef typename PsimagLite::Real<ComplexType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef Hamiltonian<ComplexType> HamiltonianType;
	typedef GroundStateParams<HamiltonianType, ComplexType> GroundStateParamsType;
	typedef typename GroundStateParamsType::MinimizerParamsType MinimizerParamsType;

	typedef GroundStateParamsType FitnessParamsType;

	GroundStateFitness(SizeType samples,
	                  EvolutionType& evolution,
	                  FitnessParamsType* fitParams)
	    : evolution_(evolution),
	      fitParams_(*fitParams)
	{
		if (evolution.numberOfInputs() != 1)
			err("QuantumOracle::ctor(): 1 input expected\n");
		if (samples != 1)
			err("Expecting samples == 1\n");

		const SizeType threadNum = 0;
		evolution.setInput(0, fitParams->inVector, threadNum);
	}

	RealType getFitness(const ChromosomeType& chromosome,
	                    long unsigned int seed,
	                    SizeType threadNum)
	{
		typedef FunctionToMinimize2<ChromosomeType, EvolutionType, GroundStateParamsType>
		        FunctionToMinimizeType;
		typedef typename PsimagLite::Minimizer<RealType, FunctionToMinimizeType> MinimizerType;
		typedef typename ChromosomeType::VectorStringType VectorStringType;

		evolution_.setInput(0, fitParams_.inVector, threadNum);

		RealType norma = PsimagLite::norm(fitParams_.inVector);
		if (fabs(norma - 1) > 1e-4) err("Input vector not normalized\n");

		FunctionToMinimizeType f(evolution_, chromosome, fitParams_, threadNum);

		if (f.size() == 0) {
			return f.fitness(nullptr,
			                 FunctionToMinimizeType::FunctionEnum::FITNESS,
			                 evolution_.verbose());
		}

		const MinimizerParamsType& minParams = fitParams_.minParams;
		MinimizerType min(f, minParams.maxIter, minParams.verbose);

		int used = 0;
		VectorRealType angles(f.size());
		FunctionToMinimizeType::initAngles(angles, chromosome.effectiveVecString(), seed);
		if (minParams.algo == MinimizerParamsType::SIMPLEX) {
			used = min.simplex(angles,
			                   minParams.delta,
			                   minParams.tol);
		} else if (minParams.algo == MinimizerParamsType::NONE) {
			used = 1;
		} else {
			used = min.conjugateGradient(angles,
			                             minParams.delta,
			                             minParams.delta2,
			                             minParams.tol,
			                             minParams.saveEvery);
		}

		int status = (min.status() == MinimizerType::GSL_SUCCESS) ? 0 : 1;

		if (status == 0) {
			VectorStringType vecStr = chromosome.vecString();
			FunctionToMinimizeType::encodeAngles(vecStr, angles);
			const ChromosomeType* chromosome2 = new ChromosomeType(chromosome.params(),
			                                                       evolution_,
			                                                       vecStr,
			                                                       threadNum);

			ChromosomeType* chromosomeNonconst = const_cast<ChromosomeType*>(&chromosome);
			*chromosomeNonconst = *chromosome2;

			delete chromosome2;
			chromosome2 = nullptr;
		}

		const bool printFooter = minParams.verbose;
		//const int returnStatus = (used > 0) ? 0 : 1;

		RealType value = f.fitness(&angles,
		                           FunctionToMinimizeType::FunctionEnum::FITNESS,
		                           evolution_.verbose());

		if (!printFooter) return value; // <--- EARLY EXIT HERE

		std::cerr<<"QuantumOracle::minimize(): ";
		if (min.status() == MinimizerType::GSL_SUCCESS) {
			std::cerr<<" converged after ";
		} else {
			std::cerr<<"NOT CONVERGED after ";
		}

		++used;
		std::cerr<<used<<" iterations. "<<toString(angles)<<"\n";

		return value;
	}

	RealType maxFitness() const { return 100; }

	PsimagLite::String info(const ChromosomeType& chromosome) const
	{
		SizeType threadNum = 0;
		evolution_.setInput(0, fitParams_.inVector, threadNum);
		return fitParams_.hamiltonian.info(chromosome);
	}

private:

	static PsimagLite::String toString(const VectorRealType& angles)
	{
		const SizeType n = angles.size();
		PsimagLite::String str;
		for (SizeType i = 0; i < n; ++i)
			str += ttos(angles[i]) + " ";
		return str;
	}

	EvolutionType& evolution_;
	const GroundStateParamsType fitParams_;
}; // class QuantumOracle
} // namespace Gep

#endif // EVENDIM_GROUND_STATE_FITNESS_H
