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
#include "BaseFitness.h"
#include "GroundStateParams.h"
#include "LinearTreeExec.hh"
#include "MersenneTwister.h"
#include "Minimizer.h"
#include "MinimizerParams.h"
#include "PsimagLite.h"

namespace Gep {

template <typename ChromosomeType, typename EvolutionType, typename GroundStateParamsType>
class FunctionToMinimize2 {
public:

	typedef typename GroundStateParamsType::ComplexType ComplexType;
	typedef typename EvolutionType::PrimitivesType PrimitivesType;
	typedef typename PrimitivesType::NodeType NodeType;
	typedef typename NodeType::ValueType QuasiVectorType;
	typedef typename PsimagLite::Real<ComplexType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef RealType FieldType;
	typedef typename ChromosomeType::VectorStringType VectorStringType;
	typedef PsimagLite::Matrix<ComplexType> MatrixType;
	typedef typename EvolutionType::NodeHelperType::NodeFactoryType NodeFactoryType;
	using LinearTreeExecType = LinearTreeExec<typename NodeType::ValueType, typename NodeType::AnglesType, NodeFactoryType>;

	enum class FunctionEnum { FITNESS,
		                  DIFFERENCE };

	FunctionToMinimize2(EvolutionType& evolution,
	                    const ChromosomeType& chromosome,
	                    const GroundStateParamsType& groundStateParams,
	                    SizeType thread)
	    : evolution_(evolution)
	    , chromosome_(chromosome)
	    , groundStateParams_(groundStateParams)
	    , outVector_(groundStateParams_.inVector.size())
	    , threadNum_(thread)
	{
		outVector_.blowUp(outVector_.size());
		numberOfAngles_ = findNumberOfAngles(chromosome.effectiveVecString());
	}

	SizeType size() const { return numberOfAngles_; }

	RealType operator()(VectorRealType& angles)
	{
		// the case without angles is handled elsewhere
		// if it reaches here, it's an internal error
		assert(numberOfAngles_ > 0);

		assert(angles.size() == numberOfAngles_);

		// false means don't be verbose here
		return fitness(&angles, FunctionEnum::DIFFERENCE, false, false);
	}

	void df(VectorRealType& dest, const VectorRealType& angles)
	{
		const SizeType geneLength = chromosome_.geneLength();
		VectorStringType vecStr = chromosome_.vecString();
		const SizeType numberOfGenes = chromosome_.size();
		encodeAngles(vecStr, angles, geneLength, numberOfGenes);
		const ChromosomeType* chromosome = new ChromosomeType(chromosome_.params(),
		                                                      evolution_,
		                                                      vecStr,
		                                                      threadNum_);

		dest.resize(angles.size());
		evolution_.nodeHelper().setInput(0, groundStateParams_.inVector, threadNum_);

		const QuasiVectorType& inVector = groundStateParams_.inVector;
		for (SizeType angleIndex = 0; angleIndex < numberOfAngles_; ++angleIndex) {
			evolution_.nodeHelper().setInput(0, inVector, threadNum_);

			computeDifferentialVector(differential_, angles, angleIndex);

			const RealType tmp = diffVectorDiff2(chromosome->exec(0),
			                                     outVector_,
			                                     differential_);
			dest[angleIndex] += tmp;
		}
	}

	RealType fitness(VectorRealType* angles,
	                 FunctionEnum functionEnum,
	                 bool useXaccOptimizer,
	                 bool verbose)
	{
		const ChromosomeType* chromosome = nullptr;
		VectorStringType vecStr = chromosome_.vecString();

		if (angles) {
			const SizeType geneLength = chromosome_.geneLength();
			const SizeType numberOfGenes = chromosome_.size();
			encodeAngles(vecStr, *angles, geneLength, numberOfGenes);
			chromosome = new ChromosomeType(chromosome_.params(),
			                                evolution_,
			                                vecStr,
			                                threadNum_);
		}
		else {
			chromosome = &chromosome_;
		}

		evolution_.nodeHelper().setInput(0, groundStateParams_.inVector, threadNum_);
		if (verbose)
			evolution_.nodeHelper().printInputs(std::cout);

		// oracle goes here
		RealType e = 0;
		if (chromosome->params().options.isSet("useLinearTreeIfPossible")
		    && chromosome->isLinearTree()) {
			const LinearTreeExecType& linearTreeExec = evolution_.nodeHelper().linearTreeExec();
			typename LinearTreeExecType::HandleType handle = linearTreeExec.getHandle(groundStateParams_.inVector, chromosome->vecString(), threadNum_);
			e = linearTreeExec.energy(handle, groundStateParams_.hamiltonian, groundStateParams_.useXaccOptimizer);
			if (useXaccOptimizer && angles) {
				linearTreeExec.fillAngles(*angles, handle);
			}
		}
		else {
			e = groundStateParams_.hamiltonian.energy(chromosome->exec(0), threadNum_);
		}

		if (angles) {
			delete chromosome;
			chromosome = nullptr;
		}

		return (functionEnum == FunctionEnum::DIFFERENCE) ? e : -e;
	}

	static void encodeAngles(VectorStringType& vecStr,
	                         const VectorRealType& angles,
	                         SizeType geneLength,
	                         SizeType numberOfGenes)
	{
		const SizeType n = vecStr.size();
		SizeType currentIndex = 0;
		SizeType genesSoFar = 0;
		bool flag = true;
		bool inCell = false;

		for (SizeType i = 0; i < n; ++i) {

			if (!inCell && numberOfGenes > 1 && i % geneLength == 0) {
				flag = true;

				if (genesSoFar == numberOfGenes) {
					inCell = true;
				}

				++genesSoFar;
			}

			if (EvolutionType::isAnInteger(vecStr[i])) {
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

	static void initAngles(VectorRealType& angles,
	                       const VectorStringType& vStr,
	                       long unsigned int seed)
	{
		const SizeType n = vStr.size();
		SizeType currentIndex = 0;
		PsimagLite::MersenneTwister rng(seed);
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

private:

	template <typename SomeRngType>
	static void initAngle(RealType& angle,
	                      PsimagLite::String str,
	                      SomeRngType& rng)
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

	void computeDifferentialVector(QuasiVectorType& differential,
	                               const VectorRealType& angles,
	                               SizeType angleIndex)
	{
		assert(angles.size() == numberOfAngles_);

		VectorStringType cString = chromosome_.effectiveVecString();

		VectorStringType tmpString = replaceOneR(cString, angleIndex, chromosome_.vecString(0).size());

		assert(tmpString.size() == chromosome_.vecString(0).size());

		// create derivative individual angle-th
		ChromosomeType newChromosome(chromosome_.params(), evolution_, tmpString, threadNum_);

		// apply to inVector
		evolution_.nodeHelper().setInput(0, groundStateParams_.inVector, threadNum_);
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
			if (m == 0)
				continue;
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
	QuasiVectorType outVector_;
	QuasiVectorType differential_;
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
template <typename ChromosomeType>
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

	using QuasiVectorType = QuasiVector<ComplexType>;
	using NodeType = typename PrimitivesType::NodeType;
	using NodeFactoryType = typename EvolutionType::NodeHelperType::NodeFactoryType;
	using LinearTreeExecType = LinearTreeExec<typename NodeType::ValueType, typename NodeType::AnglesType, NodeFactoryType>;

	enum { HAS_XACC = LinearTreeExecType::HAS_XACC };

	using GroundStateParamsType = GroundStateParams<HamiltonianType, ComplexType, HAS_XACC>;
	using MinimizerParamsType = typename GroundStateParamsType::MinimizerParamsType;
	using FunctionToMinimizeType = FunctionToMinimize2<ChromosomeType, EvolutionType, GroundStateParamsType>;
	using FitnessParamsType = GroundStateParamsType;

	GroundStateFitness(SizeType samples,
	                   EvolutionType& evolution,
	                   FitnessParamsType* fitParams)
	    : evolution_(evolution)
	    , fitParams_(*fitParams)
	{
		if (evolution.nodeHelper().numberOfInputs() != 1)
			err("QuantumOracle::ctor(): 1 input expected\n");
		if (samples != 1)
			err("Expecting samples == 1\n");

		const SizeType threadNum = 0;
		evolution.nodeHelper().setInput(0, fitParams->inVector, threadNum);
	}

	RealType getFitness(const ChromosomeType& chromosome,
	                    long unsigned int seed,
	                    SizeType threadNum)
	{
		evolution_.nodeHelper().setInput(0, fitParams_.inVector, threadNum);

		RealType norma = fitParams_.inVector.norm();
		if (fabs(norma - 1) > 1e-4)
			err("Input vector not normalized\n");

		FunctionToMinimizeType f(evolution_, chromosome, fitParams_, threadNum);

		if (f.size() == 0) {
			// if no angles
			return f.fitness(nullptr,
			                 FunctionToMinimizeType::FunctionEnum::FITNESS,
			                 fitParams_.useXaccOptimizer,
			                 evolution_.verbose());
		}
		else {
			// there are angles so we can ...
			VectorRealType angles(f.size());
			FunctionToMinimizeType::initAngles(angles, chromosome.effectiveVecString(), seed);

			// ... use XACC optimizer
			if (fitParams_.useXaccOptimizer) {

				double value = f.fitness(&angles,
				                         FunctionToMinimizeType::FunctionEnum::FITNESS,
				                         fitParams_.useXaccOptimizer,
				                         evolution_.verbose());
				ChromosomeType* chromosomeNonconst = const_cast<ChromosomeType*>(&chromosome);
				setAnglesIntoChromosome(*chromosomeNonconst, angles, threadNum);
				return value;
			}
			else { // ... or use QuantumGEP's internal optimizer
				return optimizeAndReturnFitness(angles, f, chromosome, threadNum);
			}
		}
	}

	RealType optimizeAndReturnFitness(VectorRealType& angles,
	                                  FunctionToMinimizeType& f,
	                                  const ChromosomeType& chromosome,
	                                  SizeType threadNum)
	{
		using MinimizerType = typename PsimagLite::Minimizer<RealType, FunctionToMinimizeType>;

		const MinimizerParamsType& minParams = fitParams_.minParams;
		MinimizerType min(f, minParams.maxIter, minParams.verbose);

		int used = 0;
		if (minParams.algo == MinimizerParamsType::SIMPLEX) {
			used = min.simplex(angles,
			                   minParams.delta,
			                   minParams.tol);
		}
		else if (minParams.algo == MinimizerParamsType::NONE) {
			used = 1;
		}
		else {
			used = min.conjugateGradient(angles,
			                             minParams.delta,
			                             minParams.delta2,
			                             minParams.tol,
			                             minParams.saveEvery);
		}

		int status = (min.status() == MinimizerType::GSL_SUCCESS) ? 0 : 1;

		// set angles if successful
		if (status == 0) {
			ChromosomeType* chromosomeNonconst = const_cast<ChromosomeType*>(&chromosome);
			setAnglesIntoChromosome(*chromosomeNonconst, angles, threadNum);
		}

		const bool printFooter = minParams.verbose;
		// const int returnStatus = (used > 0) ? 0 : 1;

		assert(!fitParams_.useXaccOptimizer);
		RealType value = f.fitness(&angles,
		                           FunctionToMinimizeType::FunctionEnum::FITNESS,
		                           fitParams_.useXaccOptimizer,
		                           evolution_.verbose());

		if (!printFooter)
			return value; // <--- EARLY EXIT HERE

		std::cerr << "QuantumOracle::minimize(): ";
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

	RealType maxFitness() const { return 100; }

	PsimagLite::String info(const ChromosomeType& chromosome) const
	{
		SizeType threadNum = 0;
		evolution_.nodeHelper().setInput(0, fitParams_.inVector, threadNum);
		return infoInternal(chromosome);
	}

private:

	void setAnglesIntoChromosome(ChromosomeType& chromosome, const VectorRealType& angles, SizeType threadNum)

	{
		using VectorStringType = typename ChromosomeType::VectorStringType;

		VectorStringType vecStr = chromosome.vecString();
		const SizeType numberOfGenes = chromosome.size();
		const SizeType geneLength = chromosome.geneLength();
		FunctionToMinimizeType::encodeAngles(vecStr, angles, geneLength, numberOfGenes);
		const ChromosomeType* chromosome2 = new ChromosomeType(chromosome.params(),
		                                                       evolution_,
		                                                       vecStr,
		                                                       threadNum);

		chromosome = *chromosome2;

		delete chromosome2;
		chromosome2 = nullptr;
	}

	template <typename SomeChromosomeType>
	static PsimagLite::String infoInternal(const SomeChromosomeType& chromosome)
	{
		return info_(chromosome.exec(0), 1e-4);
	}

	static PsimagLite::String info_(const QuasiVectorType& v, double epsilon)
	{
		const SizeType n = v.size();
		PsimagLite::String buffer;
		for (SizeType i = 0; i < n; ++i) {
			if (v.hasWeight(i, epsilon))
				buffer += ttos(i) + " ";
		}

		return buffer;
	}
	static PsimagLite::String toString(const VectorRealType& angles)
	{
		const SizeType n = angles.size();
		PsimagLite::String str;
		for (SizeType i = 0; i < n; ++i)
			str += ttos(angles[i]) + " ";
		return str;
	}

	EvolutionType& evolution_;
	const GroundStateParamsType& fitParams_;
}; // class QuantumOracle
} // namespace Gep

#endif // EVENDIM_GROUND_STATE_FITNESS_H
