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
#include "Engine.h"
#include "Evolution.h"
#include "Fitness/GroundStateFitness.h"
#include "Fitness/QuantumFitness.h"
#include "FloatingPoint.h"
#include "InputCheck.h"
#include "InputNg.h"
#include "Primitives/QuantumCircuit.h"
#include "Primitives/QuasiVector.hh"
#include "ProgramGlobals.h"
#include <unistd.h>

typedef double RealType;
typedef std::complex<RealType> ComplexType;
typedef Gep::QuasiVector<ComplexType> QuasiVectorType;
typedef Gep::QuantumCircuit<QuasiVectorType> PrimitivesType;
typedef typename PrimitivesType::CanonicalFormType CanonicalFormType;
typedef Gep::Evolution<PrimitivesType> EvolutionType;
typedef Gep::ParametersEngine<RealType> ParametersEngineType;
typedef PsimagLite::Tree<PrimitivesType> TreeType;
typedef Gep::Chromosome<TreeType, EvolutionType, ParametersEngineType> ChromosomeType;

template <template <typename> class FitnessTemplate>
RealType getFitness2(PsimagLite::InputNg<Gep::InputCheck>::Readable& io,
                     const Gep::ParametersEngine<double>& params,
                     EvolutionType& evolution,
                     const ChromosomeType& chromosome)
{
	typedef Gep::Engine<FitnessTemplate, EvolutionType> EngineType;
	typedef typename EngineType::FitnessType FitnessType;
	typedef typename FitnessType::FitnessParamsType FitnessParamsType;

	constexpr SizeType samples = 1;
	constexpr SizeType seed = 1234;
	constexpr SizeType threadNum = 0;
	FitnessParamsType fitParams(io, params.threads);
	FitnessType fitness(samples, evolution, &fitParams);

	return fitness.getFitness(chromosome, seed, threadNum);
}

RealType getFitness(PsimagLite::InputNg<Gep::InputCheck>::Readable& io,
                    const Gep::ParametersEngine<double>& params,
                    EvolutionType& evolution,
                    const ChromosomeType& chromosome)
{
	PsimagLite::String runType;
	io.readline(runType, "RunType=");

	if (runType == "FunctionFit") {
		return getFitness2<Gep::QuantumFitness>(io, params, evolution, chromosome);
	}
	else if (runType == "GroundState") {
		return getFitness2<Gep::GroundStateFitness>(io, params, evolution, chromosome);
	}

	throw PsimagLite::RuntimeError("RunType=FunctionFit or GroundState, but not " + runType + "\n");
}

int main(int argc, char* argv[])
{
	typedef typename ChromosomeType::VectorStringType VectorStringType;

	PsimagLite::String filename;
	PsimagLite::String vectorFilename;
	bool verbose = false;
	PsimagLite::FloatingPoint::enableExcept();
	SizeType randomSize = 0;
	PsimagLite::String splitString = " ";
	int opt = 0;

	PsimagLite::String strUsage(argv[0]);
	strUsage += " -f filename [-i filenameForVector] [-v] individual | -r size\n";
	strUsage += "\t-f filename similar to the one used by quantumGep driver\n";
	strUsage += "\t-i filenameForVector is an ASCII file with number of entries first " + PsimagLite::String("followed by entries separated by C++ whitespace\n");
	strUsage += "\tindividual is a comma-separated list of gates ending in 0\n";
	strUsage += "\t-r size will generate a random vector of norm 1\n";

	while ((opt = getopt(argc, argv, "f:i:r:c:p:v")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		case 'i':
			vectorFilename = optarg;
			break;
		case 'r':
			randomSize = PsimagLite::atoi(optarg);
			break;
		case 'v':
			verbose = true;
			break;
		case 's':
			splitString = optarg;
			break;
		case 'p':
			std::cout.precision(PsimagLite::atoi(optarg));
			break;
		default:
			throw PsimagLite::RuntimeError(strUsage);
			return 1;
		}
	}

	if (randomSize > 0) {
		QuasiVectorType rVector;
		PsimagLite::MersenneTwister rng(12345);
		rVector.randomize(randomSize, rng, 1., 0.);
		rVector.print(std::cout);
		return 0;
	}

	if (filename.empty())
		throw PsimagLite::RuntimeError(strUsage);

	Gep::InputCheck inputCheck;
	PsimagLite::InputNg<Gep::InputCheck>::Writeable input(filename, inputCheck);
	PsimagLite::InputNg<Gep::InputCheck>::Readable io(input);

	Gep::ParametersInput gepOptions(io);

	if (vectorFilename.empty()) {
		io.readline(vectorFilename, "InVectorFile=");
	}

	SizeType numberOfBits = 0;
	io.readline(numberOfBits, "NumberOfBits=");

	if (numberOfBits == 0)
		err("You need numberOfBits > 0\n");

	PsimagLite::String gates;
	if (gepOptions.primitives == "" || gepOptions.primitives == "?")
		gates = "C,H,P";
	else
		gates = gepOptions.primitives;

	if (gepOptions.primitives == "?") {
		std::cout << "Default gates are: " << gates << "\n";
		return 0;
	}

	SizeType seed = 12345;
	try {
		io.readline(seed, "RngSeed=");
	}
	catch (std::exception&) {
	}

	VectorStringType tokens;

	if (argc < 2)
		throw PsimagLite::RuntimeError(strUsage);

	PsimagLite::split(tokens, argv[argc - 1], splitString);

	if (tokens.size() == 0)
		err("No individual specified\n");

	SizeType last = tokens.size() - 1;
	if (tokens[last] != "0")
		err("Last node of quantum individual must be zero\n");

	gepOptions.head = tokens.size() - 1;

	ParametersEngineType params(gepOptions);
	PrimitivesType primitives(numberOfBits, gates, io);
	EvolutionType evolution(primitives, seed, verbose);

	constexpr SizeType threadNum = 0;
	ChromosomeType chromosome(params,
	                          evolution,
	                          tokens,
	                          threadNum);

	QuasiVectorType inVector(vectorFilename);

	const SizeType x = (1 << numberOfBits);
	if (x != inVector.size())
		err("File " + vectorFilename + " should contain " + ttos(x) + " entries.\n");

	std::cout << "Norm of input state= " << inVector.norm() << "\n";
	evolution.nodeHelper().setInput(0, inVector, threadNum);

	QuasiVectorType outVector = chromosome.exec(0);
	std::cout << "Norm of output state= " << outVector.norm() << "\n";

	outVector.print(std::cout);

	RealType f = getFitness(io, params, evolution, chromosome);
	std::cout << "Fitness= " << f << "\n";

	VectorStringType vecStr = chromosome.effectiveVecString();
	CanonicalFormType canonicalForm(vecStr, evolution.nodeHelper().nodeFactory());
	canonicalForm.changeIfNeeded(vecStr);
	for (SizeType i = 0; i < vecStr.size(); ++i)
		std::cout << vecStr[i] << " ";
	std::cout << "\n";
}
