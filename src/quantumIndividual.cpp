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
#include "Evolution.h"
#include "Primitives/QuantumCircuit.h"
#include "Engine.h"
#include <unistd.h>
#include "Functions/QuantumOracle.h"
#include "InputNg.h"
#include "InputCheck.h"
#include "FloatingPoint.h"
#include "ProgramGlobals.h"

template<typename SomeType, typename SomeRngType>
void randomVector(std::vector<SomeType>& outVector, SomeRngType& rng)
{
	typedef typename PsimagLite::Real<SomeType>::Type RealType;

	const SizeType n = outVector.size();
	RealType sum = 0;
	for (SizeType i = 0; i < n; ++i) {
		SomeType value = rng();
		outVector[i] = value;
		sum += PsimagLite::real(PsimagLite::conj(value)*value);
	}

	assert(sum > 0);
	RealType factor = 1/sqrt(sum);
	for (SizeType i = 0; i < n; ++i)
		outVector[i] *= factor;

}

template<typename SomeType>
void writeVector(std::ostream& os, const std::vector<SomeType>& outVector)
{
	const SizeType n = outVector.size();
	os<<n<<"\n";
	for (SizeType i = 0; i < n; ++i)
		os<<outVector[i]<<" ";
	os<<"\n";
}

int main(int argc, char* argv[])
{
	typedef double RealType;
	typedef std::complex<RealType> ComplexType;
	typedef PsimagLite::Vector<ComplexType>::Type VectorType;
	typedef Gep::QuantumCircuit<VectorType> PrimitivesType;
	typedef Gep::Evolution<PrimitivesType> EvolutionType;
	typedef Gep::ParametersEngine<RealType> ParametersEngineType;
	typedef Gep::Tree<PrimitivesType> TreeType;
	typedef Gep::Chromosome<TreeType, EvolutionType, ParametersEngineType> ChromosomeType;
	typedef typename ChromosomeType::VectorStringType VectorStringType;

	PsimagLite::String filename;
	PsimagLite::String vectorFilename;
	bool verbose = false;
	PsimagLite::FloatingPoint::enableExcept();
	SizeType randomSize = 0;
	int opt = 0;

	PsimagLite::String strUsage(argv[0]);
	strUsage += " -f filename -i filenameForVector [-v] individual | -r size\n";
	strUsage += "\t-f filename similar to the one used by quantumGep driver\n";
	strUsage += "\t-i filenameForVector is an ASCII file with number of entries first " +
	        PsimagLite::String("followed by entries separated by C++ whitespace\n");
	strUsage += "\tindividual is a comma-separated list of gates ending in 0\n";
	strUsage += "\t-r size will generate a random vector of norm 1\n";

	while ((opt = getopt(argc, argv,"f:i:r:v")) != -1) {
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
		default:
			throw PsimagLite::RuntimeError(strUsage);
			return 1;
		}
	}

	if (randomSize > 0) {
		VectorType rVector(randomSize);
		PsimagLite::MersenneTwister rng(12345);
		randomVector(rVector, rng);
		writeVector(std::cout, rVector);
		return 0;
	}

	if (filename == "" || vectorFilename == "")
		throw PsimagLite::RuntimeError(strUsage);

	Gep::InputCheck inputCheck;
	PsimagLite::InputNg<Gep::InputCheck>::Writeable input(filename, inputCheck);
	PsimagLite::InputNg<Gep::InputCheck>::Readable io(input);

	Gep::ParametersInput gepOptions(io);

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
		std::cout<<"Default gates are: "<<gates<<"\n";
		return 0;
	}

	SizeType seed = 12345;
	try {
		io.readline(seed, "RngSeed=");
	} catch (std::exception&) {}

	VectorStringType tokens;

	if (argc < 2)
		throw PsimagLite::RuntimeError(strUsage);

	PsimagLite::split(tokens, argv[argc - 1], ",");

	if (tokens.size() == 0)
		err("No individual specified\n");

	SizeType last = tokens.size() - 1;
	if (tokens[last] != "0")
		err("Last node of quantum individual must be zero\n");

	gepOptions.head = tokens.size() - 1;

	ParametersEngineType params(gepOptions);
	PrimitivesType primitives(numberOfBits, gates, 1); // threads == 1
	EvolutionType evolution(primitives, seed, verbose);

	ChromosomeType  chromosome(params,
	                           evolution,
	                           tokens);

	VectorType inVector;
	Gep::ProgramGlobals::readVector(inVector, vectorFilename);
	const SizeType x = (1 << numberOfBits);
	if (x != inVector.size())
		err("File " + vectorFilename + " should contain " + ttos(x) + " entries.\n");

	constexpr SizeType threadId = 0;
	evolution.setInput(0, inVector, threadId);

	VectorType outVector = chromosome.exec(0);

	writeVector(std::cout, outVector);
}
