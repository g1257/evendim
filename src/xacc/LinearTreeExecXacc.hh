#ifndef EVENDIM_QUANTUMGEPXACC_H_
#define EVENDIM_QUANTUMGEPXACC_H_

#include "../Engine/UnderlyingType.hh"
#include "AllocatorCpu.h"
#include "Complex.h"
#include "HamiltonianXacc.hh"
#include "QuantumGEPGate.hh"
#include "xacc.hpp"
#include <string>
#include <vector>

namespace Gep {

template <typename T1, typename T2>
struct TypesEqual {
	static const bool value = false;
};

template <typename T>
struct TypesEqual<T, T> {
	static const bool value = true;
};

template <bool b, typename T1, typename T2>
struct FirstOrSecondType {
	using type = T2;
};

template <typename T1, typename T2>
struct FirstOrSecondType<true, T1, T2> {
	using type = T1;
};

template <typename VecComplexType, typename AnglesType, typename CtorParamType>
class LinearTreeExec {

public:

	using ComplexType = typename UnderlyingType<VecComplexType>::Type;
	using RealType = typename PsimagLite::Real<ComplexType>::Type;
	using VecStringType = std::vector<std::string>;
	using ProgramType = std::shared_ptr<xacc::CompositeInstruction>;
	using InstructionType = std::shared_ptr<xacc::Instruction>;
	using BogusFirstType = typename FirstOrSecondType<TypesEqual<VecComplexType, std::vector<ComplexType>>::value, int*, double*>::type;
	using BogusSecondType = typename FirstOrSecondType<!TypesEqual<VecComplexType, std::vector<ComplexType>>::value, int*, double*>::type;
	using HamiltonianType = Hamiltonian<ComplexType>;

	struct HandleType {
		ProgramType program;
		SizeType numberOfBits = 0;
		SizeType threadNum = 0;
		SizeType number_of_params = 0;
	};

	// Ctor not needed in the xacc version of LinearTreeExec
	LinearTreeExec(const CtorParamType& ctorParam)
	{
	}

	HandleType getHandle(const VecComplexType& initVector,
	                     const VecStringType& circuit,
	                     SizeType threadNum) const
	{
		VecStringType circuit2;
		pureVectorToXgates(circuit2, initVector);

		circuit2.insert(circuit2.end(), circuit.begin(), circuit.end());

		std::pair<ProgramType, SizeType> programAndNparams = createProgram(circuit2);
		SizeType numberOfBits = log2Exact(initVector.size());
		return HandleType { programAndNparams.first, numberOfBits, threadNum, programAndNparams.second };
	}

	RealType energy(const HandleType& handle, const HamiltonianType& hamiltonian) const
	{
		if (handle.numberOfBits != hamiltonian.numberOfSites()) {
			throw std::runtime_error("Hamiltonian size incorrect\n");
		}

		auto buffer = xacc::qalloc(handle.numberOfBits);
		std::vector<double> vector_of_params;

		if (handle.number_of_params > 0) {
			// set all angles to zero
			vector_of_params.resize(handle.number_of_params, 0.0);
			std::cerr << "NumberOfParams= " << handle.number_of_params << "\n";
		}

		auto evaled = handle.program->operator()(vector_of_params);
		auto accelerator = xacc::getAccelerator("tnqvm");

		auto rotatedCircuits = hamiltonian.observe(evaled);
		accelerator->execute(buffer, rotatedCircuits);

		auto energy = hamiltonian.postProcess(buffer);
		return energy;
	}

private:

	static std::pair<ProgramType, SizeType> createProgram(const VecStringType& circuit)
	{
		// Get the IRProvider and create an
		// empty CompositeInstruction
		auto provider = xacc::getIRProvider("quantum");
		std::vector<InstructionType> instructions;
		SizeType ngates = circuit.size();
		SizeType param_counter = 0;
		std::vector<std::string> total_params;
		for (SizeType i = 0; i < ngates; ++i) {
			QuantumGEPGate gate(circuit[i]);
			if (!gate.isParametric()) {
				auto someGate = provider->createInstruction(gate.name(), gate.bits());
				instructions.push_back(someGate);
			}
			else {
				SizeType nparams = gate.numberOfParams();
				if (nparams != 1) {
					throw std::runtime_error(std::string(__FILE__) + " I can only deal with 1 param for now\n");
				}

				std::cerr << "parametric " << gate.name() << "\n";
				std::string params = "t" + ttos(param_counter++);
				total_params.push_back(params);
				auto someGate = provider->createInstruction(gate.name(), gate.bits(), { params });
				instructions.push_back(someGate);
			}
		}

		// Create X, Ry, CX, and Measure gates
		// auto x = provider->createInstruction("X", { 0 });
		// auto ry = provider->createInstruction("Ry", { 1 }, { "t" });
		// auto cx = provider->createInstruction("CNOT", { 1, 0 });
		// auto m0 = provider->createInstruction("Measure", { 0 });
		// instructions.push_back(m0);

		// create program
		auto program = provider->createComposite("foo", total_params);

		// Add them to the CompositeInstruction
		program->addInstructions(instructions);
		return std::pair<ProgramType, SizeType>(program, param_counter);
	}

	// Avoid overload if second function exists

	static SizeType findPureState(const VecComplexType& initVector, BogusFirstType = 0)
	{
		return findPureState(initVector.toVector());
	}

	static SizeType findPureState(const std::vector<ComplexType>& initVector, BogusSecondType = 0)
	{
		SizeType n = initVector.size();
		bool hasSeenNonZero = false;
		SizeType x = 0;
		for (SizeType i = 0; i < n; ++i) {
			if (std::norm(initVector[i]) > 0.0) {
				if (hasSeenNonZero) {
					dieVectorNotPure(initVector, "1");
				}

				hasSeenNonZero = true;
				if (std::imag(initVector[i]) != 0) {
					dieVectorNotPure(initVector, "2");
				}

				if (std::abs(std::real(initVector[i]) - 1) > 1e-4) {
					dieVectorNotPure(initVector, "3");
				}

				x = i;
			}
		}

		if (!hasSeenNonZero) {
			throw std::runtime_error("initVector is zero\n");
		}

		return x;
	}

	static void pureVectorToXgates(VecStringType& circuit,
	                               const VecComplexType& initVector)

	{
		typename FirstOrSecondType<!TypesEqual<VecComplexType, std::vector<ComplexType>>::value, BogusFirstType, BogusSecondType>::type bogus = 0;
		SizeType x = findPureState(initVector, bogus);
		SizeType i = 0;
		while (x > 0) {
			if (x & 1) {
				circuit.push_back("Sx" + ttos(i));
			}

			x >>= 1;
			++i;
		}
	}

	// If n is 2^x, this function returns x
	// Else it throws
	static SizeType log2Exact(SizeType nn)
	{
		SizeType n = nn;
		SizeType x = 0;
		while (n > 0) {
			if (n & 1) {
				break;
			}

			n >>= 1;
			++x;
		}

		SizeType mustBeN = (1 << x);
		if (mustBeN != nn) {
			throw std::runtime_error("n is not a power of 2\n");
		}

		return x;
	}

	static void dieVectorNotPure(const std::vector<ComplexType>&, const std::string& msg)
	{
		throw std::runtime_error("initVector must be pure " + msg + "\n");
	}
};
}
#endif // LINEARTREEEXEC_DUMMY_HH
