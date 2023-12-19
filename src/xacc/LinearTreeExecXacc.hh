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
		SizeType numberOfBits;
		SizeType threadNum;
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

		ProgramType program = createProgram(circuit2);
		SizeType numberOfBits = log2Exact(initVector.size());
		return HandleType { program, numberOfBits, threadNum };
	}

	RealType energy(const HandleType& handle, const HamiltonianType& hamiltonian) const
	{
		if (handle.numberOfBits != hamiltonian.numberOfSites()) {
			throw std::runtime_error("Hamiltonian size incorrect\n");
		}

		auto buffer = xacc::qalloc(handle.numberOfBits);
		double angle = 0.;
		auto evaled = handle.program->operator()({ angle });
		auto accelerator = xacc::getAccelerator("tnqvm");

		// TODO: implement Hamiltonian::observe()
		auto rotatedCircuits = hamiltonian.observe(evaled);
		accelerator->execute(buffer, rotatedCircuits);

		// TODO: implement Hamiltonian::postProcess()
		auto energy = hamiltonian.postProcess(buffer);
		return energy;
	}

private:

	static ProgramType createProgram(const VecStringType& circuit)
	{
		// Get the IRProvider and create an
		// empty CompositeInstruction
		auto provider = xacc::getIRProvider("quantum");
		auto program = provider->createComposite("foo", { "t" });
		std::vector<InstructionType> instructions;
		SizeType ngates = circuit.size();
		for (SizeType i = 0; i < ngates; ++i) {
			QuantumGEPGate gate(circuit[i]);
			if (!gate.isParametric()) {
				auto someGate = provider->createInstruction(gate.name(), gate.bits());
				instructions.push_back(someGate);
			}
			else {
				auto someGate = provider->createInstruction(gate.name(), gate.bits(), gate.params());
				instructions.push_back(someGate);
			}
		}

		// Create X, Ry, CX, and Measure gates
		// auto x = provider->createInstruction("X", { 0 });
		// auto ry = provider->createInstruction("Ry", { 1 }, { "t" });
		// auto cx = provider->createInstruction("CNOT", { 1, 0 });
		auto m0 = provider->createInstruction("Measure", { 0 });
		instructions.push_back(m0);

		// Add them to the CompositeInstruction
		program->addInstructions(instructions);
		return program;
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
