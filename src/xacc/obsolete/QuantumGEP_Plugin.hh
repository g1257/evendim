#ifndef QUANTUM_GEP_PLUGIN_HH_
#define QUANTUM_GEP_PLUGIN_HH_

#include "Algorithm.hpp"
#include "AlgorithmGradientStrategy.hpp"
#include "xacc.hpp"
#include "xacc_service.hpp"
#include <chrono>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

namespace xacc {

namespace algorithm {

	class QuantumGEP : public Algorithm {

	protected:

		Optimizer* optimizer;
		Accelerator* accelerator;
		HeterogeneousMap parameters;
		std::shared_ptr<AlgorithmGradientStrategy> gradientStrategy;

	public:

		bool initialize(const HeterogeneousMap& parameters) override;
		const std::vector<std::string> requiredParameters() const override;
		void execute(const std::shared_ptr<AcceleratorBuffer> buffer) const override;
		const std::string name() const override { return "mc-vqe"; }
		const std::string description() const override { return ""; }
		DEFINE_ALGORITHM_CLONE(QuantumGEP)
	};
} // namespace algorithm
} // namespace xacc
#endif
