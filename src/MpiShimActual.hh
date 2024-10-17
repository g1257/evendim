#ifndef EVENDIM_MPI_SHIM_ACTUAL_HH
#define EVENDIM_MPI_SHIM_ACTUAL_HH
#include "PsimagLite.h"
#include <mpi.h>

namespace Gep {

class MpiShim {

public:

	MpiShim(int argc, char** argv)
	{
		// Initialize MPI
		MPI_Init(&argc, &argv);

		// Get the rank of the process
		MPI_Comm_rank(MPI_COMM_WORLD, &rank_);

		// Get the total number of processes
		MPI_Comm_size(MPI_COMM_WORLD, &size_);
	}

	~MpiShim()
	{
		// Finalize MPI
		MPI_Finalize();
	}

	unsigned int rank() const { return rank_; }

	unsigned int size() const { return size_; }

	bool isMPI() const { return true; }

	std::string buildInput(const std::string& name)
	{
		bool is_open = false;
		std::string buffer;
		std::string var;

		for (char c : name) {
			if (c == ';' && is_open) {
				is_open = false;
				buffer += getVarValue(var);
				var = "";
				continue;
			}

			if (c == '%' && !is_open) {
				is_open = true;
				var = "";
				continue;
			}

			if (is_open) {
				var += c;
				continue;
			}

			buffer += c;
		}

		return buffer;
	}

	PsimagLite::String buildOutput(PsimagLite::String filename)
	{
		PsimagLite::String rootname = PsimagLite::basename(filename);
		size_t index = rootname.find(".", 0);
		if (index != PsimagLite::String::npos) {
			rootname.erase(index, filename.length());
		}

		return "runFor" + rootname + ".cout";
	}

private:

	std::string getVarValue(const std::string& name)
	{
		if (name == "i") {
			return ttos(rank_);
		}
		else if (name == "n") {
			return ttos(size_);
		}
		else {
			return "";
		}
	}

	int rank_;
	int size_;
};

}

#endif
