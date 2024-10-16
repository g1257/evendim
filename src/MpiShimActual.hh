#ifndef EVENDIM_MPI_SHIM_ACTUAL_HH
#define EVENDIM_MPI_SHIM_ACTUAL_HH
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

private:

	int rank_;
	int size_;
};

}

#endif
