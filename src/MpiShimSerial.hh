#ifndef EVENDIM_MPI_SHIM_SERIAL_HH
#define EVENDIM_MPI_SHIM_SERIAL_HH

namespace Gep {

class MpiShim {

public:

	MpiShim(int argc, char** argv)
	{
	}

	unsigned int rank() const { return 0; }

	unsigned int size() const { return 1; }

	bool isMPI() const { return false; }
};

}

#endif
