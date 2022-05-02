# Preliminaries
## Disclaimer and Licensing

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
The full software license for evendim version 1.0.0
can be found in
file LICENSE.

## Please cite this work

evendim is a free and open source genetic expression programming code. 
The full software license for evendim version 0.
can be found in
file LICENSE.
You are welcomed to use it and publish data
obtained with evendim. If you do, please cite this
work. Explain How To Cite This Work. FIXME. TBW.

## References

	@book{ferreira2006gene,
	title={Gene Expression Programming: Mathematical Modeling by an Artificial Intelligence},
	author={Ferreira, C.},
	isbn={9783540328490},
	series={Studies in Computational Intelligence},
	url={https://books.google.com/books?id=NkG7BQAAQBAJ},
	year={2006},
	publisher={Springer Berlin Heidelberg}
	}

## Hash of the latest commit

Hash of the latest commit is also posted at
FIXME

# Building and Running evendim

## Required Software

* GNU C++
* The LAPACK and BLAS libraries
* The GSL library
* PsimagLite (see below)

## Optional Software

* make or gmake (only needed to use the Makefile)
* perl (may be needed to run some auxiliary script)

## Quick Start

1. Use your distribution repository tool to install gcc with support for C++,
the LAPACK and BLAS libraries, the gsl library, make, perl, and git
if you don't have them.

2. Issue

    cd someDirectory/

    git clone https://github.com/g1257/PsimagLite.git

    git clone https://github.com/g1257/evendim.git

3. Compile PsimagLite

4. Now issue

    cd evendim/src

    cp Config.make.sample Config.make

    make

5. You can run the arithmetic primitives with
    ./gep2 -i 1 -h 5 -p 100 -t 10
	which will run an arithmetic test with one
	input, head size of a maximum of 5, population
	100 individuals, and for t=10 generations.
	The function is in src/Functions/Example1.h
	and is f(x) =  x * (x - 1) * (x + 1).
