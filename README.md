## marathon

This C++ library is designed to support the analysis of Markov chain based sampling methods. It provides functions for the analysis of so-called state graphs. For an introduction into its functionality, see the article at http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0147935.

Current ToDo's:

	Working on Documentation
	Major Refactoring and Redesign due to big expansion plans
	Integration of network analysis
	Smaller Modifications to avoid bad memory allocations
	
## Requirements:

	g++ 4.8 or later
	cuda 7.0 or later (optional)
	boost headers (libboost-dev in ubuntu)
	cblas headers (libblas-dev in ubuntu)

## Installation:

It is possible to build the library in C++ native mode without CUDA support.
To do so, just run `make`. To compile with cuda support, run `make GPU=true`.

## Link against marathon:

To use the functionality build in marathon, link against the libmarathon.so file.
You will need the following other libraries installed:

	gomp
	pthread
	openblas
	arpack++
	arpack
	superlu
	cublas	(only when built with cuda support)

## Example:

The [example directory](https://github.com/srechner/marathon/blob/master/examples/) contains several example programs which show how to use marathon. To install, first compile the library via `make` or `make GPU=true`. To run the example, do the following. (The make option has to match the option when building the library.)
	
	cd examples/MixingBounds/
	make
	export LD_LIBRARY_PATH=../:$LD_LIBRARY_PATH
	./MixingBounds swapBip "2,1,1;1,2,1" 1e-3
