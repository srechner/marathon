## marathon

This C++ library is designed to support the analysis of Markov chain based sampling methods. It provides functions for the analysis of so-called state graphs. For an introduction into its functionality, see my article at arXiv http://arxiv.org/abs/1508.04740.

Current ToDo's:

	Working on Documentation
	
## Requirements:

	g++ 4.8 or later
	cuda 7.0 or later (optional)
	boost headers (libboost-dev in ubuntu)
	cblas headers (libblas-dev in ubuntu)

## Installation:

It is possible to build the library in C++ native mode without CUDA support.
To do so, just run `make cpp`. To compile with cuda support, run `make cuda`.

## Link against marathon:

To use the functionality build in marathon, link against the libmarathon.so file.
You will need the following other libraries:

	gomp
	pthread
	openblas
	arpack++
	arpack
	superlu
	cublas	(only when built with cuda support)

## Example:

The example directory contains an [example](https://github.com/srechner/marathon/blob/master/examples/totalMixingTime.cpp) program which shows how to use marathon to compute the total mixing time. To install, first compile the library via `make cpp` or `make cuda`. To run the example, do the following:
	
	cd examples
	make
	export LD_LIBRARY_PATH=../:$LD_LIBRARY_PATH
	./totalMixingTime
