## marathon library

version 0.1

## Requirements:

	g++ 4.8 or later
	cuda 7.0 or later (optional)
	boost headers (libboost-dev)
	cblas headers (libblas-dev)

## Installation:

It is possible to build the library in cpu native mode without cuda support.
To do so, just run 'make cpp'. To compile with cuda support, run 'make cuda'.

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
