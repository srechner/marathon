## marathon

This C++ library is designed to support the analysis of Markov chain based sampling methods. It provides functions for the analysis of so-called state graphs. For an introduction into its functionality, see the article at http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0147935.

## Current Development Status:
* Working on Documentation
* Integrating of network analysis
* Integrating additional chains

## Requirements:

This library is developed and tested at Linux systems (primarily Ubuntu). In principle it should possible to migrate the library to other operating systems.

Compiling the library requires the following software packages:
 * `g++` in version 4.8 or later
 * `boost` headers 
 * `cblas` headers
 * `cuda` in version 7.0 or later (optional)

## Installation:

First Step: Build the `marathon` library.

1. Download the sources as zip archive.
2. Extract the archive to a directory of your choice.
3. Enter the `marathon` directory.
4. Run `make` to build `marathon` in standard C++ mode (without CUDA support). To compile with CUDA support, run instead `make CUDA=true`.

Second Step: Build an application and link to the `marathon` library.

The [example directory](https://github.com/srechner/marathon/blob/master/examples/) contains several example programs that demonstrate how `marathon` applications can look like. To compile an example program:

1. Enter an example directory.
2. Run `make` or `make CUDA=true`. The command has to match the mode you have used to build the library. Please see next section if something goes wrong.
3. Tell your system where to find the `marathon` library: `export LD_LIBRARY_PATH=../:$LD_LIBRARY_PATH`
4. Run the application.

### Link Options:

Linking your application requires several third-party libraries, which have to exist at your system. (Names are those of Ubuntu packages.)
 * `gomp`
 * `pthread`
 * `openblas` (or an equivalent `CBLAS` implementation)
 * `arpack++`
 * `arpack`
 * `superlu`
 * `cublas`	(only when built with CUDA support)

## Example:

This instruction shows how to build the [MixingBounds](https://github.com/srechner/marathon/blob/master/examples/MixingBounds/) example.

First Step: Build the library.

	git clone https://github.com/srechner/marathon.git
	cd marathon
	make

Second Step: Compile the application.

	cd examples/MixingBounds/
	make
	export LD_LIBRARY_PATH=../:$LD_LIBRARY_PATH
	./MixingBounds swapBip "2,1,1;1,2,1" 1e-3

The output should look like:

	number of states:          5
	number of transition arcs: 21
	total mixing time:         72
	lower spectral bound:      34.1803
	upper spectral bound:      102.206
	upper congestion bound:    183.971

