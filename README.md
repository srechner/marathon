# marathon 1.0

This C++ library is designed to support the analysis of Markov chain Monte Carlo methods. 
It provides functions for the construction and analysis of so-called state graphs. 
For an introduction into its functionality, see its introductional 
[article](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0147935).
If you want to use this library for your work, please feel free to contact me and cite

Steffen Rechner and Annabell Berger. *marathon: An open source software library for the
analysis of Markov-Chain Monte Carlo algorithms*. PLOS ONE **11** (2016). DOI: 10.1371/journal.pone.0147935.

## Main Features ##

* Construction and analysis of a Markov chain's state graph. 
This includes the calculation of the associated total mixing time or
its lower and upper bounds.

* Enumeration, uniform sampling, and counting of
  * binary matrices with prescribed row and column sums,
  * perfect matchings in bipartite graphs.


## Installation:

This library is developed and tested at Linux systems (primarily Ubuntu). 
However, it should be manageable to migrate the library to other operating systems.

The marathon software consists of a central library (header files) and a couple of  example applications that demonstrate how to use the library.
Some parts of marathon depend on various third party libraries and can therefore only be built when all dependencies are fulfilled.
 * [Boost](www.boost.org) 
 * [OpenBLAS](http://www.openblas.net/) (or another BLAS implemtation)
 * [Armadillo](arma.sourceforge.net/)
 * [Eigen](http://eigen.tuxfamily.org/)
 
The CMake installation script will automatically identify and build the components that can be build on your system. 

This instruction shows how to build the library and run the [transitionMatrix](./examples/transitionMatrix/) example at a fresh Ubuntu 16.04 system.

1. Install package requirements.

        sudo apt-get install git cmake g++ libboost-all-dev libblas-dev libarmadillo-dev libeigen3-dev

2. Build examples.

        git clone https://github.com/srechner/marathon.git
        cd marathon
        mkdir build
        cd build
        cmake ..
        make
        
4. Run an example. 

        ./examples/transitionMatrix/transtitionMatrix classical-switch "4,4,2,1;3,3,3,1,1" 
        
    The output should look like:
        
          19/20  1/60  1/60  1/60  0  0
          1/60  19/20  0  0  1/60  1/60
          1/60  0  19/20  1/60  1/60  0
          1/60  0  1/60  19/20  0  1/60
          0  1/60  1/60  0  19/20  1/60
          0  1/60  0  1/60  1/60  19/20
        
## Usage 

There are several examples that demonstrate how our library can be used.

 * [sample](./examples/sample/): Construct uniformly distributed binary matrices with bounded row and column sums.
 * [count](./examples/count/): Calculate the number of objects.
 * [enumerate](./examples/enumerate/): Enumerate the set of binary matrices with bounded row and column sums.
 * [transitionMatrix](./examples/transitionMatrix/): Construct the transition matrix of a Markov chain. 
 * [totalMixingTime](./examples/totalMixingTime/): Calculate the total mixing time of a Markov chain.
 * [spectralBound](./examples/spectralBound/): Calculate the lower and upper spectral bound on the total mixing time of a Markov chain.
