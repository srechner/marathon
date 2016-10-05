## marathon 0.4

This C++ library is designed to support the analysis of Markov chain Monte Carlo methods. It provides functions for the analysis of so-called state graphs. For an introduction into its functionality, see its introductional [article](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0147935).
If you want to use this library for your work, please feel free to contact me and cite

Steffen Rechner and Annabell Berger:
marathon: An open source software library for the analysis of Markov-Chain Monte Carlo algorithms   , PLoS ONE 2016.

An API description can be found [here](./doc/marathon.pdf) (It's a start and will be gradually improved.)

### Installation:

This library is developed and tested at Linux systems (primarily Ubuntu). However, it should be manageable to migrate the library to other operating systems.

The marathon software consists of a central library (header files) and a couple of applications that demonstrate how to use the library.
Some parts of marathon depend on various third party libraries and can therefore only be built when all dependencies are fulfilled.
Building the full library requires the following libraries:
 * `Boost` 
 * a BLAS implementation (e.g. `OpenBLAS`)
 * `Arpack++`
 * `SuperLU`
 * `Eigen`
 * `CUDA`
 * `Eigen`
 
However, the CMAKE installation script will automatically identify and build the components that can be build on your system. 

This instruction shows how to build the library and run the [MixingBounds](./src/apps/MixingBounds/) example at a fresh Ubuntu 16.04 system.

1. Install packages.

        sudo apt-get install git cmake g++ libboost-all-dev libblas-dev libarpack++2-dev libsuperlu-dev

2. Build the software.

        git clone https://github.com/srechner/marathon.git
        cd marathon
        mkdir build
        cd build
        cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME ..
        make
        
3. Install (optional).

        make install

4. Run an application.

        MixingBounds swapBip "2,1,1;1,2,1" 0.001

   The output should look like:

        number of states:          5
        number of transition arcs: 21
        total mixing time:         72
        lower spectral bound:      34.1803
        upper spectral bound:      102.206
        upper congestion bound:    183.971
