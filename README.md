## marathon 0.4

This C++ library is designed to support the analysis of Markov chain based sampling methods. It provides functions for the analysis of so-called state graphs. For an introduction into its functionality, see its introductional [article](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0147935).
If you use this library for your work, please cite

Steffen Rechner and Annabell Berger:
marathon: An open source software library for the analysis of Markov-Chain Monte Carlo algorithms   , PLoS ONE 2016


An API description can be found [here](./doc/marathon.pdf) (It's a start and will be gradually improved.)

### Current ToDo's:
* Improving Documentation / API Description
* Integrating of network analysis
* Integrating additional chains
* Setting up new experiments

### Requirements:

This library is developed and tested at Linux systems (primarily Ubuntu). In principle it should possible to migrate the library to other operating systems.

Building the marathon library requires the following libraries:
 * `Boost` 
 * a BLAS implementation (e.g. `OpenBLAS`)
 * `Arpack++`
 * `SuperLU`
 * `CUDA` (optional)

### Installation:

This instruction shows how to build the library and run the [MixingBounds](./applications/MixingBounds/) example at a fresh Ubuntu 16.04 system.

1. Install packages.

        sudo apt-get install git cmake g++ libboost-all-dev libblas-dev libarpack++2-dev libsuperlu-dev

2. Build the library.

        git clone https://github.com/srechner/marathon.git
        cd marathon
        mkdir build
        cd build
        cmake..
        
3. Install the library to the system path.

        make
        sudo make install

4. Run an application.

        MixingBounds swapBip "2,1,1;1,2,1" 1e-3

   The output should look like:

        number of states:          5
        number of transition arcs: 21
        total mixing time:         72
        lower spectral bound:      34.1803
        upper spectral bound:      102.206
        upper congestion bound:    183.971

