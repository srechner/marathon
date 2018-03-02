# transitionMatrix

This example demonstrates how to construct the transition matrix of a specified Markov chain.

## Usage

    transitionMatrix CHAIN INSTANCE
    
The parameter `CHAIN` must be substituted by one of the following strings: 

    'classical-switch':
       Markov chain defined by 'Kannan et al. Simple Markov-chain  algorithms
       for generating bipartite graphs and tournaments. Random Structures and
       Algorithms 14 (1997), 293–308'.
       
    'edge-switch':
       Variant of the Markov-chain suggested by 'Kannan et al.' based on an
       informed edge selection at the cost of a larger memory consumption.
              
    'curveball':
       Markov chain defined by 'Strona et al. A fast and unbiased procedure to
       randomize ecological binary matrices with fixed row and column totals.
       Nature communications 5 (2014).'

An `INSTANCE` has the form `r,r,...,r;c,c,...,c`, where `r` and `c` are positive integers. 
While the i-th occurrence of `r` defines the sum of row i, the j-th
occurrence of `c` is the sum of column j. 

For example, the instance `2,2,2;1,2,1,2` corresponds to the

            row sums:    (2,2,2)
            column sums: (1,2,1,2)
