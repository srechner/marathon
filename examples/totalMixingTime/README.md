# totalMixingTime

This example demonstrates how to calculate the total mixing time from a state graph.

## Usage

    totalMixingTime CHAIN EPSILON INSTANCE
    
The parameter `CHAIN` must be substituted by one of the following strings: 

* `classical-switch`: Use the Markov chain defined in
  > Kannan et al. Simple Markov-chain  algorithms for generating bipartite graphs and tournaments. Random Structures and Algorithms 14 (1997), 293–308.
       
* `edge-switch`: Use a variant of the Markov chain defined in 
  > Kannan et al. Simple Markov-chain  algorithms for generating bipartite graphs and tournaments. Random Structures and Algorithms 14 (1997), 293–308.
  
  The transition rules of this Markov chain are based on a more informed edge selection at the cost of a large memory consumption.
  
* `curveball`: Use the Markov chain defined in
  > Strona et al. A fast and unbiased procedure to randomize ecological binary matrices with fixed row and column totals.        Nature communications 5 (2014).

The parameter `EPSILON` is a real number between zero and one. It specifies the target distance to the uniform distribution. 

An `INSTANCE` has the form `r(,r)*;c(,c)*`, where `r` and `c` are positive integers. 
While the i-th occurrence of `r` defines the sum of row i, the j-th
occurrence of `c` is the sum of column j. 

For example, the instance `2,2,2;1,2,1,2` corresponds to the

            row sums:    (2,2,2)
            column sums: (1,2,1,2)
