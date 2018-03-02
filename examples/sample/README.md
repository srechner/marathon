# sample

This example demonstrates how to use `RandomGenerator` classes to produce uniformly distributed binary matrices.

## Usage

    sample METHOD NUMBER [OPTION] INSTANCE
    
The parameter `METHOD` must be substituted by one of the following:

* `exact`: Use an exact sampling algorithm (not based on MCMC).

* `simple`: Use the first Markov chain defined by 
  > Steffen Rechner, Linda Strowick, and Matthias Müller-Hannemann. Uniform sampling of bipartite graphs with degrees in prescribed intervals. Journal of Complex Networks (2017). doi:10.1093/comnet/cnx059.
                
* `informed`: Use the second Markov chain defined by 
  > Steffen Rechner, Linda Strowick, and Matthias Müller-Hannemann. Uniform sampling of bipartite graphs with degrees in prescribed intervals. Journal of Complex Networks (2017). doi:10.1093/comnet/cnx059.

If `METHOD` is `simple` or `informed`, the optional parameter

    -s STEPS

must be used to specify a random walk of length `STEPS`. (Larger is better.)

The number of samples to produce is specified by the parameter `NUMBER`.

The parameter `INSTANCE` is a string of the form `l-u(,l-u)*;l-u(,l-u)*`. 
The semicolon separates the row sums from the column sums. The i-th pair on the left
side of the semicolon prescribes a lower and upper bound on the sum of
row i. In contrast, the j-th pair on the right side of the semicolon
defines the sum of column j. For convenience, a pair `l-u` can be
replaced by an integer `l` if `l=u`.
  
For example, the string `"1-2,2,2-3;0-2,0-1,1-1,1-3"` corresponds to

        lower bounds on row sums:    (1,2,2)
        upper bounds on row sums:    (2,2,3)
        lower bounds on column sums: (0,0,1,1)
        upper bounds on column sums: (2,1,1,3)
