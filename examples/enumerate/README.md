# enumerate

This example demonstrates how to use  an `Enumerator` class to list the set of binary matrices whose row and column sums lie in prescribed intervals.
## Usage

    Usage: enumerate INSTANCE

The parameter `INSTANCE` is a string of the form `l-u(,l-u)*;l-u(,l-u)*`.

The semicolon separates the row sums from the column sums. The i-th pair on the
left side of the semicolon defines lower and upper bounds on the sum of row i.
In contrast, the j-th pair on the right side of the semicolon defines the sum
of column j. For convenience, a pair `l-u` can be replaced by an integer `l` if `l` equals `u`.

For example, the string `1-2,2,2-3;0-2,0-1,1-1,1-3` corresponds to

        lower bounds on row sums:    (1,2,2)
        upper bounds on row sums:    (2,2,3)
        lower bounds on column sums: (0,0,1,1)
        upper bounds on column sums: (2,1,1,3)
