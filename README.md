# bipartite-partitions
Code to generate integer partitions of bipartite numbers

J.T. Kindt, Emory University, based on an algorithm developed in collaboration with C. Weeden

Generates integer partitions of bipartite numbers (NA,NB), motivated by statistics of
simulations of molecular clusters with two components.  

These are represented as matrices with row and column indices ranging from 0
to NA and 0 to NB; each matrix element represents the number of clusters whose
composition is given by its row and column.  The sum of the product of all
matrix elements with their row index will always equal NA, and likewise for
the column index and NB.  

Algorithm works by starting with the integer
partitions of Nmax=NA+NB and working through all the ways of distributing A
molecules among the distribution of cluster sizes corresponding to that
integer partition.  

Input is Nmax, NA.

Output is not in a regular "alphabetical", is grouped by the "parent" integer
partitions of Nmax.
