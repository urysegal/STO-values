# STO-values
Calculate highly accurate sample values of 3- and 4- center electronic integrals using STO functions.

This package provides three programs. 

1. _s2g_: create an accurate estimate of an STO as a sum of GTOs
1. _investigate_ : Calculate the accuracy of the result of (1), with emphasis on cusp and tail
1. _integrate_ : Calculate 3- and 4- center STO integrals using the output of (1)

Program input and output is in JSON format.

The programs can be compiled and run or accessed as a cloud service at <> with the same input and output expected.

###  s2g - accurate estimate of an STO as a sum of GTOs

The goal is to find, given a positive number ```a```, a solution to 

```exp(-a*x) = sum ( i=1...N )  { Ci * exp(-Bi*x^2) }```

That is, finding  the following:

1. number of elements to sum N
1. coefficients list Ci , i=1..N
1. exponent list Bi, i=1..N

As an exact solution does not exist for finite N, the program tries to get an approximation that is good
up to a given criteria. There is no guarantee the program will succeed in finding a solution.

The accuracy of the approximation is given in root sum of squares of errors. The test 
points are chosen with emphasis on the cusp and tail of the STO.

In addition, you can specify the maximum error allowed in any test point. If this is not achieved, the proposed estimate
in considered inaccurate and the program will continue to iterate.

#### s2g input

The program works by picking a trail N, and then using an iterative method to find Ci and Bi until
they do not improve any more. If the desired accuracy has not reached, the program
increases N and tries to find a new set of Ci and Bi.

You can control the maximum N and also the maximum number of iteration
per a trial N. 

The input is given in a JSON document as follows:
```
{
  "sto_exponent": A number representing 'a' above
  "max_number_of_terms" : An integer limiting the upper value of N to try
  "max_iterations" : An integer limiting the maximum number of iterations per test value of N
  "accuray" : accuracy desired, in root sum of squares.
  "test_points" : Integer number of test points
  "max_test_error" : maximum error allowed in any test point
}
```

Sample and default values are:

In addition, the JSON fields can be specified on the command line.

###  investigate - explore the accuracy of the s2g output


###  integrate - Calculate 3- and 4- center STO integrals