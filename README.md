# STO-values
Calculate highly accurate sample values of 3- and 4- center electronic integrals using STO functions.

This package provides three programs. 

1. _s2g_: create an accurate estimate of an STO as a sum of GTOs
1. _investigate_ : Calculate the accuracy of the result of (1)
1. _integrate_ : Calculate 3- and 4- center STO integrals using the output of (1)

Program input and output is in JSON format.

The programs can be compiled and run or accessed as a cloud service at <> with the same input and output expected.

###  s2g - accurate estimate of an STO as a sum of GTOs

The goal is to find a solution to 

```exp(-x) = sum ( i=1...N )  { Ci * exp(-Bi*x^2) }```

That is, finding for a given N the following:

1. coefficients list Ci , i=1..N
1. exponent list Bi, i=1..N

As an exact solution does not exist for a finite N, the program tries to get an approximation that is good
up to a given criteria. There is no guarantee the program will succeed in finding a solution.

The accuracy of the approximation is given in root sum of squares of errors. 

In addition, you can specify the maximum error allowed in any test point. If this is not achieved, the proposed estimate
in considered inaccurate and the program will continue to iterate.

#### s2g input

The program works by using an iterative method to find Ci and Bi until
they do not improve any more. When solving for a given N, the program needs the output
of a previous run with N _prev_=N-1, in order to generate a good initial guess, unless N=1, which has a builtin initial guess.

The input is given in a JSON document as follows:
```
{
  "number_of_terms" : The integer N above
  "max_iterations" : An integer limiting the maximum number of iterations per test value of N
  "accuracy" : average accuracy desired, in root integral of error squared.
  "test_points" : Integer number of test points
  "max_test_error" : maximum error allowed at any test point
  "guess": A URL to the output of s2g used with N=number_of_terms-1
}
```

Most fields have default values:

```
{
  "max_iterations" : 1024,
  "accuracy" : 3.2e-8,
  "test_points" : 1024,
  "max_test_error" : 1e-10
}
```


In addition, the JSON fields can be specified on the command line.

#### s2g output

The output is printed to standard output as JSON with three sections:
```
{
    "input": {},
    "program_info": {},
    "terms": [
       { C : numer, b : number },
       ...
       { C : number, b: Number },
    ],
    "accuracy" : number
    "test_points": [
       {
        x: number,
        sto: number,
        estimate: number,
        error: number
        },
        ...
    ] 
}    
```

The input section is just the input variables. The program_info contains
various information about the program - which compiler was used, which CPU
was used, versions of libraries used, running time, etc. The accuracy field 
contains the root-squared-errors of all test point. The test_points contain
the tests done - the x value, the sto value at that x, and the estimated 
sto value. The difference is in the error field.





###  investigate - explore the accuracy of the s2g output

#### input

Using the output of s2g, this program produces values of 

```exp(-alpha*x)```

and

```sum ( i=1...N )  { Ci * exp(-Bi*x^2) }```

And the error between them. You can ask for multiple points
of comparison in various distributions. 


The input is given in a JSON document as follows:
```
{
  "sto_exponent": A number representing 'alpha' above,
  "estimator": A URL to the output of s2g ,
}
```
#### output

The output is printed to standard output as JSON with three sections:
```
{
    "input": {},
    "program_info": {},
    "N" : integer
    "terms": [
       { C : numer, b : number },
       ...
       { C : number, b: Number },
    ],
    "test_points": [
       {
        x: number,
        sto: number,
        estimate: number,
        error: number
        },
        ...
    ] 
}    


###  integrate - Calculate 3- and 4- center STO integrals