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
of a previous run with N _prev_>=N-1, in order to generate a good initial guess, unless N=1, which has a builtin initial guess.
Two different approaches are tried, in one the set of parameters C is optimized with an iterative method and in the
other they are directly calculated using the recent set of betas.

The input is given in a JSON document as follows:
```
{
  "number_of_terms" : The integer N above
  "max_iterations" : An integer limiting the maximum number of iterations per test value of N
  "guess": A URL to the output of s2g used with N>=number_of_terms-1
}
```

Most fields have default values:

```
{
  "max_iterations" : 1024,
}
```


In addition, the JSON fields can be specified on the command line.

#### s2g output

The output is printed to standard output as JSON with three sections:
```
{
    "input": {},
    "program_info": {},
       "C_conjugate_grandient" :{ 
            error: number, 
            "terms": [
               { C : numer, beta : number },
               ...
               { C : number, beta: Number },
            ]
       },
       "C_implied": { 
        error: number, 
        "terms": [
           { C : numer, beta: number },
           ...
           { C : number, beta: Number },
        ]
       }
    }
    "best_method": name of the best method
}    
```

The input section is just the input variables. The program_info contains
various information about the program - which compiler was used, which CPU
was used, versions of libraries used, running time, etc. 
The program uses two methods to calculate C, output of both is 
given with the final error function result, plus a determination which one was better
based on that.




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
   "test_points" : Integer number of test points
  "max_test_error" : maximum error allowed at any test point

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

The test_points array contain
the tests done - the x value, the sto value at that x, and the estimated 
sto value. The difference is in the error field.
###  integrate - Calculate 3- and 4- center STO integrals