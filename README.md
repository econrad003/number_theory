# number_theory
relatively small programs for studying elementary number theory

## Description

This is a place for me to save programs that I write as experiments in elementary number theory.  (Don't expect anything particularly deep, but if you happen to find something useful, then great!)  Most of what goes here will probably be written in Python 3.

## License

The programs are licensed under an open MIT-style license.  The full license may be found in the *LICENSE* file in this repository and in program and module documentation.  If you do find something to be useful or fun or surprising, you are not required to share your efforts and modifications, but I encourage you to do so.

## Python programs and modules

* *primality.py* - a module for studying relatively small prime numbers (say four or fewer decimal digits, depending on time and memory).  To test the module, run it as the main program.  The test code can be used as guide to importing this as a module.  The module maintains a list of consecutive "small" primes in the class Primes -- the list is extended (by sieving) as needed.  In addition to the class, the module defines a decorator called *multiplicative* which turns a function f(p,e), where p is a positive prime and e is a positive integer representing an exponent into a multiplicative number-theoretic (or "arithmetic") function defined on the integers.  A few standard arithmetic functions are defined in the module.

* *aliquot.py* - a module for studying aliquot sequences - *coming soon!*
