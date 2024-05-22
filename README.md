# number_theory
relatively small programs for studying elementary number theory

## Description

This is a place for me to save programs that I write as experiments in elementary number theory.  (Don't expect anything particularly deep, but if you happen to find something useful, then great!)  Most of what goes here will probably be written in Python 3.

## License

The programs are licensed under an open MIT-style license.  The full license may be found in the *LICENSE* file in this repository and in program and module documentation.  If you do find something to be useful or fun or surprising, you are not required to share your efforts and modifications, but I encourage you to do so.

## Python programs and modules

* *primality.py* - a module for studying relatively small prime numbers (say four or fewer decimal digits, depending on time and memory).  To test the module, run it as the main program.  The test code can be used as guide to importing this as a module.  The module maintains a list of consecutive "small" primes in the class Primes -- the list is extended (by sieving) as needed.  In addition to the class, the module defines a decorator called *multiplicative* which turns a function f(p,e), where p is a positive prime and e is a positive integer representing an exponent into a multiplicative number-theoretic (or "arithmetic") function defined on the integers.  A few standard arithmetic functions are defined in the module.

* *aliquot.py* - a program for studying aliquot sequences - usage information is in the program docstring and can be accessed using in help using the command "`python3 aliquot.py -h`".  For a typical example, try something like `python3 aliquot.py 360` for a plot of the aliquot sequence for 360.  There are some checks analogous to recursion limits to help avoid accidentally using all your computer memory, as well as ways to soften or circumvent these checks.

## Examples

* *plot_aliquot(1264460)-sociable-period4.png* - a plot of the aliquot sequence for 1264460, a sociable number with a period of 4.  To reproduce this example, you need to soften one of the checks:
  ```
  python3 aliquot.py 1264460 --largest_value 2000000
  ```
  The two million argument allows terms in the sequence to be as large as two million.  By default, the program will terminate with an exception if the largest value or the largest prime in the sieve exceeds 100,000.  The --largest_value and --largest_prime optional arguments allow these limits to be increased.  The following display shows the console output:
  ```
  ~/Projects/Number Theory> python aliquot.py 1264460 --largest-value 2000000
  WARNING (aliquot): sieved up to 3947 from 19
  WARNING (aliquot): sieved up to 43783 from 3947
  [1547860, 1727636, 1305184, 1264460, 1547860]
  The sequence is periodic.
  Period 4: repeats with sociable number 1547860
  ```
