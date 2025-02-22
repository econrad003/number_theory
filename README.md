# number_theory - elementary number theory

This repository is a collection of relatively small programs for studying elementary number theory.

## Description

This is a place for me to save programs that I write as experiments in elementary number theory.  (Don't expect anything particularly deep, but if you happen to find something useful, then great!)  Most of what goes here will probably be written in Python 3.

## Recent Updates

2 Feb 2025 (one module added) - Module *harmonic.py* approximates positive real numbers by series based on greedily chosen subsequences of the harmonic sequence.  For numbers between 0 and 1, the result is as for Sylvester's Egyptian fraction algorithm.  For large positive numbers, the sequences become unmanageable.  This is apt to be of more theoretical than practical significance.  As with *cfrac.py* and *sylvester.py* this contains a *self\_test()* method and a few useful support methods.

31 Jan 2025 (two modules added) - Module *cfrac.py* finds regular continued fraction representations for real numbers.  Module *sylvester.py* finds Egyptian fraction representations for real numbers.  See the program and function documentation for details.  Both modules contain a method *self\_test()* which is executed as a main routine -- the *self\_test()* method provides some usage examples.

27 May 2024 - The Division Algorithm and the Euclidean Algorithm for Gaussian integers have been incorporated in *gaussian_int.py*.  A demo which lists Gaussian primes (*list_gaussian_primes.py*) and some sample output have also been added to this repository.

## License

The programs are licensed under an open MIT-style license.  The full license may be found in the *LICENSE* file in this repository and in program and module documentation.  If you do find something to be useful or fun or surprising, you are not required to share your efforts and modifications, but I encourage you to do so.

## Python programs and modules

* *primality.py* - a module for studying relatively small prime numbers (say four or fewer decimal digits, depending on time and memory).  To test the module, run it as the main program.  The test code can be used as guide to importing this as a module.  The module maintains a list of consecutive "small" primes in the class Primes -- the list is extended (by sieving) as needed.  In addition to the class, the module defines a decorator called *multiplicative* which turns a function f(p,e), where p is a positive prime and e is a positive integer representing an exponent into a multiplicative number-theoretic (or "arithmetic") function defined on the integers.  A few standard arithmetic functions are defined in the module.

* *aliquot.py* - a program for studying aliquot sequences - usage information is in the program docstring and can be accessed using in help using the command "`python3 aliquot.py -h`".  For a typical example, try something like `python3 aliquot.py 360` for a plot of the aliquot sequence for 360.  There are some checks analogous to recursion limits to help avoid accidentally using all your computer memory, as well as ways to soften or circumvent these checks.

* *gaussian_int.py* - a module which implements Gaussian integers (class GaussianInt -- complex numbers whose real and imaginary parts are both integers) and Gaussian rational numbers (class GaussianFrac -- quotients of Gaussian integers).  Basic operations (addition, subtraction, multiplication, division, and exponentiation) are implemented.  In addition, the GaussianInt class incorporates rudimentary primality testing using some of the features implemented in module *primality.py*.  Running the module as a main program executes a number of program checks.  Since this module uses *math.isqrt()", it requires Python 3.8 or later.  The Gaussian integer division algorithm (division with remainder) and the Euclidead GCD algorithm have been programmed and seem to be working.

* *list_gaussian_primes.py* - a demo program to list Gaussian primes.  Run it with the *-h* option for more details.  Sample output (console output for Gaussian primes less than 100 in complex norm and CSV output for Gaussian primes less than 500 in complex norm, both sorted in increasing norm order) is in this repository.  (Note: I count zero as a prime, so it appears first in both lists.  The Gaussian units \[1, -1, i, -i\] are not counted as prime.)

```
        python3 list_gaussian_pries.py -h
```

## Examples

* *Gaussian_primes_to_100.txt* - Gaussian primes (including 0, but not units) with complex norm less than 100 in display format, with headings.

* *Gaussian_primes_to_500.csv* - Gaussian primes (including 0, but not units) with complex norm less than 500 in comma-separated variable (*i.e.* simple spreadsheet) format, with headings.

* *plot_aliquot(1264460)-sociable-period4.png* - a plot of the aliquot sequence for 1264460, a sociable number with a period of 4.  To reproduce this example, you need to soften one of the checks:

```
        python3 aliquot.py 1264460 --largest_value 2000000
```

>  The two million argument allows terms in the sequence to be as large as two million.  By default, the program will terminate with an exception if the largest value or the largest prime in the sieve exceeds 100,000.  The --largest_value and --largest_prime optional arguments allow these limits to be increased.  The following display shows the console output:

```
        ~/Projects/Number Theory> python aliquot.py 1264460 --largest-value 2000000
        WARNING (aliquot): sieved up to 3947 from 19
        WARNING (aliquot): sieved up to 43783 from 3947
        [1547860, 1727636, 1305184, 1264460, 1547860]
        The sequence is periodic.
        Period 4: repeats with sociable number 1547860 
```
 
* *plot_aliquot(276)-unknown.png* - a plot for the aliquot sequence for 276.  It is not known whether this sequence is bounded. According to ProofWiki, the 469th term is 45 digits long.  Our plot only found the first 26 digits.  This is the smallest unsolved example.  Here is the console output:

```
        ~/Projects/Number Theory> python aliquot.py 276 -l 25 --log -L 1000000
        WARNING (aliquot): sieved up to 31 from 19
        WARNING (aliquot): sieved up to 577 from 31
        WARNING (aliquot): sieved up to 929 from 577
        WARNING (aliquot): sieved up to 2441 from 929
        WARNING (aliquot): sieved up to 16703 from 2441
        WARNING (aliquot): Maximum length 25 exceeded.
        [396, 696, 1104, 1872, 3770, 3790, 3050, 2716, 2772, 5964, 10164, 19628, 19684, 22876, 26404, 30044, 33796, 38780, 54628, 54684, 111300, 263676, 465668, 465724, 465780]
        The sequence does not repeat or terminate after 25 terms.
```

>  The -l parameter told the program to stop after 25 terms.  (It finds 26.)  The -L parameter says that terms should not exceed ome million.  The --log flag says that the plot should be logarithmic. The console output tells us that the last three evaluated terms were between 400,000 and 500,000.

* *plot_aliquot(360)-terminating.png* and *plot_aliquot(360)-terminating-log.png* - a plot and a logarithmic plot of the aliquot sequence for 360, and abundant number as the sum s(360)=396 of its proper divisors exceeds 360.  Its aliquot sequence terminates after 29 terms.  It reaches a peak of 6490 (s(6464)=6490) before collapsing.  The console output tells the story:

```
        ~/Projects/Number Theory$ python aliquot.py 360
        WARNING (aliquot): sieved up to 211 from 19
        WARNING (aliquot): sieved up to 1499 from 211
        WARNING (aliquot): sieved up to 2441 from 1499
        [810, 1368, 2532, 3404, 2980, 3320, 4240, 5804, 4360, 5540, 6136, 6464, 6490, 6470, 5194, 4040, 5140, 5696, 5734, 3194, 1600, 2337, 1023, 513, 287, 49, 8, 7, 1, 0]
        The sequence is terminating.
```
