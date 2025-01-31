"""
sylvester.py - Sylvester's greedy Egyptian fraction algorithm
Eric Conrad

BACKGROUND

    A unit fraction is simply a positive fraction whose numerator is one.  In
    other words, x is a unit fraction if and only if 1/x is a positive integer.
    In ancient Egypt, a standard way of writing a number between zero and one
    was to write it as a sum consisting of distinct unit fractions, the
    fraction 2/3 and sometimes also the fraction 3/4.  In fact, every positive
    real number can be expressed as a convergent sum of distinct unit
    fractions.  (This representation is not unique.)

    For example:

        1 = 1/2 + 1/3 + 1/6
            Proof: 1/2 + 1/3 + 1/6 = (3+2+1)/6 = 6/6 = 1
        2/3 = 1/3 + 1/4 + 1/12
            Proof: 1/3 + 1/4 + 1/12 = (4+3+1)/12 = 8/12 = 2/3

    In number theory, an Egyptian fraction is defined as a sum of distinct
    unit fractions.  The terms are usually written as a descending sequence.

SYLVESTER'S ALGORITHM

    J J Sylvester described a greedy algorithm for finding an Egyptian fraction
    representation of any real number between zero and one.  The algorithm
    seems to have been known to some ancient Egyptian scribes and (much later)
    to Fibonacci, but Sylvester is almost certainly the first to write down a
    description of the algorithm.  It is a modification of the algorithm used
    to convert a real number between 0 and 1 into a regular continued fraction.
    We can describe the algorithm in terms of a state machine which contains
    a numberin the interval [0,1].  We press a button, and the machine returns
    the next term.

    INITIALIZATION:
        let x be a number in the interval [0,1]

    GET_NEXT_TERM:
        if x is 0:
            return 0;
        let y = ceil(1/x);
        x = x - 1/y
        return 1/y

    If the initial value of x is a rational number, then the algorithm will
    terminate, eventually returning zeros.  If the initial value is irrational,
    the sequence does not terminate.

    With rational input, the trailing terms of the fraction often have huge
    denominators, but the number of terms (excluding the integer part) never
    exceeds the numerator of the proper part of the fraction.

    The result is a best rational approximation of the input.

IMPLEMENTATION

    In our implementation, we provide a function
        sylvester(x, err)
    where x is a non-negative real number and err is a positive real upper
    bound for the error in the representation.  It returns a finite sequence
    consisting of an integer and a series of unit fractions whose sum will be
    less than or equal to the input x. Specifically
        |x - sum| <= err

LICENSE
    Copyright 2024 by Eric Conrad.

    Permission is hereby granted, free of charge, to any person obtaining
    a copy of this software and associated documentation files (the
    “Software”), to deal in the Software without restriction, including
    without limitation the rights to use, copy, modify, merge, publish,
    distribute, sublicense, and/or sell copies of the Software, and to
    permit persons to whom the Software is furnished to do so, subject
    to the following conditions:

    The above copyright notice and this permission notice shall be included
    in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
    OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
    IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
    CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
    TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
from math import floor, ceil
from numbers import Real
from fractions import Fraction

def sylvester(x:Real, epsilon:Real=0, debug=True) -> list:
    """Sylvester's representation

    Returns a list [u, a, b, ..., n] of integers.  The first term is the
    greatest integer in (or floor of) x.  The remaining terms are the
    denominators of the unit fractions in the representation.

    The result satisfies the following inequality:

        |      1     1           1     |
        | u + --- + --- + ... + --- - x| <= epsilon
        |      a     b           n     |

    BUGS
        The result is subject to the typical disclaimers involving floating
        point rounding errors.

    SEE ALSO:
        cfrac.py - for continued fraction estimates
    """
    if debug:
        assert epsilon > 0, "This value should be positive for floating point"
    if epsilon < 0:
        raise ValueError("The error must be non-negative")
    result = list()
    u = floor(x)
    result.append(u)
    x -= u
    while x > epsilon:
        n = ceil(1/x)
        result.append(n)
        x -= Fraction(1, n)
    return result

def value_of(result:list) -> Fraction:
    """return the value of the decomposition"""
    x = result[0]
    for y in result[1:]:
        x += Fraction(1, y)
    return x

def convergents_of(result:list) -> list:
    """return the convergents of the decomposition"""
    convergents = list()
    x = result[0]
    convergents.append(x)
    for y in result[1:]:
        x += Fraction(1, y)
        result.append(x)
    return result

def texform(result:list) -> str:
    """return the result as a TeXForm expression"""
    s = str()
    s += str(result[0])
    for y in result[1:]:
        s += "+\\frac{1}{" + str(y) + "}"
    return s

def self_test():
    """a simple self test"""
    from math import sqrt, pi, e
    from os.path import basename

    def display(x, err, result, note):
        """display a result"""
        print("-"*72)
        xstar = value_of(result)
        if err > 0:
            print(x, xstar, result, x-xstar, note)
        else:
            print(x, xstar, result, note)
        convergents = convergents_of(result)
        print("  convergents:", *convergents)
        print("      TeXForm:", texform(result))
        assert abs(x-xstar) <= err

    def nota_bene():
        """display a little observation"""
        print("-"*72)
        print("NB: In the course of these investigations, we unconvered a few")
        print("    estimates of π that deserve some mention.")
        print("         22/7 ≈", 22/7)
        print("      355/113 ≈", 355/113, "(my favorite estimate)")
        print("         25/8 =", 25/8, "(from Sylvester's algorithm)")
        print("     1533/488 ≈", 1533/488, "(from Sylvester's algorithm)")
        print("            π ≈", pi, "(Pythonic estimate)")
        print("-"*72)

    print("Test cases for", basename(__file__))
        # 5/3 = 1 + 2/3
    x = Fraction(5, 3)
    result = sylvester(x, debug=False)
    display(x, 0, result, "1+2/3")
        # 355/113 (a very nice rational approximation of π)
    x = Fraction(355, 113)
    result = sylvester(x, debug=False)
    display(x, 0, result, "≈ π")
    print("\tπ - 355/113 ≈", pi - 355/113)
        # √2
    x = sqrt(2)
    err = 1e-4
    result = sylvester(x, err)
    display(x, err, result, "≈ √2")
        # (1+√5)/2
    x = (1+sqrt(5))/2
    result = sylvester(x, err)
    display(x, err, result, "≈ (1+√5)/2")
        # e
    x = e
    result = sylvester(x, err)
    display(x, err, result, "≈ e")
        # π
    x = pi
    result = sylvester(x, err)
    display(x, err, result, "≈ π")
    nota_bene()

if __name__ == "__main__":
    self_test()
