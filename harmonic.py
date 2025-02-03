"""
harmonic.py - harmonic decomposition
Eric Conrad

BACKGROUND

    The harmonic series (1+1/2+1/3+1/4+1/5+...) diverges to infinity.  Many
    calculus textbooks today cite a beautiful proof of this fact that goes
    back to the French monk Nicole Oresme (1325-1382).  The terms are monotone
    decreasing and the partial sums are monotone increasing.

    The terms of the series are positive and decrease strictly to zero.  Since
    the series diverges, we can approximate any positive real number with some
    subsequence.  One way of doing this is using a greedy algorithm -- namely
    a modification of Sylvester's Egyptian fraction algorithm.  That is the
    algorithm that we implement here.

DESCRIPTION

    Given a nonnegative number x and a positive number ε we want to find
    a decreasing sequence of unit fractions whose sum is in the interval
        (x-ε, x].

    Let s = [], the empty sequence, and let S=0.  (S is the sum of the sequence
    s.)

    m = 0
    While S ≤ x-ε:
        let y = x - S;
        let n = max(m+1, ceil(1/y))
        s.append(1/n)
        S += 1/n
        m = n + 1

    Return s.

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

def harmonic(x:Real, epsilon:Real=0, debug=False) -> list:
    """harmonic representation of a non-negative real number

    Returns a an increasing list [a, b, ..., n] of integers.  The sum of
    reciprocals:
        1/a + 1/b + ... + 1/n
    is greater than x-epsilon and less than or equal to x.

	For values between 0 and 1, the result is the same as that given by
	Sylvester's algorithm.  For values larger than or equal to 2, the first
	few terms will agree with the harmonic sequence.  Convergence will be
	very slow for values much larger than 2.

    DEBUG

        Setting debug to False admits a zero error.  This should only be
        used if x is an instance of an exact rational type.  (Don't try
        this with type float!)

    BUGS
        The result is subject to the typical disclaimers involving floating
        point rounding errors.

    SEE ALSO:
        cfrac.py - for continued fraction estimates
        sylvester.py - for Sylvester's algorithm estimates
    """
    if x < 0:
        raise ValueError("The number must be non-negative")
    if debug:
        assert epsilon > 0, "This value should be positive for floating point"
    if epsilon < 0:
        raise ValueError("The error must be non-negative")
    result = list()
    m = 0                           # the last denominator

    while x >= epsilon and x > 0:
        # print(x)
        n = max(m+1, ceil(1/x))
        result.append(n)
        x -= Fraction(1, n)
        m = n

    return result

def value_of(result:list) -> Fraction:
    """return the value of the decomposition"""
    x = 0
    for y in result:
        x += Fraction(1, y)
    return x

def convergents_of(result:list) -> list:
    """return the convergents of the decomposition"""
    convergents = list()
    x = 0
    for y in result:
        x += Fraction(1, y)
        convergents.append(x)
    return convergents

def texform(result:list) -> str:
    """return the result as a TeXForm expression"""
    if len(result) == 0:
        return "0"
    s = str()
    s += "\\frac{1}{" + str(result[0]) + "}"
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
        print("N.B.:")
        print("  Note that 3, π, and 4 use respectively the first 10, 12, and")
        print("  30 terms of the harmonic series before homing in with small")
        print("  terms.  With larger numbers, it only gets worse!")

    print("Test cases for", basename(__file__))
        # 5/3 = 1 + 2/3
    x = Fraction(5, 3)
    result = harmonic(x, debug=False)
    display(x, 0, result, "1+2/3")
        # 355/113 (a very nice rational approximation of π)
    x = Fraction(355, 113)
    result = harmonic(x, debug=False)
    display(x, 0, result, "355/113")
        # 3
    x = Fraction(3, 1)
    result = harmonic(x, debug=False)
    display(x, 0, result, "3 N.B.")
        # 3
    x = Fraction(4, 1)
    result = harmonic(x, debug=False)
    display(x, 0, result, "4 N.B.")
        # √2
    x = sqrt(2)
    err = 1e-4
    result = harmonic(x, err)
    display(x, err, result, "≈ √2")
        # (1+√5)/2
    x = (1+sqrt(5))/2
    result = harmonic(x, err)
    display(x, err, result, "≈ (1+√5)/2")
        # e
    x = e
    result = harmonic(x, err)
    display(x, err, result, "≈ e")
        # π
    x = pi
    result = harmonic(x, err)
    display(x, err, result, "≈ π")
    nota_bene()

if __name__ == "__main__":
    self_test()
