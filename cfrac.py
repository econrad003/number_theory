"""
cfrac.py - greedy continued fraction algorithm
Eric Conrad

BACKGROUND

    A regular continued fraction is a number written in the following
    form:

        b0 +  1
             -------------------------
              b1 + 1
                  --------------------
                   b2 + 1
                       ---------------
                        b3 + 1
                            ----------
                             ...

    The number b0 is an integer and the sequence (b1, b2, b3, ...) is a
    sequence of positive integers.  The sequence may be terminating or
    infinite.  It can be shown that the sequence always converges and, for
    for each real number, the continued fraction is essentially unique.

    There is a technical condition: if the last term of the sequence is 1,
    we have a simpler representation, for example:
        1 + 1
           ---  = 2
            1

    The denominator sequence is intimately tied in with the Euclidean GCD
    algorithm.

CONVERGENCE

    Given a sequence like the sequence [1;2,2,2...] which converges to the
    square root of 2, the standard convergents are the finite continued
    fractions:
        1
        [1;2] = 1+1/2
        [1;2,2] = 1+1/(2+1/2)
        [1;2.2,2] = 1+1/(2+1/(2+1/2))
        etc.

    These can be shown to produce a sequence which alternates between low
    estimates and high estimates.  Using the square root of 2, we can simolify
    the example by adding 1 to extend the repeat to the integer term:
        1+√2 = [2;2,2,2,...]
    For this sequence, we have:
        x(n) = [2;2,2,2,...,2]                       (n terms)
        x(n+1) = [2;2,2,2,...,,2,2]                  (n+1 terms)
               = [2;1/x(n)] = 2+1/x(n)
    Looking at the first few terms:
        2 < 1+√2                                     2
        [2;2] = 2+1/2 = 5/2 > 1+√2                   2.5
        [2;2,2] = 2+2/5 = 12/5 < 1+√2                2.4
        [2;2,2,2] = 2+5/12 = 29/12 > 1+√2            2.41666...
        [2;2,2,2,2] = 2+12/29 = 70/29 < 1+√2         2.41379...
    It can be proven that this alternating pattern of low estimates and
    high estimates holds for all regular continued fractions.  In addition,
    it can be shown that, by taking enough terms, the difference between
    consecutive convergents can be made as close to zero as we wish.

    The convergents of a regular continued fraction are all rational numbers.
    So finite regular continued fractions represent (and are in fact equal to
    rational numbers.  Moreover, every rational number has a regular continued
    fraction representation.  (If the last term of the reprentation is 1, there
    is another representation.  If we don't admit finite regular continued
    fractions whose lst term is 1 -- except for the integer [1], then the
    representation is unique.)

    If the representation is repeating, then the continued fraction represents
    and converges to a quadratic surd (i.e. a root of a quadratic equation).
    Conversely, every quadratic surd has a unique repeating continued fraction
    representation.  In the examples in the self-test, we have:
        √2 = [1;2,2,2,...]                  square root of 2
        (1+√5)/2 = [1;1,1,1,...]            golden ratio

    If a number is not the root of a quadratic equation, for example, the
    transcendental numbers π and e, or the cube root of 3 (an algebraic
    number), then its regular continued fraction is non-repeating.  Again, the
    representation is unique.

RATE OF CONVERGENCE

    Another continued fraction of historical interest is Bombelli's continued
    fraction for the square root of thirteen.  Bombelli presented it in a
    different form in a treatise (Bombelli's Algebra) which was published
    posthumously in 1599.  As a regular continued fraction its representation
    is:
        √13 = [3;1,1,1,1,6,...]
    The digits [1,1,1,1,6] (four ones and a six) are repeating in the
    representation.  Python displays its estimate of √13 as 3.605551275463989.
    Let's look at some convergents:
        3/1 = 3                             LOW     [3]
        4/1 = 4                             HIGH    [3;1] = [4]
        7/2 = 3.5                           LOW     [3;1,1] = [3;2]
        11/3 = 3.66666...                   HIGH    [3:1,1,1] = [3;1,2]
        18/5 = 3.6                          LOW     [3;1,1,1,1]             (*)
        119/33 = 3.6060606...               HIGH    [3;1;1;1,1,6]
        137/38 = 3.605263157...             LOW     [3;1,1,1,1,6,1]
        256/71 = 3.605633802...             HIGH    [3;1,1,1,1,6,1,1]
        393/109 = 3.605504587...            LOW     [3;1,1,1,1,6,1,1,1]
        649/180 = 3.605555555...            HIGH    [3;1,1,1,1,6,1,1,1,1]   (*)
        4287/1189 = 3.605550883...          LOW     [3;1,1,1,1,6,1,1,1,1,6]
    Now look at the two starred lines.  They come immediately before the lines
    whose convergents end with a six.  Note that there is a significant
    improvement in convergence at these points.

    Now note the first few terms in the representation of π:
        π ≈ [3;7,15,1,292,...]          (non-repeating!)
    Except for the fourth term, the first five terms are larger than 1. The
    fifth term (292) is quite large, so four terms gives the memorable
    representation:
        π ≈ [3;7,15,1] = 355/113
    This representation gives the value of π with eight significant digits of
    precision using a ratio of two three-digit integers.  The next convergent
    is:
        π ≈ [3;7,15,1,292] = 103993/33102
    This gives us nine signficant digits, but the ratio is six digits over five
    digits.  355/113 requires a total of six digits, less than a local telephone
    number in the US.  103993/33102 has a total of 11 digits making it akin to
    an international telephone number in complexity.

    In general, large terms are omens of jumps in the rate of convergence.
    One jump occurs before the large term is added.  The second jump occurs
    with the addition of the term, but comes with an increase in the complexity
    of the rational approximation.

    The slowest convergence is exhibited by the golden ratio:
        (1+√5)/2 = [1;1,1,1,1,...]          (repeating ones)

CAVEAT EMPTOR

    Python typically implements the float type using "double precision floating
    point" which gives about 14 signficant digits of precision.  If floating
    point is used as input, the resulting continued fraction will show a small
    loss of precision which will increase with the number of terms in the
    estimate.  In any case, the number of significant digits in the last
    convergent will be at most 14.

    The moral: use exact arithmetic when you can.  When you cannot use exact
    arithmetic, make sure that you know the limits of the precision in your
    estimates.

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

def cfrac(x:Real, max_terms:int=50) -> list:
    """Sylvester's representation

    Returns a list [u, a, b, ..., n] of integers.  The first term is the
    greatest integer in (or floor of) x.  The remaining terms are the
    denominators in the continued fraction representation.

    The result satisfies the following inequality:

        | [u;a,b,...,n] - x| <= epsilon

    BUGS
        The result is subject to the typical disclaimers involving floating
        point rounding errors.

    SEE ALSO:
        sylvester.py - for Egyption fraction estimates
    """
    result = list()
    u = floor(x)
    result.append(u)
    r = x - u
    while r > 0 and len(result) < max_terms:
        r = 1/r
        q = floor(r)
        # print(result, q, r)
        result.append(q)
        r = r - q
    return result

def value_of(result:list) -> Fraction:
    """return the value of the decomposition"""
    a0, b0 = 0, 1
    a1, b1 = 1, 0
    for q in result:
        a = q*a1 + a0
        b = q*b1 + b0
        a0, b0 = a1, b1
        a1, b1 = a, b
    return Fraction(a, b)

def convergents_of(result:list) -> list:
    """return the convergents of the decomposition"""
    a0, b0 = 0, 1
    a1, b1 = 1, 0
    convergents = list()
    for q in result:
        a = q*a1 + a0
        b = q*b1 + b0
        convergents.append(Fraction(a, b))
        a0, b0 = a1, b1
        a1, b1 = a, b
    return convergents

def self_test():
    """a simple self test"""
    from math import sqrt, pi, e
    from os.path import basename

    def display(x, result, note):
        """display a result"""
        print("-"*72)
        xstar = value_of(result)
        print(x, xstar, result, note)
        convergents = convergents_of(result)
        print("  convergents:", *convergents)
        err = x - xstar
        print("        error:", err)

    print("Test cases for", basename(__file__))
        # 1870/741 = 2*5*11*17/3*13*19
    x = Fraction(1870, 741)
    result = cfrac(x)
    display(x, result, "2*5*11*17/3*13*19")
        # 355/113 (a very nice rational approximation of π)
    x = Fraction(355, 113)
    result = cfrac(x)
    display(x, result, "≈ π")
    print("\tπ - 355/113 ≈", pi - 355/113)
        # √2
    x = sqrt(2)
    result = cfrac(x, 10)
    display(x, result, "≈ √2")
        # (1+√5)/2
    x = (1+sqrt(5))/2
    result = cfrac(x, 12)
    display(x, result, "≈ (1+√5)/2")
        # ∛2
    x = 2**(1/3)
    result = cfrac(x, 10)
    display(x, result, "≈ ∛2")
        # e
    x = e
    result = cfrac(x, 12)
    display(x, result, "≈ e")
        # π
    x = pi
    result = cfrac(x, 10)
    display(x, result, "≈ π")
    print("Consider some subsequences and look at improvements...")
    print("Where is improvement dramatic? Where is it inconsequential?")
    display(x, result[:2], "≈ π")
    display(x, result[:3], "≈ π")
    display(x, result[:4], "≈ π")
    display(x, result[:5], "≈ π")
    display(x, result[:6], "≈ π")
    display(x, result[:7], "≈ π")

if __name__ == "__main__":
    self_test()
