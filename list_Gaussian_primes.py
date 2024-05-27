"""
list_Gaussian_primes.py - list primes up to some norm
Eric Conrad

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
from math import sqrt, isqrt      # isqrt (int sqrt) requires Python >= 3.8

from primality import primes
from gaussian_int import GaussianInt, normsq

def make_comma(csv=False):
    """build a comma"""
    return "," if csv else "\t"

def make_heading(csv=False, quadrant1=False, primes_only=False):
    """prepare the heading"""
    s = ""
    comma = make_comma(csv)
    if not primes_only:
        s += "norm" + comma
    s += "I"
    if not quadrant1:
        s += comma + "II"
        s += comma + "III"
        s += comma + "IV"
    return s

def make_detail(p, csv=False, quadrant1=False, primes_only=False,
                round_to=3):
    """prepare the line for Gaussian prime p"""

    s = ""
    comma = make_comma(csv)
    if not primes_only:                 # |p| is an integer
        p_dot_p = normsq(p)
        n = isqrt(p_dot_p)
        if n * n == p_dot_p:
            s += str(n) + comma
        else:                           # |p| is not an integer
            fmt = "%%.%df%s" % (round_to, comma)
            s += fmt % sqrt(p_dot_p)
    s += str(p)
    if not quadrant1:
        s += comma + str(p * GaussianInt(0,-1)) # second quadrant
        s += comma + str(-p)                    # third quadrant
        s += comma + str(p * GaussianInt(0,1))  # fourth quadrant
    return s

def build_list(maximum, csv=False, unsorted=False, quadrant1=False,
               primes_only=False, round_to=3, to_string=False):
    """
    build the requested list of primes

    REQUIRED ARGUMENTS

        maximum - an upper bound on the complex norms

    KEYWORD ARGUMENTS

        csv         if True, the output will be in comma-separated variable
                format.

        unsorted    if True, the output will not be sorted.

        quadrant1   if True, the second, third and fourth quadrant primes
                the negative primes, and the pure imaginary primes will
                not be listed

        primes_only if True, the first column, i.e. the complex norm,
                will not be listed.

        round_to    the number of places to the right of the decimal
                point to round the norm (default: 3).  If the default
                is taken, then rounding is to the nearest thousandth.

                    *** PROGRAMMING NOTE ***
            The proper working of function main depends on the order
            of the above arguments.  Except for round_to, the default
            values are all False.

        to_string   if True, output will be return as a string result.
                If False (default), output will be displayed to standard
                output.
    """
    items = []
    items.append((0, 0, 0))

    n = 1
    maxsq = maximum ** 2
    for b in range(maximum):                
        for a in range(1, maximum+1):       # real part > 0
            c = GaussianInt(a, b)
            enorm = normsq(c)
            if enorm > maxsq:
                break
            if c.is_prime:
                items.append((enorm, n, c))
                n += 1

    if not unsorted:
        items.sort()

    lines = make_heading(csv=csv, quadrant1=quadrant1,
                         primes_only=primes_only)
    if not to_string:
        print(lines)
        lines = ""

    for item in items:
        enorm, _, p = item          # unpack
        lines += make_detail(p, csv=csv, quadrant1=quadrant1,
                             primes_only=primes_only, round_to=round_to)
        if not to_string:
            print(lines)
            lines = ""
    return lines

def main(argv, description=None, epilogue=None):
    """main entry point"""
    import argparse

    if not description:
        description = "list Gaussian primes up to some norm"
    if not epilogue:
        epilogue = "First quadrant primes are listed, in ascending" \
            + " norm order.  Associates in other quadrants are listed" \
            + " in the same row.  The list always starts with 0.  " \
            + "The first entry in each row is the complex norm."

    parser = argparse.ArgumentParser(description=description, epilog=epilogue)
    parser.add_argument("max", type=int, \
        help="An upper bound on the complex norm of the primes")
    parser.add_argument("-c", "--csv", action="store_true", \
        help="output the list as a comma-separated variable file")
    parser.add_argument("-u", "--unsorted", action="store_true", \
        help="do not sort the output")
    parser.add_argument("-I", "--quadrant1", action="store_true", \
        help="list only the first quadrant values")
    parser.add_argument("-P", "--primes-only", action="store_true", \
        help="do not display the norm")
    parser.add_argument("-r", "--round_to", type=int, default=3, \
        help="the number of places to round the norm when it is not" \
        + " an integer.")
    args = parser.parse_args(argv)

    if args.max < 0:
        raise ValueError

    build_list(args.max, args.csv, args.unsorted,
               args.quadrant1, args.primes_only, args.round_to)

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])

