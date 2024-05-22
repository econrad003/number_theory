"""
primes.py - simple prime factorization
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

    # We maintain a list of small primes which can (and will) be extended
    # using the sieve of Eratosthenes (also known as the prime number sieve).
    # Initially the list includes all single-digit positive primes.

from numbers import Real
from math import sqrt, log, gcd

    # We define a "maybe" value for situations where testing is incomplete.
    # Programmers should take care in its use as "Maybe", used by itself as
    # a condition, will be treated as true.

Maybe = float("inf")            # True, False, or Maybe
phi = (1+sqrt(5))/2             # mean-extreme ratio (≈1.618)

class Primes(object):
    """
    manage a list of known "small" primes

    By "small", we mean that this list includes all positive primes up to
    some given value.  The list can be extended using the sieve of
    Eratosthenes (aka: the prime number sieve).

    Three class variables are used to maintain the list, namely:

        odd_primes_list - the ascending list of small odd positive primes;
        small_primes_set - the set of small positive primes;
        last_tested - the last value tested.

    The last tested value will always be at least as large as the largest
    prime in the list.  It will also be odd, since integer primes other
    than 0, 2 and -2 are necessarily odd.

    We use the ring ideal definition of prime:

        an integer p is prime if and only if p is not a unit and, for all
        integers a and b, if p divides ab, then either p divides a or p
        divides b.

    Because 0 is not a unit and 0 satisfies the zero-product property,
    this definition implies that 0 should be treated as prime.  (Note
    that 0 is not, however, irreducible.)
    """
    odd_primes_list = [3, 5, 7]
    small_primes_set = {2, 3, 5, 7}
    last_tested = 9

    @classmethod
    def sieve(cls, stop:int):
        """
        continue sieving for primes

        ARGUMENTS
            stop - a stop value for sieving

        DESCRIPTION

            Sieving continues until the last sieved value exceeds the stop
            value.
        """
        if not stop:
            return                          # nothing to do if stop is False

        if not isinstance(stop, Real):
            raise TypeError("The stop value must be a real number")
        stop = abs(stop)

            # the sieve of Eratosthenes

        while stop > cls.last_tested:
            cls.last_tested += 2            # odd values only!
            is_prime = True                 # assume prime
            for p in cls.odd_primes_list:
                if cls.last_tested % p == 0:
                    is_prime = False            # composite
                    break
            if is_prime:
                cls.odd_primes_list.append(cls.last_tested)
                cls.small_primes_set.add(cls.last_tested)

    def __init__(self, sieve_to=1000):
        """
        constructor

        ARGUMENTS

            sieve_to - a stop value for sieving
        """
        self.sieve(sieve_to)

    def is_irreducible(self, n:int, sieve=True) -> (bool, float):
        """
        test whether n is irreducible
        """
        if not isinstance(n, int):
            raise TypeError
        n = abs(n)              # ignore the sign

            # if we are allowed to sieve, the answer is simple

        if sieve:
            self.sieve(n)
            return n in self.small_primes_set

            # if we are not allowed to sieve, things are more
            # complicated.

        if n in self.small_primes_set:
            return True

        if n <= self.last_tested:
            return False

        if n % 2 == 0:
            return False            # even and absolutely larger than 2

        for p in self.odd_primes_list:
            if n % p == 0:
                return False

            # if the square root of n is smaller than the largest tested value,
            # then if no tested prime divides n, then n must be prime.

        if n < self.last_tested * self.last_tested:
            return True

        return Maybe

    def is_prime(self, n:int, sieve=True) -> (bool, float):
        """
        test whether n is prime

        The only difference between this test and the test for irrecucibility
        is when 0 is tested.  0 is prime (under ourdefinition) but not
        irreducible.

        Note that Python treats the boolean False as 0.
        """
        if n == 0:
            return True

        return self.is_irreducible(n, sieve)

    def is_unit(self, n) -> bool:
        """
        test whether n is a unit in the ordinary integers

        Returns True if n is 1 or -1, and False otherwise.

        Note that Python treats the boolean True as 1.
        """
        return n in {-1, 1}

    def is_composite(self, n:int, sieve=True) -> (bool, float):
        """
        test whether n is composite

        For the ordinary integers, we have:

            0 - prime, but not irreducible
            -1, 1 - units
            2, 3, 5, 7, 11, 13, ... - prime
            -2, -3, -5, -7, -11, -13, ... - prime

        Every other ordinary integer is composite.
        """
        if self.is_unit(n):
            return False
        f = {True:False, False:True, Maybe:Maybe}
        return f[self.is_prime(n, sieve)]

    @staticmethod
    def extract_power(n, p) -> tuple:
        """returns the exponent of p in n and the quotient

        RETURN VALUE
            (e, q) where p^e divides n and q = n / p^e.  The value
            of e is the largest e such that p^e divides n
        """
        e = 0
        while n % p == 0:
            e += 1
            n //= p
        return (e, n)

    def factor_slow(self, n) -> (list, int):
        """
        returns the unique factorization of n

        RETURN VALUE
            A list of prime power factors in ascending order. The
            list starts with 1 if n is positive and -1 if n is
            negative.  The remaining entries in the list are pairs
            of the form [p, e] where p is a positive prime, e is a
            positive integer and p^e divides n.  There is one such
            entry for each positive prime divisor of n.  If n is zero,
            the value 0 is returned.
        """
        if not isinstance(n, int):
            raise TypeError
        if n == 0:
            return 0            # zero

            # extract the unit
        unit = 1 if n>0 else -1
        factors = [unit]
        n = abs(n)
        if n == 1:
            return factors

            # extract powers of two
        e, n = self.extract_power(n, 2)
        if e > 0:
            factors.append([2,e])
        if n == 1:
            return factors          # complete

            # extract powers of small odd primes
            # Since most exponents will be 0, this is slow
        for p in self.odd_primes_list:
            e, n = self.extract_power(n, p)
            if e > 0:
                factors.append([p,e])
            if n == 1:
                return factors          # complete

            # at this point we will need to find more factors...
            # the largest possible proper factor is the square root
            # of the current value of n, but we might not need to
            # go that far. Instead, we will try to slowly increase the
            # size of the largest tested value.
            # This is the slowest part of the factorization process.

        m = len(self.odd_primes_list)
        while n > 1:
            stop = int(self.last_tested * phi)
            self.sieve(stop)                # extend the small primes
            while m < len(self.odd_primes_list):
                p = self.odd_primes_list[m]
                m += 1
                e, n = self.extract_power(n, p)
                if e > 0:
                    factors.append([p,e])
                if n == 1:
                    break

        return factors

    @staticmethod
    def multiply(factors:(list,int)) -> int:
        """multiply a list of factors"""
        if factors == 0:
            return 0
        if not isinstance(factors, list):
            raise TypeError("Require factor list")
        n = 1
        for factor in factors:
            if isinstance(factor, list):
                p, e = factor
                if not isinstance(e, int):
                    raise TypeError("Require integer exponent")
                if e < 0:
                    raise ValueError("Exponent may not be negative")
                for _ in range(e):
                    n *= p
            else:
                n *= factor
        return n

    def d(self, n:int) -> int:
        """the number of positive divisors of n

        This is the same as $\sigma_0(n)$, the sum of positive divisors raised
        the zeroth power.  It is a multiplicative number-theoretic function.

        For n=0, the value $\infty$ is returned since every positive integer
        is a divisor of zero.

        If n is a unit (i.e. n is 1 or -1), then 1 is returned.
        """
        if not isinstance(n, int):
            raise TypeError

        if n == 0:
            return float('inf')

        n = abs(n)          # the result is unsigned
        if n == 1:
            return 1        # units 1 and -1 each have one positive divisor

        factors = self.factor_slow(n)
        divisors = 1
        for factor in factors:
                # here we use the multiplicative property.
                #   (1) the set of divisors of p^e, where p is prime, is
                #       {p^0, p^1, p^2, ..., p^e}, so the number of divisors
                #       is e+1; and
                #   (2) if m and n are relatively prime, then:
                #           d(mn) = d(m) d(n)
            if isinstance(factor, list):
                p, e = factor       # we use the fact that p is prime
                divisors *= e+1     # (1) e+1 divisors of p^e; (2) multiply
        return divisors

    def sigma(self, n:int, k=1) -> (int, float, complex):
        """
        the sum of k^th powers of the positive divisors of n (but see notes!)

        BASE CASES

            For n=0, the result is undefined (float('nan').

            For n=1, the result is 1.

            For n=-1 (the other unit), we return -1 if k is odd, 1 if k is even,
            and a complex principal value for (-1)^k if  k is not an integer.

            For n=p^e (where p is a positive prime), we have:
                1 + p^k + p^(2k) + ... p^(ek)

        DESCRIPTION

            This is a multiplicative number-theoretic function.  If k is a
            non-negative integer, the result will be an integer.  If k is a
            negative integer, the result will be a floating point value.
            If k is not an integer, then the result will be a complex value.
            In this last case, the Pythonic choice of principal value will
            govern the result.
        """
        if not isinstance(n, int):
            raise TypeError

            # Base cases
        if k == 0:
            return self.d(n)

        if n == 0:
            return float('nan')

        s = 1 if n>0 else (-1)**k       # the unit result
        n = abs(n)
        if n == 1:                      # a unit was entered
            return s

        factors = self.factor_slow(n)
        for factor in factors:
            t, a = 0, 1             # t is the sum, a is the next term
            if isinstance(factor, list):
                p, e = factor       # we use the fact that p is prime
                m = p**k            # geometric sequence multiplier
                for f in range(e+1):
                    t += a          # add p^f
                    a *= m          # now a=p^(f+1)
                s *= t              # multiplicative property
        return s

    def s(self, n:int) -> int:
        """
        for a positive integer, the sum of its proper divisors

        The result is meaningless if n is not a positive integer.
        (Well, not really.  See the test cases below for hints.)
        """
        return self.sigma(n) - n

primes = Primes(20)             # an instance

def multiplicative(minus1=-1, zero=0):
    """
    decorator to define a multiplicative function

    USAGE

        @multiplicative(**kwargs)
        def f(p:int, e:int) -> int:
            '''docstring for f'''
            return the value of f at p^e

    OPTIONAL ARGUMENTS

        minus1  the value of the result when the input is -1 -- the legal
                values are -1, 0, or 1. The default is -1.

        zero    the value of the result when the input is 0 -- the legal
                values are 0, float('inf'), and float('nan'). The default
                is 0.
    """
        # validate the decorator arguments
    if not minus1 in {-1, 0, 1}:
        raise ValueError("The value at -1 should be -1, 0 or 1")

    if not zero in {0, float('inf'), float('nan')}:
        raise ValueError("The value at 0 should be zero, infinite," \
            + " or undefined")

        # define the wrapper
    def wrap(f):
        """wrapper"""

        docstring = f.__doc__[:]

        def wrapped_f(n):
            """wrapped multiplicative function

            Documentation for the wrapped function, if any, is saved
            in the attribute doc.  You can display this using the
            help atribute.

            For example, for the multiplicative function square_free
            (which returns 1 if n is square_free and 0 if there
            is a prime p such that p^2 is a divisor of n):
                (1) you can access the docstring as:
                        square_free.doc
                (2) you can display it in a more user-friendly way as:
                        square_free.help()
                (3) given an integer n, the wrapped function is
                    evaluated as:
                        square_free(n)
            """

            if n == 0:
                return zero
            retval = 1 if n>0 else minus1
            n = abs(n)
            factors = primes.factor_slow(n)
            for factor in factors:
                if isinstance(factor, list):    # factor = [p, e]
                    retval *= f(*factor)
            return retval

        wrapped_f.doc = docstring
        def helper():
            """display the docstring"""
            print(docstring)

        wrapped_f.help = helper
        return wrapped_f

    return wrap

@multiplicative(1, 0)
def square_free(p:int, e:int) -> int:
    """
    test whether the input is squarefree

    DESCRIPTION

        Undecorated, this function returns 1 if e<2 and 0 otherwise.

        After decoration with the multiplicative(1, 0) decorator, this
        becomes a multiplicative arithmetic function which returns 1
        if its integer argument is squarefree or 0 otherwise.

        This is also known as the Moebius function.

    USAGE

        primality.square_free(n)

    ARGUMENTS

        n   an integer

    RETURNS

        1 if n is squarefree and 0 if there is a positive ptime p such
        that p^2 divides n.

    REMARKS ON WRAPPING

        Since 4 divides 0, 0 is not squarefree, hence we set zero=0.
        A negative integer is squarefree if and only if its absolute
        value is squarefree.  Hence we set minus1=1.  So squarefree
        as defined here is an even multiplicative function which vanishes
        at 0.

        Programmers who (like me) can be excessively pedantic might want
        to treat this as a True-False function (or boolean).  For example,
        in a program which deals with quadratic field extensions where we
        frequently require that a given integer be squarefree, we might
        be pedantic and say:

            if not bool(primality.square_free(n)):
                raise ValueError(f"n={n} is not squarefree")

        or we might stick with integer functions and relations and say:

            if primality.square_free(n) != 1:
                raise ValueError(f"n={n} is not squarefree")

        or we might be more relaxed and simply say:

            if not primality.square_free(n):
                raise ValueError(f"n={n} is not squarefree")

        Although there are some technical differences between the three
        statements, the difference in practice is largely a matter of style.

    MOEBIUS INVERSION

        This arithmetic function is important in the theory of arithmetic
        functions on the positive integers.  Search the literature for
        connections with the Moebius function and convolution.
    """
        # this code is used inside the wrapper to handle
        # prime powers
    return 1 if e<2 else 0

@multiplicative(-1, 0)
def totient(p:int, e:int) -> int:
    """
    Euler's totient function

    DESCRIPTION

        For a positive integer n, this is the number of positive integers
        which are less than or equal to n and relatively prime to n.

        For n=1, since gcd(1,1)=1, the totient of 1 is 1.  For all other
        positive integers n, the totient of n is less than n.

        Now consider the prime p=3 with some small exponents e:
            (e=0) the totient of 3^0=1 is 1;
            (e=1) then p^e=3 -- the totient is 2 since 1 and 2 are coprime but
                3 is not;
            (e=2) then p^e=9 -- multiples of 3, namely 3, 6 and 9 are not
                coprime, so the totient is 9-3 or 6;
            (e=3) then p^e=27 -- there are 9 multiples of 3, i.e. 3, 6, 9, 12,
                15, 18, 21, 24 and 27, giving a totient of 27-9=18
        For n=3^e when e is at least 1, there are n/3 multiples of 3 leaving a
        totient of n-n/3.

        This generalizes nicely to other primes p:

            for n=p^e when e is at least 1, there are n/p multiples of p
            leaving a totient of n-n/p.

    EXTENSION

        For negative integers and zero, there is some freedom in extending the
        the totient.  The choices are not completely arbitrary, for example,
        we cannot have totient(-1)=5 as this would not preserve
        multiplicativity.  And even among the available choices, we should
        have at least some guiding principles.

        For n=0, we will take the value as 0 as there are no positive integers
        less than or equal to 0.

        For n<0, we will take totient(n)=-totient(n), making totient an odd
        multiplicative function.

        The wrapper is multiplicative(-1, 0).
    """
        # note that, when wrapped, e is a positive integer
    n = p**e
    return n - n // p

if __name__ == "__main__":
    print(f"Testing {__file__}.")
    # print("primes =", sorted(list(primes.small_primes_set)))
    assert primes.odd_primes_list == [3, 5, 7, 11, 13, 17, 19]
    assert primes.small_primes_set == {2, 3, 5, 7, 11, 13, 17, 19}
    assert primes.last_tested == 21

    assert primes.is_unit(1) == True
    assert primes.is_unit(-1) == True
    assert primes.is_unit(0) == False

    assert primes.is_irreducible(0) == False
    assert primes.is_prime(0) == True
    assert primes.is_irreducible(1) == False
    assert primes.is_irreducible(2) == True
    assert primes.is_irreducible(-2) == True
    assert primes.is_irreducible(7) == True
    assert primes.is_irreducible(-7) == True
    assert primes.is_irreducible(9) == False
    assert primes.is_irreducible(-9) == False
    assert primes.is_irreducible(18, False) == False
    assert primes.is_irreducible(-18, False) == False
    assert primes.is_irreducible(19, False) == True
    assert primes.is_irreducible(-19, False) == True
    assert primes.is_irreducible(22, False) == False
    assert primes.is_irreducible(-22, False) == False
    assert primes.is_irreducible(25, False) == False
    assert primes.is_irreducible(-25, False) == False
    assert primes.is_irreducible(23, False) == True
    assert primes.is_irreducible(-23, False) == True
    assert primes.is_irreducible(23*23, False) == Maybe
    assert primes.is_irreducible(29) == True
    assert primes.is_irreducible(-29, False) == True
    final = sorted(list(primes.small_primes_set))
    # print("primes =", final)
    assert final == [2] + primes.odd_primes_list

        # test factorization
    factors = primes.factor_slow(360)
    # print(360, "factors =", factors)
    assert primes.multiply(factors) == 360
    n = 997*997*360
    factors = primes.factor_slow(n)
    # print(n, "factors =", factors)
    assert primes.multiply(factors) == n

    # final = sorted(list(primes.small_primes_set))
    # print("primes =", final)

        # test divisor sums
    assert primes.d(6) == 4         # 1, 2, 3, 6
    assert primes.d(360) == 24      # 360=2^3 3^2 5; 4x3x2=24

        # perfect numbers
    assert primes.s(6) == 6         # 1 + 2 + 3 = 6
    assert primes.s(28) == 28       # 1 + 2 + 4 + 7 + 14 = 28
        # deficient numbers
    assert primes.s(1) == 0         # no proper divisors so s(1)=0
    assert primes.s(3) == 1         # s(p)=1 iff p is prime
    assert primes.s(9) == 4         # 1 + 3 = 4
    assert primes.s(10) == 8        # 1 + 2 + 5 = 8
        # abundant numbers
    assert primes.s(12) == 16       # 1 + 2 + 3 + 4 + 6 = 16
    assert primes.s(30) == 42       # 1 + 2 + 3 + 5 + 6 + 10 + 15 = 42
        # amicable pairs
        #   220, 284
    assert primes.s(220) == 284
    assert primes.s(284) == 220
        #   1184, 1210
    assert primes.s(1184) == 1210
    assert primes.s(1210) == 1184
        # sociable numbers
    a = 1264460
    msg = f"{a} is a sociable number with an aliquot sequence of period 4."
    print("Claim:", msg)
    n = a
    for _ in range(4):
        factors = primes.factor_slow(n)
        print(f"factors of {n}:", factors)
        m = primes.s(n)
        print(f"s({n}) = {m}")
        n = m
    assert n == a

        # the undocumented stuff
    assert primes.s(-1) == 0        # -1 - (-1) = 0
    assert primes.s(-3) == -1       # -1 (1 + 3) - (-3) = -4 + 3 = -1
    assert primes.s(-6) == -6       # -1 (1+2) (1+3) - (-6) = -12 + 6 = -6

        # number-theoretic functions
    print("Testing implementation of square_free...")
    test_values = [1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21, 22, 23]
    for n in test_values:
        assert square_free(n), f"ERROR: {n} is squarefree."
        assert square_free(-n), f"ERROR: {-n} is squarefree."
    test_values = [4, 8, 9, 12, 16, 18, 20, 24, 25]
    for n in test_values:
        assert not square_free(n), f"ERROR: {n} is not squarefree."
        assert not square_free(-n), f"ERROR: {-n} is not squarefree."
    assert not square_free(0)

    print("Testing implementation of Euler's totient...")
    test_values = [[1,1], [2,1], [3,2], [4,2], [5,4], [6,2], [7,6], [8,4],
                   [9,6], [10,4], [11,10], [12,4], [13,12], [14,6], [15,8]]
    for pair in test_values:
        x, y = pair
        assert totient(x) == y, f"EXPECTED: φ({x})={y}, GOT {totient(x)}"
        assert totient(-x) == -y, f"EXPECTED: φ({-x})={-y}, GOT {totient(-x)}"
    assert totient(0) == 0

        # brute force calculation of φ(360)
    m = 0                   # counter
    for k in range(1,361):
        if gcd(k, 360) == 1:
            m += 1              # k and 360 are relatively prime
        # test
    assert totient(360) == m, f"EXPECTED {m}"
    assert totient(-360) == -m, f"EXPECTED {-m}"

        # usage statistics
    print()
    print("Table Space Statistics:")
    n = primes.last_tested
    print(f"\tn={n} largest number tested")
    print(f"\tlog(n) = {log(n)}")
    print(f"\tn / log(n) = {n/log(n)} asymptotic prime density")
    print(f"\t{len(primes.odd_primes_list)+1} positive prime numbers found")

    print("SUCCESS!")