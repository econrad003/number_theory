"""
gaussian_int.py - Gaussian integers and rationals
Eric Conrad

DEFINITION

    A Gaussian integer is a complex number a+bi whose real part a and
    whose imaginary part b are both integers.  The set of Gaussian
    integers is usually called Q[i] ("queue bracket eye") with the Q
    written either in bold face or in blackboard bold.

    A Gaussian rational number is the quotient u/v of two Gaussian integers
    u and v, with v not equal to zero.  The set of Gaussian rationals
    is usually called Q(i) ("queue paren eye"), again with Q in either
    bold or blackboard bold.

CONVENTIONS

    We will try to observe the following conventions in this module's
    documentation:
        (1) Majescules denote sets

            C - the complex numbers
            N - the non-negative integers (aka: natural numbers)
            Q[i] - the Gaussian integers
            Q(i) - the Gaussian rationals
            Q - the ordinary rationals
            R - the real numbers
            Z - the ordinary integers

            Note that here we count zero as a natural number.  Also, as
            is customary among matheticians, we use the minuscule i to
            denote one of the square roots of -1.  In a graph of the Argand
            plane (aka: Gaussian plane), the value i i plotted 1 unit above
            the value 0, while the other square root (-i) is plotted 1 unit
            below.  In a right-handed coordinate system, 1 is plotted 1 unit
            to the right of zero, and -1 to the left:

                                 |
                                 2i                      X
                                 |
                                 |
                                 i
                                 |
                                 |
                    --- -1 ----- 0 ----- 1 ----- 2 ----- 3 ---
                                 |
                                 |
                                 -i
                                 |
                                 |

                Figure 1: Argand plane (right-handed coordinates)
                    The point marked X represents the value 3+2i.

        (2) Minuscules

            i - the square root of -1 in the upper half-plane
            a, b, c, ... - real numbers, typically integers and rationals
            u, v, w, ... - complex numbers, typically Gaussian integers
                and rationals

DESCRIPTION

    Two classes are defined here:

        GaussianInt - representing Gaussian integers
        GaussianFrac - representing Gaussian rationals

    These are both subclasses of complex numbers.

    If GaussianInt is used for multiplicative number theory, then the
    module primality.py is required.

    Python 3.8 or greater is required for math.isqrt().

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

from math import sqrt, gcd, atan, pi, isnan, isqrt
    # NOTE: math.isqrt() requires Python 3.8 or later

from numbers import Complex
from fractions import Fraction

    # NOTE: primality.py is imported by class GaussianInt if it is
    #   needed.

def sgn(x):
    """the signum function

    sgn(x) is 1 if x>0, 0 if x=0, and -1 if x<0
    """
    if x>0:
        return 1
    return -1 if x<0 else 0

class GaussianFrac(Complex):
    """the class Q(i) of Gaussian rational numbers

    Q(i) = Q+Qi = {a+bi: a in Q, b in Q}.
    """

    def __init__(self, a:int=0, b:int=0, c:int=1):
        """constructor

        u = (a+bi)/c
        """
        for n in (a, b, c):
            if not isinstance(n, int):
                raise TypeError("The values a and b must be integers.")
        if c == 0:
            raise ZeroDivisionError("The value of c cannot be zero.")

            # normalization
        d = gcd(a, b, c)        # okay since c != 0
        if sgn(c) != sgn(d):
            d = -d
        a //= d
        b //= d
        c //= d         # c>0

        self._a = a
        self._b = b
        self._c = c

    def __repr__(self):
        """representation"""
        if self._c == 1:          # a Gaussian integer
            if self._a == 0:
                if self._b == 0:
                    return "0"
                return f"{self._b}i"
            else:
                if self._b > 0:
                    return f"{self._a}+{self._b}i"
                if self._b < 0:
                    return f"{self._a}-{-self._b}i"
                return f"{self._a}"
        if self._b == 0:
            return f"{self._a}/{self._c}"
        if self._a == 0:
            return f"{self._b}i/{self._c}"
        if self._b > 0:
            return f"({self._a}+{self._b}i)/{self._c}"
        return f"({self._a}-{-self._b}i)/{self._c}"

    def __hash__(self):
        """return the hash"""
        if self._b == 0:            # imaginary part is zero
            if self._c == 1:
                return hash(self._a)
            return hash(Fraction(self._a, self._c))
        return hash((self._a, self._b, self._c))

    @property
    def numerator(self):
        """the numerator

        This will allways be a Gaussian integer...  We return an
        int if the imaginary part is zero.
        """
        if self._b == 0:
            return self._a
        return GaussianInt(self._a, self._b)

    @property
    def denominator(self) -> int:
        """the denominator"""
        return self._c

    @property
    def real(self) -> (Fraction, int):
        """return the real part"""
        result = Fraction(self._a, self._c)
        if result.denominator == 1:
            return result.numerator
        return result

    @property
    def imag(self) -> (Fraction, int):
        """return the imaginary part"""
        result = Fraction(self._b, self._c)
        if result.denominator == 1:
            return result.numerator
        return result

    def normsq(self) -> (Fraction, int):
        """square of the Euclidean norm"""
        result = Fraction(self._a ** 2 + self._b ** 2, self._c ** 2)
        if result.denominator == 1:
            return result.numerator
        return result

    def __abs__(self):
        """complex absolute value (or 2-norm or Euclidean norm)"""
        return sqrt(self.normsq())

    def arg(self) -> float:
        """angle of inclination in radians

        The result is in the interval [0, 2*pi) and is subject to
        rounding errors and difficulties for points near the imaginary
        axis (a+bi with |a| >> |b|).
        """
        if self.is_zero:
            return float('nan')         # undefined at zero
        if self._a == 0:                # vertical tangent
            return pi/2 if self._b > 0 else 3*pi/2

        alpha = atan(abs(self._b / self._a))    # reference angle
        if self._b >= 0:                # upper half plane
            if self._a > 0:                 # first quadrant
                return alpha
            else:                           # second quadant
                return pi - alpha
        # else:                           # lower half plane
        if self._a < 0:                     # third quadrant
            return pi + alpha
        #   else:                           # fourth quadrant
        return 2*pi - alpha

    def argd(self) -> float:
        """angle of inclination in degrees

        Problems similar to those with arg can occur.  Exact values
        are returned for a+bi whenever either a=0 or b=0 or |a|=|b|.
        The returned value is in the interval [0,360).
        """
            # We have a small number of exact solutions
        if self.is_zero:
            return float('nan')         # undefined at zero
        if self._b == 0:                # horizontal tangent
            return 0 if self._a > 0 else 180
        if self._a == 0:                # vertical tangent
            return 90 if self._b > 0 else 270
        if abs(self._a) == abs(self._b):    # 45 degree solutions
            if self._b > 0:             # upper half plane
                if self._a > 0:             # first quadrant
                    return 45
                else:                       # second quadrant
                    return 135
            else:                       # lower half plane
                if self._a < 0:             # third quadrant
                    return 225
                else:                       # fourth quadrant
                    return 315

        return self.arg() * 180 / pi

    def conjugate(self):
        """complex conjugate"""
        return GaussianFrac(self._a, - self._b, self._c)

    def __complex__(self):
        """cast to complex"""
        return (self._a + self._b * 1j) / self._c

    def __eq__(self, other):
        """test equality"""
        if isinstance(other, GaussianFrac):
            return self._a == other._a and self._b == other._b \
                and self._c == other._c
        return complex(self) == complex(other)

    @property
    def is_zero(self):
        """test if the value is zero (the additive identity)"""
        return self._a == 0 and self._b == 0

    @property
    def is_one(self):
        """test if the value is one (the multiplicative identity)"""
        return self._a == 1 and self._b == 0 and self._c == 1

    def __add__(self, other):
        """add a value to a Gaussian integer"""
        if self.is_zero:
            return other
        if isinstance(other, GaussianFrac):
            a = self._a * other._c + self._c * other._a
            b = self._b * other._c + self._c * other._b
            c = self._c * other._c
            result = GaussianFrac(a, b, c)
            if result.denominator == 1:
                return result.numerator
            return result
        if isinstance(other, int):
            other = GaussianFrac(other)
            return self + other
        if isinstance(other, Fraction):
            other = GaussianFrac(other.numerator, 0, other.denominator)
            return self + other
        if isinstance(other, (float, complex)):
            a, b = other.real, other.imag
            if a == int(a) and b == int(b):
                return self + GaussianInt(int(a), int(b))
            return complex(self) + other
        raise NotImplementedError

    def __radd__(self, other):
        """reverse addition"""
        return self + other         # addition is commutative

    def __pos__(self):
        """pass-through unchanged"""
        return self

    def __neg__(self):
        """additive inverse"""
        return GaussianFrac(-self._a, -self._b, self._c)

    def __mul__(self, other):
        """multiplication"""
        if self.is_zero:
            return 0
        if self.is_one:
            return other
        if isinstance(other, GaussianFrac):
            a = self._a * other._a - self._b * other._b
            b = self._a * other._b + self._b * other._a
            c = self._c * other._c
            result = GaussianFrac(a, b, c)
            if result.denominator == 1:
                return result.numerator
            return result
        if isinstance(other, int):
            other = GaussianFrac(other)
            return self * other
        if isinstance(other, Fraction):
            other = GaussianFrac(other.numerator, 0, other.denominator)
            return self * other
        if isinstance(other, (float, complex)):
            a, b = other.real, other.imag
            if a == int(a) and b == int(b):
                return self * GaussianInt(int(a), int(b))
            return complex(self) * other
        raise NotImplementedError

    def __rmul__(self, other):
        """reverse multiplication"""
        return self * other         # multiplication is commutative

    @property
    def reciprocal(self):
        """the multiplicative inverse

        1 / ((a+bi)/c) = c / (a+bi)
                       = c (a-bi) / ((a+bi)(a-bi))
                       = (ac - bci) / (a^2 + b^2)
        """
        if self.is_zero:
            raise ZeroDivisionError("divide by zero")
        a, b, c = self._a, self._b, self._c
        result = GaussianFrac(a*c, -b*c, a**2 + b**2)
        if result.denominator == 1:
            return result.numerator
        return result

    def __truediv__(self, other):
        """multiplication"""
        if isinstance(other, GaussianFrac):
            return self * other.reciprocal
        if other == 1:
            return self
        if isinstance(other, int):
            return self / GaussianInt(other)
        if isinstance(other, Fraction):
            return self / GaussianFrac(other.numerator, 0, other.denominator)
        return self * (1 / other)

    def __rtruediv__(self, other):
        """reverse multiplication"""
        return other * self.reciprocal

    def integer_power(self, e:int):
        """raise to an integer power"""
        if not isinstance(e, int):
            raise TypeError("This method requires an integer")

        if e == 0:
            if self.is_zero:
                raise ZeroDivisionError("0**0 is undefined")
            return 1

        if e < 0:
            return self.reciprocal ** -e

        stack = []
        f = 1
        a, b, c = self._a, self._b, self._c
        while f <= e:
                # update the stack
            stack.append((f, GaussianFrac(a, b, c)))
            f *= 2              # double f
                # square the top of the stack
            a2 = a**2 - b**2
            b2 = 2*a*b
            c2 = c**2
                # update the current value
            a, b, c = a2, b2, c2
        result = 1
        while e > 0:
            f, u = stack.pop()
            if f <= e:          # 1 in the binary expansion of e
                result *= u
                e -= f
        return result

    def __pow__(self, other):
        """exponentiation"""
        if isinstance(other, int):
            return self.integer_power(other)

        return complex(self) ** complex(other)

    def __rpow__(self, other):
        """reverse exponentiation"""
        if self._b == 0:
            if self._c == 1:
                return other ** self._a
            return other ** Fraction(self._a, self._b)
        return other ** complex(self)

class GaussianInt(GaussianFrac):
    """Gaussian integers"""

    PRIMALITY = False
    _pr = None
    _primes = None

    def __init__(self, a:int=0, b:int=0):
        """constructor"""
        super().__init__(a, b)

    @classmethod
    def _set_primality(cls):
        """private class method called for multiplicative number theory"""
        if cls.PRIMALITY:
            return                  # primes has been imported
        cls.PRIMALITY = True
        import primality
        cls._pr = primality
        cls._primes = primality.primes

    @property
    def is_unit(self) -> bool:
        """
        test whether the number is a Gaussian unit

        Returns True if n is 1 or -1 or i 0r -i, and False otherwise.

        Note that Python treats the boolean True as 1.
        """
        i = GaussianInt(0, 1)
        return self in {-1, 1, i, -i}

    def is_associate(self, other) -> bool:
        """test whether the numbers are associates"""
        u = self / other
        i = GaussianInt(0, 1)
        return u in {-1, 1, i, -i}

    @property
    def is_prime(self) -> bool:
        """test whether the number is a Gaussian prime

        To show that a+bi is a Gaussian prime, there are several steps --

            1) Calculate the square of the norm:
                    N = (a+bi)(a-bi) = a^2 + b^2.

            2) If N is an ordinary prime, then a+bi is prime.
               Examples:
                    (1+i)(1-i) = 1 + 1 = 2, so 1+i and 1-i are prime;
                    (2+i)(2-i) = 4 + 1 = 5, so 2+i and 2-i are prime;
                    (3+2i)(3-2i) = 9 + 4 = 13 so 3+2i and 3-2i are prime.
               Note that if N is an ordinary prime, then N % 4 is either
               1 or 2.  Note also that 2, 5 and 13 are reducible as
               Gaussian integers and hence not Gaussian primes.  It can
               be shown that an ordinary prime p is a Gaussian prime
               if and only if p % 4 = 3.

               If N is not an ordinary prime, we must continue...

            3) If N is not a perfect square, then a+bi is not prime.
               Example:
                    (3+5i)(3-5i) = 9 + 25 = 34 = 2 * 17.
                    Above we found that (1+i) and (1-i) are prime divisors
                    of 2, so one of these must be a proper divisor of 3+5i.
                        (3+5i) / (1+i) = (3+5i)(1-i) / 2
                                       = (3+5i-3i+5) / 2
                                       = (8+2i)/2 = 4+i
                    So:
                        3+5i = (i+i)(4+i)
                    Note also that 17=(4+i)(4-i).

                If N is a perfect square, we must continue.

            4) Let n be the square root of N, i.e. the complex norm
               of a+bi.  At this point we know n is an integer.  If
               n is not prime, then a+bi is not prime.  If n % 4 = 3,
               then a+bi is prime.  Otherwise a+bi is not prime.
               Examples:
                    (Case i)
                        (13+84i)(13-84i) = 13^2 + 84^2 = 7225 = 85^2
                        Since 85 = 5x17 is an ordinary composite,
                        13+84i is a Gaussian composite.
                        In fact:
                            (13+84i)/(2+i) = (13+84i)(2-i) / 5
                                           = (110+155i) / 5 = 22+31i
                        The trial divisor 2+i came from factoring 5
                        as indicated above.

                    (Case ii)
                        (3+0i)(3-0i) = 3^2 + 0^2 = 9 = 3^2.
                        Since 3 is an ordinary prime and 3 % 4 = 3, it
                        follows that 3 is a Gaussian prime.

                    (Case iii)
                        17x17 is composite because 17 % 4 = 1.
                        (4+i)(4-i) = 17
        """
        self._set_primality()           # primality.py is required
        normsq = self.normsq()

            # STEP 1
            # if the norm-square is prime, then the number is prime
        if self._primes.is_prime(normsq):
            return True

            # STEP 2
            # we must look at the square root...
        norm = isqrt(normsq)            # integer square root (Python >= 3.8)
        if norm * norm != normsq:
            return False                # not a perfect square so not prime

            # STEP 3
            # normsq is a perfect square
            # for the number to be a Gaussian prime, the norm must be
            # BOTH congruent to three modulo 4 AND an ordinary prime
        return norm % 4 == 3 and self._primes.is_prime(norm)

if __name__ == "__main__":
        # TEST CASES

    print("Test GaussianFrac and GaussianInt:")
    
    print("\tTest __repr__")
    assert str(GaussianFrac()) == "0"
    assert str(GaussianFrac(3)) == "3"
    assert str(GaussianFrac(0, 5)) == "5i"
    assert str(GaussianFrac(7,11)) == "7+11i"
    assert str(GaussianFrac(13,-17)) == "13-17i"
    assert str(GaussianFrac(19,-23,29)) == "(19-23i)/29"

    print("\tTest __eq__")
    assert GaussianFrac() == 0
    assert GaussianFrac(2) == 2
    assert GaussianFrac(b=3) == 3j
    assert GaussianFrac(5,7) == 5+7j
    assert GaussianFrac(11,13) == GaussianInt(11,13)
    assert GaussianFrac(17,19,23) == GaussianFrac(34,38,46)
    assert GaussianFrac(29,31,37) == GaussianFrac(-29,-31,-37)

    print("\tTest normsq:")
    assert GaussianFrac().is_zero
    assert not GaussianFrac(1).is_zero
    assert GaussianFrac(3,4).normsq() == 25
    assert GaussianFrac(3,4,5).normsq() == 1

    print("\tTest arg and argd:")
    assert isnan(GaussianFrac().arg())
    assert isnan(GaussianFrac().argd())
        # some exact values
    assert GaussianFrac(1).arg() == 0
    assert GaussianFrac(0,1).arg() == pi/2
    assert GaussianFrac(-1).arg() == pi
    assert GaussianFrac(0,-1).arg() == 3*pi/2
        # some estimates
    theta = GaussianFrac(1,2).arg()     # Quadrant I above y=x
    assert pi/4 < theta and theta < pi/2, theta
    theta = GaussianFrac(-3,5).arg()    # Quadrant II above y=-x
    assert pi/2 < theta and theta < 3*pi/4, theta
    theta = GaussianFrac(-3,-5).arg()   # Quadrant III below y=x
    assert 5*pi/4 < theta and theta < 3*pi/2, theta
    theta = GaussianFrac(3,-5).arg()    # Quadrant IV below y=-x
    assert 3*pi/2 < theta and theta < 7*pi/4, theta
        # some exact values
    assert GaussianFrac(1).argd() == 0
    assert GaussianFrac(1,1).argd() == 45
    assert GaussianFrac(0,1).argd() == 90
    assert GaussianFrac(-1,1).argd() == 135
    assert GaussianFrac(-1).argd() == 180
    assert GaussianFrac(-1,-1).argd() == 225
    assert GaussianFrac(0,-1).argd() == 270
    assert GaussianFrac(1,-1).argd() == 315
        # some estimates
    theta = GaussianFrac(1,2).argd()     # Quadrant I above y=x
    assert 45 < theta and theta < 90, theta
    theta = GaussianFrac(-3,5).argd()    # Quadrant II above y=-x
    assert 90 < theta and theta < 135, theta
    theta = GaussianFrac(-3,-5).argd()   # Quadrant III below y=x
    assert 225 < theta and theta < 270, theta
    theta = GaussianFrac(3,-5).argd()    # Quadrant IV below y=-x
    assert 270 < theta and theta < 315, theta

    print("\tTest __add__, __neg__ and __radd__:")
    assert GaussianInt(2,3) + GaussianInt(5,7) == GaussianInt(7,10)
    u = GaussianFrac(11,13,17)
    v = GaussianFrac(19,-23,29)
    w1 = u + v
    w2 = GaussianFrac(642,-14,493)
    assert w1 == w2, f"{u}+{v}={w2} (got {w1})"
    w1 = u - v
    w2 = GaussianFrac(-4,768,493)
    assert w1 == w2, f"{u}-{v}={w2} (got {w1})"
    assert 5 + GaussianFrac(3,4,5) == GaussianFrac(3+25,4,5)
    assert 5j + GaussianFrac(3,4,5) == GaussianFrac(3,4+25,5)
    assert type(5j + GaussianFrac(3,4,5)) == GaussianFrac
    # print(GaussianFrac(3,4,5).real)
    # print((5.2 + GaussianFrac(3,4,5)).real)
    assert abs((5.2 + GaussianFrac(3,4,5)).real - (5.2+3/5)) < 1e-5
    assert type(5.2 + GaussianFrac(3,4,5)) == complex
    assert abs((5.2j + GaussianFrac(3,4,5)).imag - (5.2+4/5)) < 1e-5
    assert type(5.2j + GaussianFrac(3,4,5)) == complex

    print("\ttest __mul__, __rmul__, __truediv__ and __rtruediv__")
    u = GaussianFrac(3,5,7)
    v = GaussianFrac(11,13,17)
    w1 = u * v
    w2 = GaussianFrac(-32,94,119)
    assert w1 == w2, f"({u})*({v})={w2} (got {w1})"
    w1 = u / v
    w2 = GaussianFrac(833,136,1015)
    assert w1 == w2, f"({u})*({v})={w2} (got {w1})"
    assert 5*u == GaussianFrac(15,25,7)
    assert 5j*u == GaussianFrac(-25,15,7)
    assert type(5.2 * GaussianInt(2,-3)) == complex
    assert abs(5.2 * GaussianInt(2,-3) - (10.4-15.6j)) < 1e-5
    assert type(5.2j * GaussianInt(2,-3)) == complex
    assert abs(5.2j * GaussianInt(2,-3) - (15.6+10.4j)) < 1e-5
    assert type((2+2j) / GaussianInt(1,1)) == int 
    assert (2+2j) / GaussianInt(1,1) == 2 
    assert abs((2.1+2.1j) / GaussianInt(1,1) - 2.1) < 1e-5

    print("\ttest __pow__")
    assert GaussianInt(1,1)**2 == GaussianInt(0,2)  # (1+i)^2=1+2i-1=2i
    assert GaussianInt(1,1)**3 == GaussianInt(-2,2) # (1+i)^3=2i(i+1)=-2+2i
    assert GaussianInt(1,1)**4 == -4                # (1+i)^4=(2i)^2=-4
    assert GaussianInt(1,1)**7 == GaussianInt(8,-8) # (1+i)^7=(-4)(-2+2i)=8-8i
        #   2 / (1-i) = 2(1+i)/2 = 1+i
    assert GaussianFrac(1,-1,2)**-1 == GaussianInt(1,1)
    assert GaussianFrac(1,-1,2)**-7 == GaussianInt(1,1)**7
    u = GaussianInt(0,1)**GaussianInt(0,1)
    v = 1j**1j
    assert abs(u-v) < 1e-5

    print("\ttest primality interface")
    assert GaussianInt(1,1).is_prime
    assert not GaussianInt(2,0).is_prime
    assert GaussianInt(3,0).is_prime
    assert not GaussianInt(3,5).is_prime
    assert not GaussianInt(4,0).is_prime
    assert not GaussianInt(5,0).is_prime
    assert not GaussianInt(13,84).is_prime

    print("\tGaussian primes whose norms are less than or equal to 10:")
    for i in range(0, 10):
        for j in range(0, 10):
            u = GaussianInt(i, j)
            if u.normsq() > 100:     # (sic!) norm > 10 ignored (10^2=100)
                continue
            if not u.is_prime:
                continue
            if u == 0:
                print(f"\t  u={GaussianInt(i,j)} (the zero ideal is prime)")
                continue
            if u.imag == 0:
                print(f"\t  u={GaussianInt(i)}, {GaussianInt(-i)}; " \
                    + f"\t|u|={isqrt(u.normsq())}, |u|^2={u.normsq()}")
                continue
            if u.real == 0:
                print(f"\t  u={GaussianInt(i,j)}, {GaussianInt(i,-j)}; " \
                    + f"\t|u|={isqrt(u.normsq())}, |u|^2={u.normsq()}")
                continue
            print(f"\t  u={GaussianInt(i,j)}, {GaussianInt(i,-j)}, " \
                + f"{GaussianInt(-i, j)}, {GaussianInt(-i,-j)};  " \
                + f"\t|u|={sqrt(u.normsq()):.5f}, |u|^2={u.normsq()}")
            assert GaussianInt(i, -j).is_prime
            assert GaussianInt(-i, j).is_prime
            assert GaussianInt(-i, -j).is_prime
    print("\tYou should visually check that 0, 3, -3, 3i, -3i, 7, -7, 7i and")
    print("\t-7i are listed, and units 1 and i and ordinary primes 2 and 5")
    print("\tare not listed.  For all other entries, the squares of their")
    print("\tnorms should all be ordinary primes congruent to 1 modulo 4.")

        # THIS PRINTS IF ALL TESTS PASS
    print("SUCCESS!")
