"""
aliquot.py - aliquot sequences
Eric Conrad

USAGE

    python3 aliquot.py [-h] [-l LENGTH] [-L LARGEST_VALUE] [-P LARGEST_PRIME]
        [--log] [--print] n

    evaluate the aliquot sequence for an integer

    POSITIONAL ARGUMENTS

        n
            an integer whose aliquot sequence is to be evaluated

    OPTIONS

        -h, --help
            show this help message and exit
        -l LENGTH, --length LENGTH
            the number of terms to evaluate
        -L LARGEST_VALUE, --largest-value LARGEST_VALUE
            if a term exceeds this value, raise an exception
        -P LARGEST_PRIME, --largest-prime LARGEST_PRIME
            if a prime factor exceeds this value, raise an exception
        --log
            if set, the vertical scale will be logarithmic
        --print
            if set, print the output but do not produce a plot

    The default is to plot the terms of the aliquot sequence. Options
    include a logarithmic plot or simply displaying the sequence.
    For a logarithmic plot, to avoid domain errors, a value of zero is
    mapped to -1.

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

    # get the important stuff from primality.py

from math import log
from primality import primes, multiplicative, square_free, totient

class Aliquot(object):

    status = 0

        # termination constants
    PREMATURE_TERMINATION = 90
    TERMINATING = 1
    PERIODIC = 2
    LENGTH_EXCEEDED = 3
    LARGEST_EXCEEDED = 4
    SIEVE_OVERFLOW = 5

    @classmethod
    def aliquot_sequence(cls, n:int,
                         max_length=50, largest=100000,
                         largest_prime=100000, warnings=True):
        """
        return the leading terms of the aliquot sequence of n

        REQUIRED ARGUMENTS

            n - a positive integer

        OPTIONAL ARGUMENTS

            max_length (default 50)
                the maximum length of the sequence returned

            largest (default 100,000)
                the largest permissible value in the sequence

            largest_prime (default 100,000)
                an upper bound on positive primes obtained by
                sieving

            warnings (default True)
                display warning messages

        RETURN VALUE

            a list containing the head of the aliquot sequence

        SIDE EFFECTS

            the status value is updated and may be queried to
            determine the reason for termination.
 
        """
        def warning(msg):
            """display a warning message"""
            if warnings:
                print("WARNING (aliquot):", msg)

        cls.status = cls.PREMATURE_TERMINATION
        if not isinstance(n, int):
            raise TypeError
        n = abs(n)      # force this to be non-negative
        seq = []
        domain = set()
        old_pmax = primes.odd_primes_list[-1]
        while True:
            n = primes.s(n)             # sum of proper divisors
            seq.append(n)

                # Check for normal completion
            if n == 0:
                break                       # terminating sequence
            if n in domain:
                cls.status = cls.PERIODIC
                return seq                  # periodic sequence

                # Check for abnormal completion
            if len(seq) >= max_length:
                warning(f"Maximum length {max_length} exceeded.")
                cls.status = cls.LENGTH_EXCEEDED
                return seq
            if n > largest:
                warning(f"Largest value {n} > {largest}")
                cls.status = cls.LARGEST_EXCEEDED
                return seq
            pmax = primes.odd_primes_list[-1]
            if pmax > largest_prime:
                warning(f"Largest prime {pmax} > {largest_prime}")
                cls.status = cls.SIEVE_OVERFLOW
                return seq

                # additional warnings (non-fatal)
            if pmax > old_pmax:
                warning(f"sieved up to {pmax} from {old_pmax}")
            old_pmax = pmax

                # add the new value to the domain
            domain.add(n)

        cls.status = cls.TERMINATING    # sequence terminates with zero
        return seq

def pad_terminating_seq(seq, length):
    """
    pad the result with zeros

    This will modify the original sequence.
    """
    while len(seq) < length:
        seq.append(0)

def pad_periodic_seq(seq, length):
    """
    determine the period and pad the sequence periodically

    This will modify the original sequence.
    """
    last_term = seq[-1]
    index = -2
    while seq[index] != last_term:
        index -= 1
    period = abs(index) - 1
    diagnoses = {1:"perfect", 2:"amicable"}
    diagnosis = diagnoses.get(period, "sociable")
    print(f"Period {period}: repeats with {diagnosis} number {last_term}")
    while len(seq) < length:
        seq.append(seq[index])

def plot(seq, title="aliquot sequence plot", logarithmic=False):
    """plot an aliquot sequence"""
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    foolog = lambda x: x if x > 0 else 0.1

    xs = list(range(len(seq)))
    ys = list(map(foolog, seq)) if logarithmic else seq
    plt.title(title)
    plt.xlabel("index k")
    plt.ylabel("term log(a(k),10)" if logarithmic else "a(k)")
    if logarithmic:
        plt.yscale("log")
        formatter = ticker.StrMethodFormatter('{x:.0f}')
        plt.gca().yaxis.set_minor_formatter(formatter)

    plt.plot(xs, ys)
    if logarithmic:
        NOTE = "minimum value of -1 is used if a(k)=0."
        plt.figtext(0.99, 0.01, NOTE, horizontalalignment='right') 
    plt.show()

def main(argv):
    """parse arguments"""
    import argparse

    DESCRIPTION = "evaluate the aliquot sequence for an integer"
    EPILOGUE = "The default is to plot the terms of the aliquot sequence. " \
        + "Options include a logarithmic plot or simply displaying the " \
        + "sequence.  To avoid domain errors on terminating sequences in, " \
        + "logarithmic plots, we evaluate the (common) log(0) as -1."
    parser = argparse.ArgumentParser(description=DESCRIPTION, \
        epilog = EPILOGUE)
    parser.add_argument("n", type=int, \
        help = "an integer whose aliquot sequence is to be evaluated")
    parser.add_argument("-l", "--length", type=int, default=50, \
        help = "the number of terms to evaluate")
    parser.add_argument("-L", "--largest-value", type=int, default=100000, \
        help = "if a term exceeds this value, raise an exception")
    parser.add_argument("-P", "--largest-prime", type=int, default=100000, \
        help = "if a prime factor exceeds this value, raise an exception")
    parser.add_argument("--log", action="store_true", \
        help = "if set, the vertical scale will be logarithmic")
    parser.add_argument("--print", action="store_true", \
        help = "if set, print the output but do not produce a plot")
    args = parser.parse_args(argv)

    # print(args)           ## debugging
    seq = Aliquot.aliquot_sequence(args.n, max_length=args.length,
                                   largest=args.largest_value,
                                   largest_prime=args.largest_prime)

    print(seq)
    status = Aliquot.status
    if status == Aliquot.LENGTH_EXCEEDED:
        print("The sequence does not repeat or terminate after",
              args.length, "terms.")
    elif status == Aliquot.PERIODIC:
        print("The sequence is periodic.")
        pad_periodic_seq(seq, args.length)
    elif status == Aliquot.TERMINATING:
        print("The sequence is terminating.")
        pad_terminating_seq(seq, args.length)
    else:
        if status == Aliquot.LARGEST_EXCEEDED:
            print(f"A term exceeds the maximum value ({args.largest_value}).")
        elif status == Aliquot.SIEVE_OVERFLOW:
            print(f"A prime divisor is larger than {args.largest_prime}.")
        else:
            print(f"Undefined status {status}.")
        raise RuntimeError("Abnormal termination.")

    # print("Padded:", seq)       ## debugging
    TITLE = f"aliquot sequence for {args.n}"
    plot(seq, logarithmic = args.log, title=TITLE)

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
