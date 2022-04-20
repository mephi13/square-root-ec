#!/bin/python
#for prime generation
from sympy import randprime, isprime, nextprime
#misc
import argparse
from blessed import Terminal
from random import getrandbits, randrange
import time


def is_qr(a, p):
    """Check if a number is a quadratic residue mod p"""
    exp = (p-1)//2
    assert exp * 2 + 1 == p, "Exponent calculations are wrong in QR checker"
    val = pow(a, exp, p) # a^[(p-1)/2]
\
    assert val == 1 or val == -1 % p, "QR checker didn't return 1 or -1!"

    return True if val == 1 else False

def gen_p(min_s, max_s, min_t, max_t, p):
    """   
    Generate a prime p of form 2^t * s + 1 where s is odd
    """

    # unitialized p has non prime value to be sure we enter the loop
    p = 4 if not p else p
    s, t = (None, None)
    while not isprime(p):
        # set s as a random odd number greater than 2^99
        s = getrandbits(randrange(min_s, max_s))
        # make sure s always add
        s = s + 1 if s % 2 == 0 else s
        assert s % 2 == 1, "s value is not odd!"

        # set t as a random exponent
        t = randrange(min_t, max_t)

        p = pow(2, t) * s + 1

    #print(f"Bit length: {p.bit_length()}")

    assert isprime(p), "p is not prime!"

    if s and t:
        print(f"s = {'{0:b}'.format(s)}, t = {t}")

    return p

def gen_a(a, p, real_x = None):
    """
    generate a such that a = x^2 mod p
    x will be then calculated using square root algorithm 
    """
    if not a:
        x = randrange(1, p) if not real_x else real_x
        a = pow(x, 2, p)
    else:
        x = None
    return (a, x)

def calculate_s_and_t(p):
    """Calculate s and t for given p value"""
    p_even = p - 1
    t = 0

    # if p is even that means were still calculating t
    while p_even % 2 == 0:
        t += 1
        p_even = p_even // 2

    # everything thats left is s
    s = p_even

    assert pow(2, t) * s + 1 == p, "s and t calculations are wrong, p != 2^t * s + 1"
    return (s, t)

def repeated_squaring(base, exp, p):
    """Repeated squaring algorithm for faster square root."""
    exp_bin = "{0:b}".format(exp)
    exp_len = len(exp_bin)

    x = base
    if x == -1 % p:
        return x % p, 0
    for j in range(1, exp_len):
        x = pow(x, 2, p)
        # if bit is odd, multiply
        if exp_bin[j] == "1":
            x = (x * base) % p

        # do it until the value is equal to -1
        if x == -1 % p:
            return x % p, j
    return x % p, exp_len - 1

def square_root(a, p):
    s, t = calculate_s_and_t(p)

    # set r = a ^[(s+1)/2] mod p
    r = pow(a, (s+1)//2, p)

    # n is any NQR mod p
    n = randrange(1, p)

    while is_qr(n, p):
        n = randrange(1, p)

    # set b = n^s mod p
    b = pow(n, s, p)
    # b is a generator of <g2>
    # ord(b) = 2^t
    # ord(j) = 2^t
    
    # we need to calculate j exponent such that:
    # a^-1 * (b^j *r )^2 = 1 mod p, then 
    # b^j * r is square root of a mod p

    # j = sum_0^(t-2) 2^k * j_k for j_k = {0, 1}
    
    # calculate j_0
    # calculate (a^-1 * r^2)^(2^(t-2)) mod p

    eq_3 = pow(pow(a, -1, p) * pow(r, 2, p), pow(2, t-2), p)
    # if eq_3 == -1 the bit value is 1, if eq_3 == 1 the bit value is 0
    assert eq_3 == 1 or eq_3 == -1 % p, "Square root calculations are wrong!"

    bit_value = 1 if eq_3 == -1 % p else 0
    j = pow(2, 0) * bit_value

    k = 1
    while k <= (t-2):
        # calculate j_k
        # (a^-1 * (b^(j) * r)^2)^(2^(t-k-2))
        
        # get i that is 2^i mod p== -1
        eq_3, i = repeated_squaring( pow(a, -1, p) * pow( pow(b, j, p) * r, 2 ), pow(2, t-k-2), p)

        assert eq_3 == 1 or eq_3 == -1 % p, "Square root calculations are wrong!" 

        if eq_3 == -1 % p:
            k_prim = t - i - 2

            # all of the bits between k and k_prim are 0
            k = k_prim

            # last bit is 1
            j += pow(2, k) 
        # eq_1 == 1
        else:
            pass
            # bit is 0, we don't add it to the sum
        k += 1

    # calculate x value
    x = (pow(b, j, p) * r) % p

    return x, s, t

def generate_square_root_instance(min_s_bits, max_s_bits, min_t_bits, max_t_bits, real_x, p, a):
    """
    Generate an instance of a square root problem
    """
    try:
        # generate p
        p = gen_p(min_s_bits, max_s_bits, min_t_bits, max_t_bits, p)

        # generate a
        a, real_x = gen_a(a, p, real_x)

        return (p, a, real_x)

    except Exception as ex:
        raise Exception(f"Generating square root problem failed, reason: \n{ex.args[0]}")


def main(s_bits = None, t_bits = None, real_x = None, p =None, a=None):
    try:
        #init params
        min_s_bits = s_bits if s_bits else 99
        max_s_bits = s_bits + 1 if s_bits else 200

        min_t_bits = t_bits if t_bits else 150
        max_t_bits = t_bits + 1 if t_bits else 200

        # generate square root instance
        p, a, real_x = generate_square_root_instance(min_s_bits, max_s_bits, min_t_bits, max_t_bits, real_x, p, a)

        print(p, isprime(p))

        term = Terminal()
        print(f'p is a {p.bit_length()} bit number')
        print(f'printing p = 2^{term.pink("t")} * {term.purple("s")} + 1 ...')
        print(f'{term.red("{0:b}".format(p))}')

        if real_x:
            print(f"Real x = {term.blue(str(real_x))},") 
            print("{0:b}".format(real_x))

        print(f"a = {term.blue(str(real_x))}^2"\
            f" mod {term.red(str(p))}\na = {term.green(str(a))}")\
            if real_x else print(f"a = {term.blue('x')}^2 mod {term.red('p')}\na = {term.green(str(a))}")

        print(f"Starting square-root(a, p)")
        
        start_time = time.time()
        calc_x, s, t = square_root(a, p)
        calc_time = time.time() - start_time

        print(f'calculated s = {term.pink("{0:b}".format(s))}')
        print(f'calculated t = {term.purple(str(t))}')
        print(f'Visually check if s and t are correct:')
        print(f'{term.pink("{0:b}".format(s))}{term.purple("0" * (t-1)+"1")}')
        
        print(f'Calculation time: {round(calc_time, 1)} seconds')
        if real_x:
            print(f"Real  x = \n{term.blue(str(real_x))}, \nReal -x = \n{term.blue(str(-real_x % p))},")
        print(f"Calculated  x = \n{term.turquoise(str(calc_x))}, \nCalculated -x = \n{term.turquoise(str(-calc_x % p))}") 

        print("")
        #Check if real x and calc x are the same
        if pow(calc_x, 2, p) == a:
            print(f"{term.turquoise(str(calc_x))}^2 mod {term.red(str(p))} = {term.green(str(a))}")
            print(f"{term.green('Calculations correct!')}")
        else:
            print(f"{term.turquoise(str(calc_x))}^2 mod {term.red(str(p))} != {term.green(str(a))}")
            print(f"{term.red('Calculations failed!')}")
    except Exception as ex:
        print(ex.args[0])

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Square root algorhitm.')
    parser.add_argument("-a", type=int)
    parser.add_argument("-s", type=int)
    parser.add_argument("-t", type=int)
    parser.add_argument("-p", type=int)
    parser.add_argument("-x", type=int)
    a = parser.parse_args()


    #print(fast_multiply(7, 8))
    #print(fast_multiply(EC[40].basepoint, 16, EC[40].a ))

    main(a.s, a.t, a.x, a.p, a.a)