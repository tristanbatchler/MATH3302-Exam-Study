from typing import *
import math

SHOW_WORKING: bool = True

def gcd(a: int, b: int) -> int:
    if (a, b) == (0, 0):
        raise ValueError("Can't find gcd(0, 0)")
    if a < 0 or b < 0:
        raise ValueError("a and b must be non-negative")
    
    if SHOW_WORKING: print(f"Find gcd(a, b) = gcd({a}, {b})")
    if b > a:
        if SHOW_WORKING: print(f"\tb > a. Set r1 := b = {b} and r2 := a = {a} so that r1 > r2")
        r1, r2 = b, a
    else:
        if SHOW_WORKING: print(f"\ta >= b. Set r1 := a = {a} and r2 := b = {b} so that r1 >= r2")
        r1, r2 = a, b        

    if SHOW_WORKING: print(f"\tProceeding with algorithm until r2 hits 0. gcd({a}, {b}) will be the ending r1 value:")
    while r2 != 0:
        if SHOW_WORKING: print(f"\t\tSet q := floor(r1 / r2) = floor({r1} / {r2}) = floor({round(r1 / r2, 2)}) = {r1 // r2}")
        q: int = r1 // r2

        if SHOW_WORKING: print(f"\t\tSet (r1, r2) := (r2, r1 - r2 * q) = ({r2}, {r1} - {r2} * {q}) = ({r2}, {r1 - r2 * q})")
        r1, r2 = r2, r1 - r2 * q

        if SHOW_WORKING: print(f"\t\t -")
    
    if SHOW_WORKING: print(f"\tStopping condition hit (r2 = 0). Result of gcd({a}, {b}) is r1 = {r1}")

    return r1

# TODO: Finish this
def modinverse(a: int, m: int) -> int:
    """Runs the extended Eulidean algorithm to find the inverse of a modulo m"""
    if a <= 0 or m <= 0:
        raise ValueError("a and m must be positive")
    
    if SHOW_WORKING: print(f"Find gcd(a, b) = gcd({a}, {m})")
    if m > a:
        if SHOW_WORKING: print(f"\tb > a. Set r1[0] := m = {m} and r2[0] := a = {a} so that r1[0] > r2[0")
        r1s, r2s = [m], [a]
    else:
        if SHOW_WORKING: print(f"\ta >= b. Set r1[0] := a = {a} and r2[0] := m = {m} so that r1[0] >= r2[0]")
        r1s, r2s = [a], [m]       

    if SHOW_WORKING: print(f"\tProceeding with algorithm until r2 hits 0. gcd({a}, {m}) will be the ending r1 value:")
    qs = []
    i = 0
    while r2s[-1] != 0:
        i += 1

        if SHOW_WORKING: print(f"\t\tSet q[{i - 1}] := floor(r1[{i - 1}] / r2[{i - 1}]) = floor({r1s[i - 1]} / {r2s[i - 1]}) = floor({round(r1s[i - 1] / r2s[i - 1], 2)}) = {r1s[i - 1] // r2s[i - 1]}")
        qs.append(r1s[i - 1] // r2s[i - 1])

        if SHOW_WORKING: print(f"\t\tSet (r1[{i}], r2[{i}]) := (r2[{i - 1}], r1[{i - 1}] - r2[{i - 1}] * q[{i - 1}]) = ({r2s[i - 1]}, {r1s[i - 1]} - {r2s[i - 1]} * {qs[i - 1]}) = ({r2s[i - 1]}, {r1s[i - 1] - r2s[i - 1] * qs[i - 1]})")
        r1, r2 = r2s[i - 1], r1s[i - 1] - r2s[i - 1] * qs[i - 1]
        r1s.append(r1)
        r2s.append(r2)

        if SHOW_WORKING: print("\t\t -")
    
    if SHOW_WORKING: print(f"\tStopping condition hit (r2[{i}] = 0). Result of gcd({a}, {m}) is r1[{i}] = {r1s[-1]}")

    if r1s[-1] != 1:
        if SHOW_WORKING: print(f"\t{a} has no inverse modulo {m} because gcd({a}, {m}) = {r1s[-1]} != 1 (they must be coprime)")
        return None

    if SHOW_WORKING: print(f"\n\tBegin working backwards:")

    i -= 1
    r2siminus1strgeneral = f"(r2[{i - 1}])"
    r2siminus1strexact = f"({r2s[i - 1]})"
    while i >= 0:
        if SHOW_WORKING: print(f"r2[{i}] = r1[{i - 1}] - r2[{i - 1}] * q[{i - 1}] = {r1s[i - 1]} - {r2s[i - 1]} * {qs[i - 1]} = {r2s[i]}")
        if SHOW_WORKING: print(f"r2[{i}] = r1[{i - 1}] - {r2siminus1strgeneral} * q[{i - 1}] = {r1s[i - 1]} - {r2siminus1strexact} * {qs[i - 1]} = {r2s[i]}")
        if SHOW_WORKING: print(" -")

        r2siminus1strgeneral = f"(r1[{i - 1}] - {r2siminus1strgeneral} * q[{i - 1}])"
        r2siminus1strexact = f"({r1s[i - 1]} - {r2siminus1strexact} * {qs[i - 1]})"

        i -= 1

def modexp(x: int, b: int, p: int) -> int:
    """Calculates x^b modulo p with the fast exponentiation algorithm"""
    if p <= 0:
        raise ValueError("p must be greater than zero")

    if SHOW_WORKING: print(f"Find modexp(x, b, p) = gcd({x}, {b}, {p})")

    if SHOW_WORKING: print("\tSet y := 1")
    y: int = 1

    if SHOW_WORKING: print("\tRepeat until b hits zero:")
    while b > 0:
        if b % 2 == 1:
            if SHOW_WORKING: print(f"\n\t\tb = {b} > 0 is odd so do the following: ")
            if SHOW_WORKING: print(f"\t\t\t * Set y to (y * x) mod p = ({y} * {x}) mod {p} = {y * x} mod {p} = {(y * x) % p}")
            y = (y * x) % p
            if SHOW_WORKING: print(f"\t\t\t * Set b to b - 1 = {b} - 1 = {b - 1}")
            b -= 1
        else:
            if SHOW_WORKING: print(f"\n\t\tb = {b} > 0 is even so do the following: ")
            if SHOW_WORKING: print(f"\t\t\t * Set x to x^2 mod p = {x}^2 mod {p} = {x ** 2} mod {p} = {(x ** 2) % p}")
            x = (x ** 2) % p
            if SHOW_WORKING: print(f"\t\t\t * Set b to b / 2 = {b} / 2 = {b / 2}")
            b = b / 2
        
    print(f"\tb hit zero, stop looping and return y = {y}")

    return y

def prime_factors(n: int) -> Dict[int, int]:
    """Returns a dictionary containing the prime factors of n as keys and their multiplicity as values
       e.g. prime_factors(18) -> {2: 1, 3: 2}"""
    factors = {}

    while n % 2 == 0:
        print(f"\tChecking if {n} divides 2")
        print(f"\t\tYes--Adding 2")
        if 2 in factors.keys():
            factors[2] += 1
        else:
            factors[2] = 1
        n //= 2

    checklimit: int = math.ceil(math.sqrt(n)) + 1
    for d in range(3, checklimit, 2):
        if n % d:
            print(f"\tChecking if {n} divides {d}")
            print(f"\t\tNo--moving on")
            d += 1
        else:
            while n % d == 0:
                print(f"\tChecking if {n} divides {d}")
                print(f"\t\tYes--Adding {d}")
                if d in factors.keys():
                    factors[d] += 1
                else:
                    factors[d] = 1
                n //= d
    if n > 1:
        factors[n] = 1
    return factors

def totient(n: int) -> int:
    """Returns the number of primes less than n"""
    pass

print(prime_factors(1025))