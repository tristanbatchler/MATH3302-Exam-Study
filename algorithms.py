from typing import *
import math
import sympy.simplify
import sympy.core.numbers

SHOW_WORKING: bool = True

class WrongNumberOfArgumentsError(Exception):
    pass

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

def setgcd(*ns: int) -> int:
    """Calculates the greatest common divisor of all elements in the given list"""
    if SHOW_WORKING: print(f"setgcd{ns}")

    if len(ns) < 2:
        raise WrongNumberOfArgumentsError("setgcd must supply two or more arguments")

    if len(ns) == 2:
        result = gcd(*ns)
    else:
        result: int = setgcd(gcd(*ns[:2]), *ns[2:])
    if SHOW_WORKING: print(f"\tsetgcd{ns} = {result}")
    return result

def coprime(*ns: int) -> bool:
    if SHOW_WORKING: print(f"coprime{ns}")
    resultgcd = setgcd(*ns)
    if SHOW_WORKING: print(f"coprime{ns} = {resultgcd == 1} because setgcd{ns} == {resultgcd} {'=' if resultgcd == 1 else '!'}= 1")
    return resultgcd == 1

def pairwise_coprime(*ns: int) -> bool:
    return coprime(*ns)

# TODO: Finish this
def modinverse(a: int, m: int) -> int:
    """Runs the extended Eulidean algorithm to find the inverse of a modulo m"""
    if SHOW_WORKING: print(f"modinverse(a, m) = modinverse({a}, {m})")
    if SHOW_WORKING: print(f"\tWe want to find some x & y such that {a} * x + {m} * y = 1")

    if a < 0 or m <= 0:
        raise ValueError("a must be non-negative and m must be positive")

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

    def getnestedexpressionstr(leftstr: str, nestedr1r2q: List[Union[int, List[int]]], rightstr: str) -> str:
        if SHOW_WORKING: print(f"\t\tgetnestedexpressionstr('{leftstr}', {nestedr1r2q}, '{rightstr}')")
        r1: int = nestedr1r2q[0]
        r2: Union[int, List[int]] = nestedr1r2q[1]
        q: int = nestedr1r2q[2]
        if SHOW_WORKING: print(f"\t\t\tr1 = {r1}, r2 = {r2}, q = {q}")

        if isinstance(r2, int):
            return f"{leftstr}{r1} - {r2} * {q}{rightstr}"
        
        if leftstr == rightstr == '':
            return getnestedexpressionstr(f"{r1} - (", r2, f") * {q}")

        return getnestedexpressionstr(f"{leftstr}{r1} - (", r2, f") * {q}{rightstr}")

    def backtrack(index: int, nestedr1r2q: List[Union[int, List[int]]]) -> List[Union[int, List[int]]]:
        """Provided an index and an ordered list representing the r1, r2, and q values of the equation
           r1 - r2 * q, this function returns another list where r2 has been broken down to the parts of 
           its equation on the previous indexed equation, e.g. if the 3rd and 4th equations from the GCD 
           algorithm are:
               (3): r1 - r2 * q2 = 4 - 4 * 1
               (4): r1 - r2 * q2 = 3 - 1 * 3
           then: 
               backtrack(4, [3, 1, 3]) -> [3, [4, 3, 1], 3].
           
           This also works when the middle element of the list (the r2 element) is given as a list of parts,
           e.g., if we follow the previous example where additionally equation 2 is:
               (2): r1 - r2 * q2 = 11 - 4 * 2
           then:
               backtrack(3, [3, [4, 3, 1], 3]) -> [3, [4, [11, 4, 2], 1], 3]."""
           
        if SHOW_WORKING: print(f"\t\tbacktrack({index}, {nestedr1r2q})")

        if index <= 0:
            raise ValueError("Can't backtrack from here, please supply a positive index")
        
        r1: int = nestedr1r2q[0]
        r2: Union[int, List[int]] = nestedr1r2q[1]
        q: int = nestedr1r2q[2]

        if index == 1:
            return [r1, [r1s[0], r2s[0], qs[0]], q]

        return [r1, backtrack(index - 1, [r1s[index - 1], r2s[index - 1], qs[index - 1]]), q]

    if i - 2 > 0:
        expression = backtrack(i - 2, [r1s[i - 2], r2s[i - 2], qs[i - 2]])

        nestedexpressionstr: str = getnestedexpressionstr('', expression, '')
        nestedexpressionstr = nestedexpressionstr.replace(str(a), 'a').replace(str(m), 'm')

        if SHOW_WORKING: print(f"\t\t{nestedexpressionstr}")
        if SHOW_WORKING: print(f"\t\t{sympy.simplify(nestedexpressionstr)}")

    x, y = sympy.core.numbers.igcdex(a, m)[:2]
    if SHOW_WORKING: print(f"\ta * x + m * y = 1 -> {a} * {x} + {m} * {y} = 1")

    if SHOW_WORKING: print(f"\tmodinverse({a}, {m}) = {x}\t(mod {m}) = {x % m}")
    
    return x % m
    
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
    if SHOW_WORKING: print(f"prime_factors({n})")
    original_n = n
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

    print(f"\t{original_n} has prime factorisation {' * '.join([str(p) + '^' + str(e) for p, e in factors.items()])}")
    return factors

def totient(n: int) -> int:
    """Returns the number of primes less than n"""
    if SHOW_WORKING: print(f"totient({n})")
    print(f"\ttotient(n) = p1^(e1 - 1) * (p1 - 1) * p2^(e2 - 1) * (p2 - 1) * ... * pk^(ek - 1) * (pk - 1)")

    result: int = 1
    
    working_str1: str = ''
    working_str2: str = ''
    working_str3: str = ''
    working_str4: str = ''

    n_prime_factors: Dict[int, int] = prime_factors(n)

    for p, e in n_prime_factors.items():
        result *= p ** (e - 1) * (p - 1)
        working_str1 += f"{p}^({e} - 1) * ({p} - 1) * "
        working_str2 += f"{p}^{e - 1} * {p - 1} * "
        working_str3 += f"{p ** (e - 1)} * {p - 1} * "
        working_str4 += f"{p ** (e - 1) * (p - 1)} * "
    
    working_str1 = working_str1[:-2]
    working_str2 = working_str2[:-2]
    working_str3 = working_str3[:-2]
    working_str4 = working_str4[:-2]

    n_prime_factorisation: str = ' * '.join([str(p) + '^' + str(e) for p, e in n_prime_factors.items()])

    if SHOW_WORKING: print(f"\n\ttotient({n}) = totient({n_prime_factorisation}) \n\t\t= {working_str1} \n\t\t= {working_str2} \n\t\t= {working_str3} \n\t\t= {working_str4} \n\t\t= {result}")

    return result

def chinese_remainder(ms: List[int], _as: List[int]) -> int:
    """Computes the unique solution to the r simultaneous congruences modulo the product of pairwise coprime ms
       x = a[0] % M
       x = a[1] % M
        ...
       x = a[r] % M"""

    if SHOW_WORKING: print(f"chinese_remainder(ms, _as) = chinese_remainder({ms}, {_as})")

    for m in ms:
        if m <= 0:
            raise ValueError("All ms must be positive integers")
    
    if not pairwise_coprime(*ms):
        raise ValueError("ms must be pairwise coprime")

    if len(ms) != len(_as):
        raise ValueError("Must supply same number of _as as ms")

    r: int = len(ms)

    if SHOW_WORKING:
        print(f"\tFinding a solution to the following {r} congruences:")
        for i in range(r):
            print(f"\t\tx ≡ {_as[i]}\t(mod {ms[i]})")
    
    M = 1
    for i in range(r):
        M *= ms[i]
    
    if SHOW_WORKING: print(f"\tLet M = m1 * ... * m{r} = {ms[0]} * ... * {ms[r - 1]} = {M}")

    Mis = []
    yis = []
    if SHOW_WORKING: print(f"\tWe must find each Mi = M / mi and yi = modinverse(Mi, mi)")
    for i in range(r):
        Mis.append(int(M / ms[i]))
        if SHOW_WORKING: print(f"\t\tM{i} = M / m{i} = {M} / {ms[i]} = {Mis[i]}")

        if SHOW_WORKING: print(f"\t\ty{i} = modinverse(M{i}, m{i}) = modinverse({Mis[i]}, {ms[i]}) = modinverse({Mis[i] % ms[i]}, {ms[i]})")
        yis.append(modinverse(Mis[i] % ms[i], ms[i]))

    if SHOW_WORKING: print(f"\tSo Mis = {Mis} and yis = {yis}")
    if SHOW_WORKING: print(f"\tCRT -> x ≡ a1 * M1 * y1 + ... + ar * Mr * yr)\t(mod M)")
    _sum: int = 0
    for i in range(r):
        _sum += _as[i] * Mis[i] * yis[i]

    print(f"\tThis is x ≡ {_sum} % {M} = {_sum % M}")
    return _sum % M
            

        
    
chinese_remainder([3, 5, 7], [1, 2, 3])