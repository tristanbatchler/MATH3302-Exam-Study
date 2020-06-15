from typing import *
import algorithms
import random
from collections import Counter
from Crypto.Util.number import getPrime

SHOW_WORKING = True

def affine_encrypt(a: int, b: int, *xs: int) -> Tuple[int]:
    """Encrypts a message of one or more plaintext x's in Z_{26} 
       using the given affine key (a, b) both in Z_{26} where a is 
       coprime to 26.

       y ≡ a * x + b    (mod 26)"""
    
    if SHOW_WORKING: print(f"affine_encrypt(a, b, xs) = affine_encrypt({a}, {b}, {xs})")

    ys = []
    for x in xs:
        if SHOW_WORKING: print(f"\ty ≡ a * x + b ≡ {a} * {x} + {b} ≡ {a * x + b} ≡ {(a * x + b) % 26}\t(mod 26)")
        ys.append((a * x + b) % 26)
    
    if SHOW_WORKING: print(f"\tys = {tuple(ys)}")
    return tuple(ys)

def affine_decrypt(a: int, b: int, *ys: int):
    """Decrypts an encrypted message of one or more ciphertext y's in Z_{26}
       using the given affine key (a, b) both in Z_{26} were a is coprime to 
       26.

       x ≡ a^{-1} (y - b)   (mod 26)"""

    if SHOW_WORKING: print(f"affine_decrypt(a, b, ys) = affine_decrypt({a}, {b}, {ys})")
    
    ainv = algorithms.modinverse(a, 26)
    if SHOW_WORKING: print(f"\tFound inverse of a = {a} mod 26. a^(-1) ≡ {ainv}\t(mod 26)")
    xs = []
    for y in ys:
        if SHOW_WORKING: print(f"\tx ≡ a^(-1) * (y - b) ≡ {ainv} * ({y} - {b}) ≡ {ainv * (y - b)} ≡ {(ainv * (y - b)) % 26}\t(mod 26)")
        xs.append((ainv * (y - b)) % 26)
    
    if SHOW_WORKING: print(f"\txs = {tuple(xs)}")
    return tuple(xs)

def rsa_encrypt(n: int, b: int, *xs: int) -> Tuple[int]:
    """Encrypts a message of one or more plaintext x's in Z_{n} 
       using the given RSA public key (n, b) where n is the product 
       of two primes and b is a positive invertible integer less than 
       totient(n).

       y ≡ x^b    (mod n)"""
    if SHOW_WORKING: print(f"rsa_encrypt(n, b, xs) = rsa_encrypt({n}, {b}, {xs})")

    ys = []
    for x in xs:
        if SHOW_WORKING: print(f"\ty ≡ x^b ≡ {x}^{b} (mod {n})")
        ys.append(algorithms.modexp(x, b, n))
    
    if SHOW_WORKING: print(f"\tys = {tuple(ys)}")
    return tuple(ys)

def rsa_decrypt(a: int, p: int, q: int, *ys: int) -> Tuple[int]:
    """Decrypts a message of one or more ciphertext y's in Z_{p*q} 
        using the given RSA private key (a, p, q) where a is the inverse 
        of b from the public key modulo p * q =: n.
        
        x ≡ y ^ a   (mod p * q)"""

    if SHOW_WORKING: print(f"rsa_decrypt(a, p, q, ys) = rsa_decrypt({a}, {p}, {q}, {ys})")

    n: int = p * q
    if SHOW_WORKING: print(f"\tCalculated n := p * q = {p} * {q} = {n}.")
    xs = []
    for y in ys:
        if SHOW_WORKING: print(f"\tx ≡ y^a ≡ {y}^{a} (mod {n})")
        xs.append(algorithms.modexp(y, a, n))
    
    if SHOW_WORKING: print(f"\txs = {tuple(xs)}")
    return tuple(xs)

p = getPrime(8)
q = getPrime(8)
n = p * q
t = (p - 1) * (q - 1)
b = random.randint(1, t - 1)
while algorithms.gcd(b, t) != 1:
    b = random.randint(1, t - 1)
a = algorithms.modinverse(b, t)

xs = tuple([random.randint(0, n - 1) for i in range(8)])

print(f"Public key (n, b) = ({n}, {b})")
print(f"Private key (a, p, q) = ({a}, {p}, {q})")
print(f"Message xs = {xs}")

ys = rsa_encrypt(n, b, *xs)

print(f"Encrypted message ys = {ys}")

xsback = rsa_decrypt(a, p, q, *ys)

print(f"Decrypted back again xs = {xsback}")

print(f"[(n, b) = ({n}, {b}), (a, p, q) = ({a}, {p}, {q})]:\t{xs} -> {ys} -> {xsback}")

assert xs == xsback