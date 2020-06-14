from typing import *
import algorithms
import random
from collections import Counter

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
