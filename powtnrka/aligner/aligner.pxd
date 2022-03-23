# Copyright (C) 2024   Alvin Hsu
"""
Implementing IUPAC mixed base nomenclature with bitwise operations. Each base
is represented by a different bit:
A = 0b0000 0001
C = 0b0000 0010
G = 0b0000 0100
T = 0b0000 1000
Y = C | T, R = A | G, etc.
"""
ctypedef enum Base:
    A = 0x01
    C = 0x02
    G = 0x04
    T = 0x08
    M = A | C
    R = A | G
    W = A | T
    S = C | G
    Y = C | T
    K = G | T
    B = C | G | T
    D = A | G | T
    H = A | C | T
    V = A | C | G
    N = A | C | G | T


"""
Implementing directions with bitwise operations. Each direction is represented by a different bit:
A = 0b0000 0001
C = 0b0000 0010
G = 0b0000 0100
T = 0b0000 1000
Y = C | T, R = A | G, etc.
"""
ctypedef size_t Direction
