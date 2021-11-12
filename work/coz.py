# Based on https://eprint.iacr.org/2010/309.pdf
from typing import Tuple

class GFp:
    def __init__(self, p):
        self.p = p

    def add(self, a, b):
        return (a+b)%self.p
    def sub(self, a, b):
        return (a-b)%self.p
    def mul(self, a, b):
        return (a*b)%self.p
    def sqr(self, a):
        return (a*a)%self.p
    def inv(self, a):
        return pow(a,-1,self.p)

    def __call__(self, v):
        if isinstance(v, IntGFp):
            assert(v.gfp.p == self.p)
            return v
        return IntGFp(self, v)

class IntGFp:
    def __init__(self, gfp: GFp, v: int):
        assert(v < gfp.p)
        self.v = v
        self.gfp = gfp

    def __add__(self, o):
        return IntGFp(self.gfp, self.gfp.add(self.v, o.v))
    def __sub__(self, o):
        return IntGFp(self.gfp, self.gfp.sub(self.v, o.v))
    def __mul__(self, o):
        return IntGFp(self.gfp, self.gfp.mul(self.v, o.v))
    def __eq__(self, o):
        return self.v == o.v
    def __neq__(self, o):
        return self.v != o.v

    def sqr(self):
        return IntGFp(self.gfp, self.gfp.sqr(self.v))
    def inv(self):
        return IntGFp(self.gfp, self.gfp.inv(self.v))

class Curve:
    def __init__(self, gfp: GFp, a: int, b: int, order: int):
        self.gfp = gfp
        self.a = self.gfp(a)
        self.b = self.gfp(b)
        self.order = order

    @property
    def nbits(self): return self.gfp.p.bit_length()

class AffinePt:
    def __init__(self, curve: Curve, x, y):
        self.curve = curve
        self.x = self.gfp(x)
        self.y = self.gfp(y)

    @property
    def gfp(self):
        return self.curve.gfp

    def __eq__(self, o):
        return self.x == o.x and self.y == o.y

class JacobianPt:
    def __init__(self, curve: Curve, x, y, z):
        self.curve = curve
        self.x = self.gfp(x)
        self.y = self.gfp(y)
        self.z = self.gfp(z)

    @property
    def gfp(self):
        return self.curve.gfp

    def to_affine(self) -> AffinePt:
        if (self.z.v == 0):
            return AffinePt(self.curve, 0, 0)
        invZ = self.z.inv()
        invZ2 = invZ.sqr()
        invZ3 = invZ*invZ2

        return AffinePt(self.curve, self.x * invZ2, self.y * invZ3)

    @classmethod
    def from_affine(cls, pt: AffinePt):
        return JacobianPt(pt.curve, pt.x, pt.y, 1)

def double(pt: JacobianPt) -> JacobianPt:
    X1, Y1, Z1 = pt.x, pt.y, pt.z

    N = Z1.sqr()
    E = Y1.sqr()
    B = X1.sqr()
    L = E.sqr()
    S = pt.gfp(2)*((X1+E).sqr() - B - L)
    M = pt.gfp(3)*B + pt.curve.a*N.sqr()

    x = M.sqr() - pt.gfp(2)*S
    y = M*(S-x) - pt.gfp(8)*L
    z = (Y1+Z1).sqr() -E -N 

    return JacobianPt(pt.curve, x,y,z)

def DBLU(pt: JacobianPt) -> Tuple[JacobianPt, JacobianPt]:
    assert(pt.z.v == 1)
    X1, Y1 = pt.x, pt.y

    two = pt.gfp(2)
    three = pt.gfp(3)
    B = X1.sqr()
    E = Y1.sqr()
    L = E.sqr()
    S = two*((X1 + E).sqr() - B - L)
    M = three*B + pt.curve.a

    x = M.sqr() - two*S
    Lm8 = pt.gfp(8)*L
    y = M*(S-x) - Lm8
    z = two*Y1

    ret = JacobianPt(pt.curve,x,y,z)

    # Update to get the same Z as the returned point
    update = JacobianPt(pt.curve, S, Lm8, z)

    return (ret, update)

def ZADDU(self: JacobianPt, o: JacobianPt) -> Tuple[JacobianPt, JacobianPt]:
    assert(self.z == o.z)
    Z = self.z
    X1,X2 = self.x,o.x
    Y1,Y2 = self.y,o.y

    C = (X1 - X2).sqr()
    W1 = X1*C
    W2 = X2*C
    D = (Y1-Y2).sqr()
    A1 = Y1*(W1-W2)
    X3 = D-W1-W2
    Y3 = (Y1-Y2)*(W1-X3)-A1
    Z3 = Z*(X1-X2)

    ret = JacobianPt(self.curve, X3, Y3, Z3)

    # Update current point to get co-Z
    update = JacobianPt(self.curve, W1, A1, Z3)

    return (ret, update)

def add_Z2_1(a: JacobianPt, b: JacobianPt) -> JacobianPt:
    assert(b.z.v == 1)
    gfp = a.curve.gfp
    X1, Y1, Z1 = a.x, a.y, a.z
    X2, Y2     = b.x, b.y

    Z1Z1 = Z1.sqr()
    U2 = X2*Z1Z1
    S2 = Y2*Z1*Z1Z1
    H = U2-X1
    HH = H.sqr()
    I = gfp(4)*HH
    J = H*I
    r = gfp(2)*(S2-Y1)
    V = X1*I
    X3 = r.sqr()-J-gfp(2)*V
    Y3 = r*(V-X3)-gfp(2)*Y1*J
    Z3 = (Z1+H).sqr()-Z1Z1-HH

    return JacobianPt(a.curve, X3, Y3, Z3)


def TPLU(pt: JacobianPt) -> Tuple[JacobianPt, JacobianPt]:
    dbl, pt = DBLU(pt)
    return ZADDU(pt, dbl)

def ZDAU(P: JacobianPt, Q: JacobianPt) -> Tuple[JacobianPt, JacobianPt]:
    assert(P.z == Q.z)
    assert(P.curve == Q.curve)
    X1, Y1, Z = P.x, P.y, P.z
    X2, Y2    = Q.x, Q.y
    two = P.gfp(2)
    four = P.gfp(4)

    Cp = (X1 - X2).sqr()
    W1p = X1*Cp
    W2p = X2*Cp
    Dp = (Y1-Y2).sqr()
    A1p = Y1*(W1p-W2p)
    X3pc = Dp - W1p - W2p
    C = (X3pc - W1p).sqr()
    Y3p = ((Y1-Y2) + (W1p - X3pc)).sqr() - Dp - C - two*A1p
    W1 = four*X3pc*C
    W2 = four*W1p*C
    D = (Y3p-two*A1p).sqr()
    A1 = Y3p*(W1-W2)
    X3 = D-W1-W2
    Y3 = (Y3p-two*A1p)*(W1-X3) - A1
    Z3 = Z*((X1-X2+X3pc-W1p).sqr() - Cp - C)
    Dc = (Y3p+two*A1p).sqr()
    X2 = Dc - W1 - W2
    Y2 = (Y3p + two*A1p)*(W1-X2)-A1

    R = JacobianPt(P.curve, X3, Y3, Z3)
    Q = JacobianPt(P.curve, X2, Y2, Z3)

    return (R,Q)

def scalar_mult(x: int, pt: JacobianPt) -> JacobianPt:
    assert((pt.curve.order & 1) == 1)
    assert(pt.z.v == 1)
    # Idea from https://www.iacr.org/archive/ches2007/47270135/47270135.pdf
    # We force x to be odd, and then subtract P at the end. If Z is 1, this
    # requires less operations than an additional ZDAU (required if we add the
    # order of the curve to x).
    subP = (x&1) == 0
    R = [0,0]

    b = (x >> 1) & 1
    R[1-b], R[b] = TPLU(pt)

    for i in range(2, pt.curve.nbits):
        b = (x >> i) & 1
        R[1-b], R[b] = ZDAU(R[1-b], R[b])

    ret = R[0]
    Rsub = add_Z2_1(ret, JacobianPt(pt.curve, pt.x, pt.gfp(-pt.y.v), pt.z))
    return Rsub if subP else ret

from Crypto.PublicKey import ECC

P256 = Curve(GFp(0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff), -3, 0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b, 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551)
G = AffinePt(P256, 0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296, 0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5)
JG = JacobianPt.from_affine(G)

JP, JGz = DBLU(JG)
assert(JGz.to_affine() == G)
P = JP.to_affine()
PREF = ECC.construct(curve="P-256", d=2).pointQ
assert(P.x.v == PREF.x and P.y.v == PREF.y)

JP, JGz = DBLU(JG)
JP, JGz = ZADDU(JGz, JP)
assert(JGz.to_affine() == G)
P = JP.to_affine()
PREF = ECC.construct(curve="P-256", d=3).pointQ
assert(P.x.v == PREF.x and P.y.v == PREF.y)

JP, JGz = TPLU(JG)
P = JP.to_affine()
assert(JGz.to_affine() == G)
assert(P.x.v == PREF.x and P.y.v == PREF.y)

import random
x = random.getrandbits(255)
print(hex(x))
JP = scalar_mult(x, JG)
print(hex(JP.x.v), hex(JP.y.v), hex(JP.z.v))
P = JP.to_affine()
print(hex(P.x.v), hex(P.y.v))
PREF = ECC.construct(curve="P-256", d=x).pointQ
print(hex(PREF.x), hex(PREF.y))
