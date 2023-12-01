#!/usr/bin/env python
# coding: utf-8

# ## Functions fot the Limits on Bosonic Currents operator $-$ $\psi^2 H^2 D$
#
# Specific functions fot the HNL decay widths needed fot the re-scaling of the CC+NC experimental bounds.

from math import sqrt, pi, log
from scipy import integrate


import mpmath as mp
from mpmath import *

mp.dps = 100
mp.pretty = True


def l(a, b, c):
    z = a**2 + b**2 + c**2 - 2 * a * b - 2 * a * c - 2 * b * c
    return z


def xfrac(m, M):
    z = mpf(m) / mpf(M)
    return z


def L(x):
    z = log((1 - 3 * x**2 - (1 - x**2) * sqrt(1 - 4 * x**2)) / (x**2 * (1 + sqrt(1 - 4 * x**2))))
    return z


def Gammavvv(M, Ue, Umu, Utau):
    GF = 1.16e-5
    A = (GF**2) / (192 * (pi**3))
    z = (Ue**2 + Umu**2 + Utau**2) * A * (M**5)
    return z


def f1(x):
    z = (1 - 14 * (x**2) - 2 * (x**4) - 12 * (x**6)) * sqrt(1 - 4 * (x**2)) + 12 * (x**4) * ((x**4) - 1) * L(x)
    return z


def f2(x):
    z = 4 * ((x**2) * (2 + 10 * (x**2) - 12 * (x**4)) * sqrt(1 - 4 * (x**2)) + 6 * (x**4) * (1 - 2 * (x**2) + 2 * (x**4)) * L(x))
    return z


def Gammavll1(ml, M):
    GF = 1.16e-5
    A = (GF**2) / (192 * (pi**3))
    sw2 = 2.3e-1
    C1lep = 0.25 * (1 - 4 * sw2 + 8 * (sw2**2))
    C2lep = 0.5 * (2 * (sw2**2) - sw2)
    z = A * (M**5) * ((C1lep * f1(xfrac(ml, M)) + C2lep * f2(xfrac(ml, M))))
    return z


def Gammavll2(ml, M):
    GF = 1.16e-5
    A = (GF**2) / (192 * (pi**3))
    sw2 = 2.3e-1
    C1lep = 0.25 * (1 - 4 * sw2 + 8 * (sw2**2))
    C2lep = 0.5 * (2 * (sw2**2) - sw2)
    z = A * (M**5) * sw2 * (f2(xfrac(ml, M)) + 2 * f1(xfrac(ml, M)))
    return z


def Gammavee(M, Ue, Umu, Utau):
    me = 5.11e-4
    if M >= 2 * me:
        z = (Ue**2 + Umu**2 + Utau**2) * Gammavll1(me, M) + (Ue**2) * Gammavll2(me, M)
    else:
        z = 0
    return z


def Gammavmumu(M, Ue, Umu, Utau):
    mmu = 1.05e-1
    if M >= 2 * mmu:
        z = (Ue**2 + Umu**2 + Utau**2) * Gammavll1(mmu, M) + (Umu**2) * Gammavll2(mmu, M)
    else:
        z = 0
    return z


def Gammavllp(M, m):
    GF = 1.16e-5
    A = (GF**2) / (192 * (pi**3))
    z = A * (M**5) * (1 - 8 * (xfrac(m, M) ** 2) + 8 * (xfrac(m, M) ** 6) - xfrac(m, M) ** 8 - 12 * (xfrac(m, M) ** 4) * log(xfrac(m, M) ** 2))
    return z


def Gammavemmup(M, Ue, Umu, Utau):
    me = 5.11e-4
    mmu = 1.05e-1
    if M >= (me + mmu):
        z = (Ue**2) * Gammavllp(M, mmu)
    else:
        z = 0
    return z


def Gammavmumep(M, Ue, Umu, Utau):
    me = 5.11e-4
    mmu = 1.05e-1
    if M >= (me + mmu):
        z = (Umu**2) * Gammavllp(M, mmu)
    else:
        z = 0
    return z


def Gammavemtaup(M, Ue, Umu, Utau):
    me = 5.11e-4
    mtau = 1.777
    if M >= (me + mtau):
        z = (Ue**2) * Gammavllp(M, mtau)
    else:
        z = 0
    return z


def Gammavtaumep(M, Ue, Umu, Utau):
    me = 5.11e-4
    mtau = 1.777
    if M >= (me + mtau):
        z = (Utau**2) * Gammavllp(M, mtau)
    else:
        z = 0
    return z


def Gammalep(M, Ue, Umu, Utau):
    z = (
        Gammavvv(M, Ue, Umu, Utau)
        + Gammavee(M, Ue, Umu, Utau)
        + Gammavmumu(M, Ue, Umu, Utau)
        + Gammavemmup(M, Ue, Umu, Utau)
        + Gammavmumep(M, Ue, Umu, Utau)
        + Gammavemtaup(M, Ue, Umu, Utau)
        + Gammavtaumep(M, Ue, Umu, Utau)
    )
    return z


def inte(s, x, y, z):
    z = 12 * (1 / s * (s - x - y) * (1 + z - s) * sqrt(l(s, x, y) * l(1, s, z)))
    return z


def IntCC(x, y, z):
    z = integrate.quad(inte, (sqrt(x) + sqrt(y)) ** 2, (1 - sqrt(z)) ** 2, args=(x, y, z))
    return z[0]


def GammaveeCC(M, Ue, Umu, Utau):
    me = 5.11e-4
    GF = 1.16e-5
    A = (GF**2) / (192 * (pi**3))
    if M >= 2 * me:
        z = A * (M**5) * (Ue**2) * IntCC(0, (xfrac(me, M)) ** 2, (xfrac(me, M)) ** 2)
    else:
        z = 0
    return z


def GammavmumuCC(M, Ue, Umu, Utau):
    mmu = 1.05e-1
    GF = 1.16e-5
    A = (GF**2) / (192 * (pi**3))
    if M >= 2 * mmu:
        z = A * (M**5) * (Umu**2) * IntCC(0, (xfrac(mmu, M)) ** 2, (xfrac(mmu, M)) ** 2)
    else:
        z = 0
    return z


def GammaVl(f, mv, V, M, ml):
    GF = 1.16e-5
    z = (
        (GF**2)
        * (M**3)
        / (16 * pi * (mv**2))
        * (f**2)
        * (V**2)
        * sqrt((1 + ((mv**2) - (ml**2)) / (M**2)) ** 2 - 4 * (mv**2) / (M**2))
        * ((1 - (mv**2) / (M**2)) * (1 + 2 * (mv**2) / (M**2)) + (ml**2) / (M**2) * (((mv**2) + (ml**2)) / (M**2) - 2))
    )
    return z


def GammaVnu(f, mv, g, M, Ue, Umu, Utau):
    GF = 1.16e-5
    z = (
        (Ue**2 + Umu**2 + Utau**2)
        * (GF**2)
        * (M**3)
        / (32 * pi * (mv**2))
        * (f**2)
        * (g**2)
        * (1 + 2 * (mv**2) / (M**2))
        * (1 - (mv**2) / (M**2)) ** 2
    )
    return z


def GammaPl(f, mp, V, M, ml):
    GF = 1.16e-5
    z = (
        (GF**2)
        * (M**3)
        / (16 * pi)
        * (f**2)
        * (V**2)
        * sqrt((1 + ((mp**2) - (ml**2)) / (M**2)) ** 2 - 4 * (mp**2) / (M**2))
        * (1 - (mp**2) / (M**2) - (ml**2) / (M**2) * (2 + ((mp**2) - (ml**2)) / (M**2)))
    )
    return z


def GammaPnu(f, mp, M, Ue, Umu, Utau):
    GF = 1.16e-5
    z = (Ue**2 + Umu**2 + Utau**2) * (GF**2) * (M**3) / (32 * pi) * (f**2) * (1 - (mp**2) / (M**2)) ** 2
    return z


def GammaPhi(M, Ue, Umu, Utau):
    mPhi = 1.02
    fPhi = 0.232
    sw2 = 2.3e-1
    gPhi = sqrt(2) * (0.5 - 2 * sw2 / 3)
    if M >= mPhi:
        z = GammaVnu(fPhi, mPhi, gPhi, M, Ue, Umu, Utau)
    else:
        z = 0
    return z


def GammaRho(M, Ue, Umu, Utau):
    mRho = 0.775
    fRho = 0.17
    sw2 = 2.3e-1
    gRho = 1 - 2 * sw2
    if M >= mRho:
        z = GammaVnu(fRho, mRho, gRho, M, Ue, Umu, Utau)
    else:
        z = 0
    return z


def GammaOmega(M, Ue, Umu, Utau):
    mOmega = 0.782
    fOmega = 0.155
    sw2 = 2.3e-1
    gOmega = 2 * sw2 / 3
    if M >= mOmega:
        z = GammaVnu(fOmega, mOmega, gOmega, M, Ue, Umu, Utau)
    else:
        z = 0
    return z


def GammaRhoe(M, Ue, Umu, Utau):
    mRho = 0.775
    fRho = 0.17
    Vud = 9.737e-1
    me = 5.11e-4
    if M >= (mRho + me):
        z = (Ue**2) * GammaVl(fRho, mRho, Vud, M, me)
    else:
        z = 0
    return z


def GammaRhomu(M, Ue, Umu, Utau):
    mRho = 0.775
    fRho = 0.17
    Vud = 9.737e-1
    mmu = 1.05e-1
    if M >= (mRho + mmu):
        z = (Umu**2) * GammaVl(fRho, mRho, Vud, M, mmu)
    else:
        z = 0
    return z


def GammaKse(M, Ue, Umu, Utau):
    mKs = 0.892
    fKs = 0.177
    Vus = 2.243e-1
    me = 5.11e-4
    if M >= (mKs + me):
        z = (Ue**2) * GammaVl(fKs, mKs, Vus, M, me)
    else:
        z = 0
    return z


def GammaKsmu(M, Ue, Umu, Utau):
    mKs = 0.892
    fKs = 0.177
    Vus = 2.243e-1
    mmu = 1.05e-1
    if M >= (mKs + mmu):
        z = (Umu**2) * GammaVl(fKs, mKs, Vus, M, mmu)
    else:
        z = 0
    return z


def GammaPi0(M, Ue, Umu, Utau):
    mPi0 = 0.13498
    fPi = 0.13
    if M >= mPi0:
        z = GammaPnu(fPi, mPi0, M, Ue, Umu, Utau)
    else:
        z = 0
    return z


def GammaEta(M, Ue, Umu, Utau):
    mEta = 0.55
    fEta = 0.08
    if M >= mEta:
        z = GammaPnu(fEta, mEta, M, Ue, Umu, Utau)
    else:
        z = 0
    return z


def GammaEtap(M, Ue, Umu, Utau):
    mEtap = 0.96
    fEtap = -0.0946
    if M >= mEtap:
        z = GammaPnu(fEtap, mEtap, M, Ue, Umu, Utau)
    else:
        z = 0
    return z


def GammaPie(M, Ue, Umu, Utau):
    mPi = 0.13957
    fPi = 0.13
    Vud = 9.737e-1
    me = 5.11e-4
    if M >= (mPi + me):
        z = (Ue**2) * GammaPl(fPi, mPi, Vud, M, me)
    else:
        z = 0
    return z


def GammaPimu(M, Ue, Umu, Utau):
    mPi = 0.13957
    fPi = 0.13
    Vud = 9.737e-1
    mmu = 1.05e-1
    if M >= (mPi + mmu):
        z = (Umu**2) * GammaPl(fPi, mPi, Vud, M, mmu)
    else:
        z = 0
    return z


def GammaPitau(M, Ue, Umu, Utau):
    mPi = 0.13957
    fPi = 0.13
    Vud = 9.737e-1
    mtau = 1.777
    if M >= (mPi + mtau):
        z = (Utau**2) * GammaPl(fPi, mPi, Vud, M, mtau)
    else:
        z = 0
    return z


def GammaKe(M, Ue, Umu, Utau):
    mK = 0.493
    fK = 0.156
    Vus = 2.243e-1
    me = 5.11e-4
    if M >= (mK + me):
        z = (Ue**2) * GammaPl(fK, mK, Vus, M, me)
    else:
        z = 0
    return z


def GammaKmu(M, Ue, Umu, Utau):
    mK = 0.493
    fK = 0.156
    Vus = 2.243e-1
    mmu = 1.05e-1
    if M >= (mK + mmu):
        z = (Umu**2) * GammaPl(fK, mK, Vus, M, mmu)
    else:
        z = 0
    return z


def GammaDe(M, Ue, Umu, Utau):
    mD = 1.896
    fD = 0.212
    Vcd = 0.225
    me = 5.11e-4
    if M >= (mD + me):
        z = (Ue**2) * GammaPl(fD, mD, Vcd, M, me)
    else:
        z = 0
    return z


def GammaDmu(M, Ue, Umu, Utau):
    mD = 1.896
    fD = 0.212
    Vcd = 0.225
    mmu = 1.05e-1
    if M >= (mD + mmu):
        z = (Umu**2) * GammaPl(fD, mD, Vcd, M, mmu)
    else:
        z = 0
    return z


def GammaDse(M, Ue, Umu, Utau):
    mDs = 1.97
    fDs = 0.249
    Vcs = 0.97
    me = 5.11e-4
    if M >= (mDs + me):
        z = (Ue**2) * GammaPl(fDs, mDs, Vcs, M, me)
    else:
        z = 0
    return z


def GammaDsmu(M, Ue, Umu, Utau):
    mDs = 1.97
    fDs = 0.249
    Vcs = 0.97
    mmu = 1.05e-1
    if M >= (mDs + mmu):
        z = (Umu**2) * GammaPl(fDs, mDs, Vcs, M, mmu)
    else:
        z = 0
    return z


def Gamma1mese(M, Ue, Umu, Utau):
    z = (
        GammaRhoe(M, Ue, Umu, Utau)
        + GammaPie(M, Ue, Umu, Utau)
        + GammaKe(M, Ue, Umu, Utau)
        + GammaKse(M, Ue, Umu, Utau)
        + GammaDe(M, Ue, Umu, Utau)
        + GammaDse(M, Ue, Umu, Utau)
    )
    return z


def Gamma1mesmu(M, Ue, Umu, Utau):
    z = (
        GammaRhomu(M, Ue, Umu, Utau)
        + GammaPimu(M, Ue, Umu, Utau)
        + GammaKmu(M, Ue, Umu, Utau)
        + GammaKsmu(M, Ue, Umu, Utau)
        + GammaDmu(M, Ue, Umu, Utau)
        + GammaDsmu(M, Ue, Umu, Utau)
    )
    return z


def Gamma1mesv(M, Ue, Umu, Utau):
    z = (
        GammaPhi(M, Ue, Umu, Utau)
        + GammaRho(M, Ue, Umu, Utau)
        + GammaOmega(M, Ue, Umu, Utau)
        + GammaPi0(M, Ue, Umu, Utau)
        + GammaEta(M, Ue, Umu, Utau)
        + GammaEtap(M, Ue, Umu, Utau)
    )
    return z


def Gamma1mes(M, Ue, Umu, Utau):
    z = Gamma1mese(M, Ue, Umu, Utau) + Gamma1mesmu(M, Ue, Umu, Utau) + Gamma1mesv(M, Ue, Umu, Utau) + GammaPitau(M, Ue, Umu, Utau)
    return z


def GammaPNl(f, mp, V, mN, ml):
    GF = 1.16e-5
    if mp >= (mN + ml):
        z = (
            (GF**2)
            * (mp**3)
            * (f**2)
            * (V**2)
            / (8 * pi)
            * sqrt(l(1, xfrac(mN, mp) ** 2, xfrac(ml, mp) ** 2))
            * (xfrac(mN, mp) ** 2 + xfrac(ml, mp) ** 2 - (xfrac(mN, mp) ** 2 - xfrac(ml, mp) ** 2) ** 2)
        )
    else:
        z = 0
    return z


def GammaPiNe(mN):
    mPi = 0.13957
    fPi = 0.13
    Vud = 9.737e-1
    me = 5.11e-4
    z = GammaPNl(fPi, mPi, Vud, mN, me)
    return z


def GammaPiNmu(mN):
    mPi = 0.13957
    fPi = 0.13
    Vud = 9.737e-1
    mmu = 1.05e-1
    z = GammaPNl(fPi, mPi, Vud, mN, mmu)
    return z


def GammaKNe(mN):
    mK = 0.493
    fK = 0.156
    Vus = 2.243e-1
    me = 5.11e-4
    z = GammaPNl(fK, mK, Vus, mN, me)
    return z


def GammaKNmu(mN):
    mK = 0.493
    fK = 0.156
    Vus = 2.243e-1
    mmu = 1.05e-1
    z = GammaPNl(fK, mK, Vus, mN, mmu)
    return z


def GammaDNe(mN):
    mD = 1.896
    fD = 0.212
    Vcd = 0.225
    me = 5.11e-4
    z = GammaPNl(fD, mD, Vcd, mN, me)
    return z


def GammaDNmu(mN):
    mD = 1.896
    fD = 0.212
    Vcd = 0.225
    mmu = 1.05e-1
    z = GammaPNl(fD, mD, Vcd, mN, mmu)
    return z


def GammaDsNe(mN):
    mDs = 1.97
    fDs = 0.249
    Vcs = 0.97
    me = 5.11e-4
    z = GammaPNl(fDs, mDs, Vcs, mN, me)
    return z


def GammaDsNmu(mN):
    mDs = 1.97
    fDs = 0.249
    Vcs = 0.97
    mmu = 1.05e-1
    z = GammaPNl(fDs, mDs, Vcs, mN, mmu)
    return z


def GammaVNl(f, MV, V, mN, ml):
    GF = 1.16e-5
    if MV >= (mN + ml):
        z = (
            (GF**2)
            * (V**2)
            * (f**2)
            * (MV**3)
            / (12 * pi)
            * sqrt(l(1, xfrac(mN, MV) ** 2, xfrac(ml, MV) ** 2))
            * (2 - xfrac(mN, MV) ** 2 - xfrac(ml, MV) ** 2 - (xfrac(mN, MV) ** 2 - xfrac(ml, MV) ** 2) ** 2)
        )
    else:
        z = 0
    return z


def GammaRhoNe(mN):
    mRho = 0.775
    fRho = 0.17
    Vud = 9.737e-1
    me = 5.11e-4
    z = GammaVNl(fRho, mRho, Vud, mN, me)
    return z


def GammaRhoNmu(mN):
    mRho = 0.775
    fRho = 0.17
    Vud = 9.737e-1
    mmu = 1.05e-1
    z = GammaVNl(fRho, mRho, Vud, mN, mmu)
    return z


# R-functions for the electron flavor


def ReBorexinoCC(M):
    me = 5.11e-4
    if M >= 2 * me:
        z = sqrt(Gammavee(M, 1, 0, 0) / GammaveeCC(M, 1, 0, 0))
    else:
        z = 0
    return z


def ReT2KCC(M):
    me = 5.11e-4
    if M >= 2 * me:
        z = sqrt(
            (GammaPie(M, 1, 0, 0) + Gammavee(M, 1, 0, 0) + Gammavemmup(M, 1, 0, 0) + Gammavmumu(M, 1, 0, 0))
            / (GammaveeCC(M, 1, 0, 0) + GammaPie(M, 1, 0, 0) + Gammavemmup(M, 1, 0, 0))
        )
    else:
        z = 0
    return z


def ReCHARMCC(M):
    me = 5.11e-4
    if M >= 2 * me:
        z = sqrt((Gammavee(M, 1, 0, 0) + Gammavemmup(M, 1, 0, 0) + Gammavmumu(M, 1, 0, 0)) / (GammaveeCC(M, 1, 0, 0) + Gammavemmup(M, 1, 0, 0)))
    else:
        z = 0
    return z


# BEBC blast, contains NC and CC.


def ReBEBCCC(M):
    mDs = 1.97
    mPi = 0.13957
    me = 5.11e-4
    if M > (mPi + me):
        if M < (mDs - me):
            z = sqrt(
                (GammaPie(M, 1, 0, 0) + Gammavee(M, 1, 0, 0) + Gammavmumu(M, 1, 0, 0) + Gammavemmup(M, 1, 0, 0))
                / (GammaPie(M, 1, 0, 0) + GammaveeCC(M, 1, 0, 0) + Gammavemmup(M, 1, 0, 0))
            )
    return z


def ReLHCCC(M):
    me = 5.11e-4
    if M >= 2 * me:
        z = sqrt((Gammavee(M, 1, 0, 0) + Gammavemmup(M, 1, 0, 0) + Gammavmumu(M, 1, 0, 0)) / (GammaveeCC(M, 1, 0, 0) + Gammavemmup(M, 1, 0, 0)))
    else:
        z = 0
    return z


# R-functions for the muon flavor


def RmuNUTEVCC(M):
    me = 5.11e-4
    mmu = 1.05e-1
    if M >= me + mmu:
        z = sqrt(
            (Gammavmumep(M, 0, 1, 0) + Gammavmumu(M, 0, 1, 0) + GammaPimu(M, 0, 1, 0) + GammaRhomu(M, 0, 1, 0))
            / (GammavmumuCC(M, 0, 1, 0) + Gammavmumep(M, 0, 1, 0) + GammaPimu(M, 0, 1, 0) + GammaRhomu(M, 0, 1, 0))
        )
    else:
        z = 0
    return z


def RmuT2KCC(M):
    me = 5.11e-4
    mmu = 1.05e-1
    if M >= me + mmu:
        z = sqrt(
            (GammaPimu(M, 0, 1, 0) + Gammavmumu(M, 0, 1, 0) + Gammavmumep(M, 0, 1, 0) + Gammavee(M, 0, 1, 0))
            / (GammavmumuCC(M, 0, 1, 0) + GammaPimu(M, 0, 1, 0) + Gammavmumep(M, 0, 1, 0))
        )
    else:
        z = 0
    return z


def RmuLHCCC(M):
    me = 5.11e-4
    mmu = 1.05e-1
    if M >= me + mmu:
        z = sqrt((Gammavmumu(M, 0, 1, 0) + Gammavmumep(M, 0, 1, 0) + Gammavee(M, 0, 1, 0)) / (GammavmumuCC(M, 0, 1, 0) + Gammavmumep(M, 0, 1, 0)))
    else:
        z = 0
    return z
