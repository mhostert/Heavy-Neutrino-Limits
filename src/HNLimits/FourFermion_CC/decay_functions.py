#!/usr/bin/env python
# coding: utf-8

# ## Functions for the Limits on 4-fermion charged current operators $-$ $\psi^4$ (duNe)
# 
# Specific functions fot the HNL decay widths needed fot the re-scaling of the standard mixing experimental bounds.
# 

# In[94]:


import math
from math import sqrt, sin, asin, cos, tan, atan2, pi, log


# In[95]:


import mpmath as mp
from mpmath import *
mp.dps=100
mp.pretty=True


# In[96]:


def l(a,b,c):
    z = a**2+b**2+c**2-2*a*b-2*a*c-2*b*c
    return z

def xfrac(m,M):
    z = mpf(m)/mpf(M)
    return z


# In[97]:


# l(1/7,10,3)
# xfrac(1,18)


# In[98]:


def L(x):
    z= log((1 - 3*x**2 - (1 - x**2)*sqrt(1 - 4*x**2))/(x**2*(1 + sqrt(1 - 4*x**2))))
    return z


# In[99]:


# L(1/0.9)


# In[100]:


def Gammavvv(M, Ue, Umu, Utau):
    GF = 1.16e-5
    A = (GF**2)/(192*(pi**3))
    z = (Ue**2 + Umu**2 + Utau**2)*A*(M**5)
    return z


# In[101]:


# Gammavvv(1/7,0,0,1)


# In[102]:


def f1(x):
    z=(1 - 14*(x**2) - 2*(x**4) - 12*(x**6))*sqrt(1 - 4*(x**2)) + 12*(x**4)*((x**4) - 1)*L(x)
    return z

def f2(x):
    z=4*((x**2)*(2 + 10*(x**2) - 12*(x**4))*sqrt(1 - 4*(x**2)) + 6*(x**4)*(1 - 2*(x**2) + 2*(x**4))*L(x))
    return z

def Gammavll1(ml, M):
    GF = 1.16e-5
    A = (GF**2)/(192*(pi**3))
    sw2 = 2.3e-1;
    C1lep = 0.25*(1 - 4*sw2 + 8*(sw2**2));
    C2lep = 0.5*(2*(sw2**2) - sw2);
    z= A*(M**5)*((C1lep*f1(xfrac(ml, M)) + C2lep*f2(xfrac(ml, M))))
    return z

def Gammavll2(ml, M):
    GF = 1.16e-5
    A = (GF**2)/(192*(pi**3))
    sw2 = 2.3e-1
    C1lep = 0.25*(1 - 4*sw2 + 8*(sw2**2))
    C2lep = 0.5*(2*(sw2**2) - sw2)
    z= A*(M**5)*sw2*(f2(xfrac(ml, M)) + 2*f1(xfrac(ml, M)))
    return z

def Gammavee(M, Ue, Umu, Utau):
    me = 5.11e-4
    if M >= 2*me:
        z= (Ue**2 + Umu**2 + Utau**2)*Gammavll1(me,M) + (Ue**2)*Gammavll2(me, M)
    else:
        z=0
    return z

def Gammavmumu(M, Ue, Umu, Utau):
    mmu = 1.05e-1
    if M >= 2*mmu:
        z= (Ue**2 + Umu**2 + Utau**2)*Gammavll1(mmu,M) + (Umu**2)*Gammavll2(mmu, M)
    else:
        z=0
    return z


# In[103]:


def Gammavllp(M, m):
    GF = 1.16e-5
    A = (GF**2)/(192*(pi**3))
    z= A*(M**5)*(1 - 8*(xfrac(m, M)**2) + 8*(xfrac(m, M)**6) - xfrac(m, M)**8 - 12*(xfrac(m, M)**4)*log(xfrac(m, M)**2))
    return z

def Gammavemmup(M, Ue, Umu, Utau):
    me = 5.11e-4
    mmu = 1.05e-1
    if M >= (me + mmu):
        z= (Ue**2)*Gammavllp(M, mmu) 
    else:
        z=0
    return z

def Gammavmumep(M, Ue, Umu, Utau):
    me = 5.11e-4
    mmu = 1.05e-1
    if M >= (me + mmu):
        z= (Umu**2)*Gammavllp(M, mmu) 
    else:
        z=0
    return z

def Gammavemtaup(M, Ue, Umu, Utau):
    me = 5.11e-4
    mtau = 1.777
    if M >= (me + mtau):
        z= (Ue**2)*Gammavllp(M, mtau) 
    else:
        z=0
    return z

def Gammavtaumep(M, Ue, Umu, Utau):
    me = 5.11e-4
    mtau = 1.777
    if M >= (me + mtau):
        z= (Utau**2)*Gammavllp(M, mtau) 
    else:
        z=0
    return z


# In[104]:


# Gammavtaumep(80,1,1,1)


# In[105]:


def Gammalep(M, Ue, Umu, Utau):
    z = Gammavvv(M, Ue, Umu, Utau) + Gammavee(M, Ue, Umu, Utau) + Gammavmumu(M, Ue, Umu, Utau) + Gammavemmup(M, Ue, Umu, Utau) + Gammavmumep(M, Ue, Umu, Utau) + Gammavemtaup(M, Ue, Umu, Utau) + Gammavtaumep(M, Ue, Umu, Utau)
    return z


# In[106]:


# Gammalep(120,1e-3,1,1)


# In[107]:


from scipy import integrate


# In[108]:


def inte(s, x, y, z):
    z= 12*(1/s*(s - x - y)*(1 + z - s)*sqrt(l(s, x, y)*l(1, s, z)))
    return z

def IntCC(x, y, z):
    z=integrate.quad(inte, (sqrt(x) + sqrt(y))**2, (1 - sqrt(z))**2, args=(x,y,z))
    return z[0]

def GammaveeCC(M, Ue, Umu, Utau):
    me = 5.11e-4
    GF = 1.16e-5
    A = (GF**2)/(192*(pi**3))
    if M >= 2*me:
        z=A*(M**5)*(Ue**2)*IntCC(0, (xfrac(me, M))**2, (xfrac(me, M))**2)
    else:
        z=0
    return z

def GammavmumuCC(M, Ue, Umu, Utau):
    mmu = 1.05e-1
    GF = 1.16e-5
    A = (GF**2)/(192*(pi**3))
    if M >= 2*mmu:
        z=A*(M**5)*(Umu**2)*IntCC(0, (xfrac(mmu, M))**2, (xfrac(mmu, M))**2)
    else:
        z=0
    return z


# In[109]:


# GammavmumuCC(22/100,1,1,1)


# In[110]:


def GammaVl(f, mv, V, M, ml):
    GF = 1.16e-5
    z= (GF**2)*(M**3)/(16*pi*(mv**2))*(f**2)*(V**2)*sqrt((1 + ((mv**2) - (ml**2))/(M**2))**2 - 4*(mv**2)/(M**2))*((1 - (mv**2)/(M**2))*(1 + 2*(mv**2)/(M**2)) + (ml**2)/(M**2)*(((mv**2) + (ml**2))/(M**2) - 2))
    return z

def GammaVnu(f, mv, g, M, Ue, Umu, Utau):
    GF = 1.16e-5
    z = (Ue**2 + Umu**2 + Utau**2)*(GF**2)*(M**3)/(32*pi*(mv**2))*(f**2)*(g**2)*(1 + 2*(mv**2)/(M**2))*(1 - (mv**2)/(M**2))**2
    return z
    
def GammaPl(f, mp, V, M, ml):
    GF = 1.16e-5
    z = (GF**2)*(M**3)/(16*pi)*(f**2)*(V**2)*sqrt((1 + ((mp**2) - (ml**2))/(M**2))**2 - 4*(mp**2)/(M**2))*(1 - (mp**2)/(M**2) - (ml**2)/(M**2)*(2 + ((mp**2) - (ml**2))/(M**2)))
    return z
    
def GammaPnu(f, mp, M, Ue, Umu, Utau):
    GF = 1.16e-5
    z = (Ue**2 + Umu**2 + Utau**2)*(GF**2)*(M**3)/(32*pi)*(f**2)*(1 - (mp**2)/(M**2))**2
    return z


# In[111]:


def GammaPhi(M,Ue,Umu,Utau):
    mPhi = 1.02
    fPhi = 0.232
    sw2 = 2.3e-1
    gPhi = sqrt(2)*(0.5 - 2*sw2/3)
    if M>= mPhi:
        z= GammaVnu(fPhi, mPhi, gPhi, M, Ue, Umu, Utau)
    else:
        z= 0
    return z

def GammaRho(M,Ue,Umu,Utau):
    mRho = 0.775
    fRho = 0.17
    sw2 = 2.3e-1
    gRho = 1 - 2*sw2
    if M>= mRho:
        z= GammaVnu(fRho, mRho, gRho, M, Ue, Umu, Utau)
    else:
        z= 0
    return z

def GammaOmega(M,Ue,Umu,Utau):
    mOmega = 0.782
    fOmega = 0.155
    sw2 = 2.3e-1
    gOmega = 2*sw2/3
    if M>= mOmega:
        z= GammaVnu(fOmega, mOmega, gOmega, M, Ue, Umu, Utau)
    else:
        z= 0
    return z


def GammaRhoe(M, Ue, Umu, Utau):
    mRho = 0.775
    fRho = 0.17
    Vud = 9.737e-1
    me = 5.11e-4
    if M>= (mRho + me):
        z = (Ue**2)*GammaVl(fRho, mRho, Vud, M, me)
    else:
        z = 0
    return z

def GammaRhomu(M, Ue, Umu, Utau):
    mRho = 0.775
    fRho = 0.17
    Vud = 9.737e-1
    mmu = 1.05e-1
    if M>= (mRho + mmu):
        z = (Umu**2)*GammaVl(fRho, mRho, Vud, M, mmu)
    else:
        z = 0
    return z

def GammaKse(M, Ue, Umu, Utau):
    mKs = 0.892
    fKs = 0.177
    Vus = 2.243e-1
    me = 5.11e-4
    if M>= (mKs + me):
        z = (Ue**2)*GammaVl(fKs, mKs, Vus, M, me)
    else:
        z = 0
    return z

def GammaKsmu(M, Ue, Umu, Utau):
    mKs = 0.892
    fKs = 0.177
    Vus = 2.243e-1
    mmu = 1.05e-1
    if M>= (mKs + mmu):
        z = (Umu**2)*GammaVl(fKs, mKs, Vus, M, mmu)
    else:
        z = 0
    return z

def GammaPi0(M, Ue, Umu, Utau):
    mPi0 = 0.13498
    fPi = 0.13
    if M>= mPi0:
        z= GammaPnu(fPi, mPi0, M, Ue, Umu, Utau)
    else:
        z=0
    return z

def GammaEta(M, Ue, Umu, Utau):
    mEta = 0.55
    fEta = 0.08
    if M>= mEta:
        z= GammaPnu(fEta, mEta, M, Ue, Umu, Utau)
    else:
        z=0
    return z

def GammaEtap(M, Ue, Umu, Utau):
    mEtap = 0.96
    fEtap = -0.0946
    if M>= mEtap:
        z= GammaPnu(fEtap, mEtap, M, Ue, Umu, Utau)
    else:
        z=0
    return z

def GammaPie(M, Ue, Umu, Utau):
    mPi = 0.13957
    fPi = 0.13
    Vud = 9.737e-1
    me = 5.11e-4
    if M >= (mPi + me):
        z = (Ue**2)*GammaPl(fPi, mPi, Vud, M, me)
    else:
        z = 0
    return z

def GammaPimu(M, Ue, Umu, Utau):
    mPi = 0.13957
    fPi = 0.13
    Vud = 9.737e-1
    mmu = 1.05e-1
    if M >= (mPi + mmu):
        z = (Umu**2)*GammaPl(fPi, mPi, Vud, M, mmu)
    else:
        z = 0
    return z

def GammaPitau(M, Ue, Umu, Utau):
    mPi = 0.13957
    fPi = 0.13
    Vud = 9.737e-1
    mtau = 1.777
    if M >= (mPi + mtau):
        z = (Utau**2)*GammaPl(fPi, mPi, Vud, M, mtau)
    else:
        z = 0
    return z

def GammaKe(M, Ue, Umu, Utau):
    mK = 0.493
    fK = 0.156
    Vus = 2.243e-1
    me = 5.11e-4
    if M >= (mK + me):
        z = (Ue**2)*GammaPl(fK, mK, Vus, M, me)
    else:
        z = 0
    return z

def GammaKmu(M, Ue, Umu, Utau):
    mK = 0.493
    fK = 0.156
    Vus = 2.243e-1
    mmu = 1.05e-1
    if M >= (mK + mmu):
        z = (Umu**2)*GammaPl(fK, mK, Vus, M, mmu)
    else:
        z = 0
    return z

def GammaDe(M, Ue, Umu, Utau):
    mD = 1.896
    fD = 0.212
    Vcd = 0.225
    me = 5.11e-4
    if M >= (mD + me):
        z = (Ue**2)*GammaPl(fD, mD, Vcd, M, me)
    else:
        z = 0
    return z

def GammaDmu(M, Ue, Umu, Utau):
    mD = 1.896
    fD = 0.212
    Vcd = 0.225
    mmu = 1.05e-1
    if M >= (mD + mmu):
        z = (Umu**2)*GammaPl(fD, mD, Vcd, M, mmu)
    else:
        z = 0
    return z

def GammaDse(M, Ue, Umu, Utau):
    mDs = 1.97
    fDs = 0.249
    Vcs = 0.97
    me = 5.11e-4
    if M >= (mDs + me):
        z = (Ue**2)*GammaPl(fDs, mDs, Vcs, M, me)
    else:
        z = 0
    return z

def GammaDsmu(M, Ue, Umu, Utau):
    mDs = 1.97
    fDs = 0.249
    Vcs = 0.97
    mmu = 1.05e-1
    if M >= (mDs + mmu):
        z = (Umu**2)*GammaPl(fDs, mDs, Vcs, M, mmu)
    else:
        z = 0
    return z


# In[112]:


def Gamma1mese(M, Ue, Umu, Utau):
    z = GammaRhoe(M, Ue, Umu, Utau) + GammaPie(M, Ue, Umu, Utau) + GammaKe(M, Ue, Umu, Utau) + GammaKse(M, Ue, Umu, Utau) + GammaDe(M, Ue, Umu, Utau) + GammaDse(M, Ue, Umu, Utau)
    return z

def Gamma1mesmu(M, Ue, Umu, Utau):
    z = GammaRhomu(M, Ue, Umu, Utau) + GammaPimu(M, Ue, Umu, Utau) + GammaKmu(M, Ue, Umu, Utau) + GammaKsmu(M, Ue, Umu, Utau) + GammaDmu(M, Ue, Umu, Utau) + GammaDsmu(M, Ue, Umu, Utau)
    return z

def Gamma1mesv(M, Ue, Umu, Utau):
    z = GammaPhi(M, Ue, Umu, Utau) + GammaRho(M, Ue, Umu, Utau) + GammaOmega(M, Ue, Umu, Utau) + GammaPi0(M, Ue, Umu, Utau) + GammaEta(M, Ue, Umu, Utau) + GammaEtap(M, Ue, Umu, Utau)
    return z

def Gamma1mes(M, Ue, Umu, Utau):
    z = Gamma1mese(M, Ue, Umu, Utau) + Gamma1mesmu(M, Ue, Umu, Utau) + Gamma1mesv(M, Ue, Umu, Utau) + GammaPitau(M, Ue, Umu, Utau) 
    return z


# In[113]:


# Gamma1mes(200,1,1/5,1/3)


# In[114]:


def GammaPNl(f, mp, V, mN, ml):
    GF = 1.16e-5
    if mp>= (mN+ml):
        z =  (GF**2)*(mp**3)*(f**2)*(V**2)/(8*pi)*sqrt(l(1, xfrac(mN, mp)**2, xfrac(ml, mp)**2))*(xfrac(mN, mp)**2 + xfrac(ml, mp)**2 - (xfrac(mN, mp)**2 - xfrac(ml, mp)**2)**2)
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

def GammaBNe(mN):
    mB = 5.279
    fB = 0.19
    Vub = 3.82e-3
    me = 5.11e-4
    z = GammaPNl(fB, mB, Vub, mN, me)
    return z

def GammaBNmu(mN):
    mB = 5.279
    fB = 0.19
    Vub = 3.82e-3
    mmu = 1.05e-1
    z = GammaPNl(fB, mB, Vub, mN, mmu)
    return z


# In[115]:


# GammaDsNmu(1/100)


# In[116]:


def GammaVNl(f, MV, V, mN, ml):
    GF = 1.16e-5
    if MV >= (mN + ml):
        z = (GF**2)*(V**2)*(f**2)*(MV**3)/(12*pi)*sqrt(l(1, xfrac(mN, MV)**2, xfrac(ml, MV)**2))*(2 - xfrac(mN, MV)**2 - xfrac(ml, MV)**2 - (xfrac(mN, MV)**2 - xfrac(ml, MV)**2)**2)
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


# In[117]:


# Decay of HNL into mesons

def GammaPieX(mN):
    Vud = 9.737e-1
    GF = 1.16e-5
    z = GammaPie(mN,1,0,0)/((Vud**2)*16*(GF**2))
    return z

def GammaPimuX(mN):
    Vud = 9.737e-1
    GF = 1.16e-5
    z = GammaPimu(mN,0,1,0)/((Vud**2)*16*(GF**2))
    return z

def GammaKeX(mN):
    Vus = 2.243e-1
    GF = 1.16e-5
    z = GammaKe(mN,1,0,0)/((Vus**2)*16*(GF**2))
    return z

def GammaKmuX(mN):
    Vus = 2.243e-1
    GF = 1.16e-5
    z = GammaKmu(mN,0,1,0)/((Vus**2)*16*(GF**2))
    return z

def GammaDmuX(mN):
    Vcd = 0.225
    GF = 1.16e-5
    z = GammaDmu(mN,0,1,0)/((Vcd**2)*16*(GF**2))
    return z

def GammaDsmuX(mN):
    Vcs = 0.97
    GF = 1.16e-5
    z = GammaDsmu(mN,0,1,0)/((Vcs**2)*16*(GF**2))
    return z

def GammaRhomuX(mN):
    Vud = 9.737e-1
    GF = 1.16e-5
    z = GammaRhomu(mN,0,1,0)/((Vud**2)*16*(GF**2))
    return z


# In[118]:


# Decay of mesons into HNL

def GammaPiNeX(mN):
    Vud = 9.737e-1
    GF = 1.16e-5
    z = GammaPiNe(mN)/((Vud**2)*16*(GF**2))
    return z

def GammaPiNmuX(mN):
    Vud = 9.737e-1
    GF = 1.16e-5
    z = GammaPiNmu(mN)/((Vud**2)*16*(GF**2))
    return z

def GammaKNeX(mN):
    Vus = 2.243e-1
    GF = 1.16e-5
    z = GammaKNe(mN)/((Vus**2)*16*(GF**2))
    return z

def GammaKNmuX(mN):
    Vus = 2.243e-1
    GF = 1.16e-5
    z = GammaKNmu(mN)/((Vus**2)*16*(GF**2))
    return z

def GammaDNeX(mN):
    Vcd = 0.225
    GF = 1.16e-5
    z = GammaDNe(mN)/((Vcd**2)*16*(GF**2))
    return z

def GammaDNmuX(mN):
    Vcd = 0.225
    GF = 1.16e-5
    z = GammaDNmu(mN)/((Vcd**2)*16*(GF**2))
    return z

def GammaDsNeX(mN):
    Vcs = 0.97
    GF = 1.16e-5
    z = GammaDsNe(mN)/((Vcs**2)*16*(GF**2))
    return z

def GammaDsNmuX(mN):
    Vcs = 0.97
    GF = 1.16e-5
    z = GammaDsNmu(mN)/((Vcs**2)*16*(GF**2))
    return z

def GammaRhoNmuX(mN):
    Vud = 9.737e-1
    GF = 1.16e-5
    z = GammaRhoNmu(mN)/((Vud**2)*16*(GF**2))
    return z

def GammaBNeX(mN):
    Vub = 3.82e-3
    GF = 1.16e-5
    z = GammaBNe(mN)/((Vub**2)*16*(GF**2))
    return z

def GammaBNmuX(mN):
    Vub = 3.82e-3
    GF = 1.16e-5
    z = GammaBNmu(mN)/((Vub**2)*16*(GF**2))
    return z

# ##### Semi-leptonic decays: fix this!

# def GammaDK0NmuX(mN):
#     Vcs = 0.97
#     GF = 1.16e-5
#     z = 1/((Vcs**2)*16*(GF**2)) #We should include the semi-leptonic decay of the D. I am approximating by the same "standard" width! 
#     return z

# def GammaBDNeX(mN):
#     Vcb = 40.8e-3
#     GF = 1.16e-5
#     z = 1/((Vcb**2)*16*(GF**2)) #We should include the semi-leptonic decay of the B. I am approximating by the same "standard" width! 
#     return z

# def GammaBDNmuX(mN):
#     Vcb = 40.8e-3
#     GF = 1.16e-5
#     z = 1/((Vcb**2)*16*(GF**2)) #We should include the semi-leptonic decay of the B. I am approximating by the same "standard" width! 
#     return z


# In[119]:


def Zud(FB):
    Vud = 9.737e-1
    if FB == True:
        z = 1
    else:
        z = Vud
    return z

def Zcd(FB):
    Vcd = 0.225
    if FB == True:
        z = 1
    else:
        z = Vcd
    return z

def Zus(FB):
    Vus = 2.243e-1
    if FB == True:
        z = 1
    else:
        z = Vus
    return z

def Zcs(FB):
    Vcs = 0.97
    if FB == True:
        z = 1
    else:
        z = Vcs
    return z

def Zub(FB):
    Vub = 3.82e-3
    if FB == True:
        z = 1
    else:
        z = Vub
    return z

def Zcb(FB):
    Vcb = 40.8e-3
    if FB == True:
        z = 1
    else:
        z = Vcb
    return z


# In[ ]:


#### electron case ####

# Peak searches. Kinematics of the meson decay products.

def RePSPion4fX(M,FB):
    z = (GammaPiNe(M)/GammaPiNeX(M))/(Zud(FB)**2) # re-scaling of U**2
    return z

def RePSKaon4fX(M,FB):
    z = (GammaKNe(M)/GammaKNeX(M))/(Zus(FB)**2) # re-scaling of U**2
    return z

# Decay on flight of the HNL into pion. The HNL come from K leptonic decay.

def ReT2K4fX(M,FB):
    mPi = 0.13957
    me = 5.11e-4
    if M > (mPi-me):
        z = sqrt(GammaKNe(M)*(GammaPie(M,1,0,0)+Gammavee(M,1,0,0)+Gammavmumu(M,1,0,0)+Gammavemmup(M,1,0,0))/(GammaKNeX(M)*GammaPieX(M)))/(Zus(FB)*Zud(FB))
    return z

# Re-analysis of BEBC. D and Ds (they do not say it, but the mass range of the HNL needs a Ds) leptonic decay. HNL decay to leptons, and pion channel. 

def ReBEBC4fX(M,FB):
    mDs = 1.97
    mPi = 0.13957
    me = 5.11e-4
    if M > (mPi+me):
        if M < (mDs-me):
            z = sqrt((GammaDNe(M)+GammaDsNe(M))*(GammaPie(M,1,0,0)+Gammavee(M,1,0,0)+Gammavmumu(M,1,0,0)+Gammavemmup(M,1,0,0))/((GammaDNeX(M)*Zcd(FB)+GammaDsNeX(M)*Zcs(FB))*GammaPieX(M)))/Zud(FB)
    # else:
    #     z = 1
    return z

def ReBEBC4fXlow(M,FB):
    # mD = 1.896
    mDs = 1.97
    mPi = 0.13957
    me = 5.11e-4
    if M > (mPi+me):
        if M < (mDs-me):
            z = (GammaPie(M,1,0,0)+Gammavee(M,1,0,0)+Gammavmumu(M,1,0,0)+Gammavemmup(M,1,0,0))/GammaPieX(M)/(Zud(FB)**2) ## the upper line goes with the width
    # else:
    #     z = 1
    return z

# Belle. B to XlN, with X = D,D*,pi,rho..., nothing and N > pi e. We approximate the width of the semileptonic B decay to that of the standard mixing. Fix this!

# def ReBelle4fYZ(M,FB):
#     mB = 5.279
#     mPi = 0.13957
#     me = 5.11e-4 
#     if M > (mPi+me):
#         if M < (mB-me):
#             #z = sqrt(GammaPie(M,1,0,0)/(GammaPieYZ(M)))/(Zcb(FB)*Zud(FB)) #+ 
#             z = sqrt(GammaBNe(M)*GammaPie(M,1,0,0)/(GammaBNeYZ(M)*GammaPieYZ(M)))/(Zub(FB)*Zud(FB)) # re-scaling of U**2
#     else:
#         z = 1
#     return z

# def ReBelle4fYZlow(M,FB):
#     mB = 5.279
#     mPi = 0.13957
#     me = 5.11e-4
#     if M > (mPi+me):
#         if M < (mB-me):
#             z = GammaPie(M,1,0,0)/GammaPieYZ(M)/(Zud(FB)**2) # re-scaling of U**2
#     else:
#         z = 1
#     return z

#### muon case ####

# Peak searches. Kinematics of the meson decay products.

def RmuPSPion4fX(M,FB):
    z = (GammaPiNmu(M)/GammaPiNmuX(M))/(Zud(FB)**2) # re-scaling of U**2
    return z

def RmuPSKaon4fX(M,FB):
    z = (GammaKNmu(M)/GammaKNmuX(M))/(Zus(FB)**2) # re-scaling of U**2
    return z

# Decay on flight of the HNL into pion. The HNL come from K leptonic decay.

def RmuT2K4fX(M,FB):
    mPi = 0.13957
    mmu = 1.05e-1
    if M > (mPi+mmu):
        z = sqrt((GammaKNmu(M))*(GammaPimu(M,0,1,0)+Gammavee(M,0,1,0)+Gammavmumu(M,0,1,0)+Gammavmumep(M,0,1,0))/(GammaKNmuX(M)*GammaPimuX(M)))/(Zus(FB)*Zud(FB))
    else:
        z = 1
    return z

# NuTeV: Kaon, D and Ds to mu HNL, HNL to mu pi, mu rho, and leptonic channels. Only CC channels!

def RmuNuTeV4fX(M,FB):
    mDs = 1.97
    mPi = 0.13957
    mmu = 1.05e-1
    if M > (mPi+mmu):
        if M < (mDs-mmu):
            z = sqrt((GammaKNmu(M)+GammaDNmu(M)+GammaDsNmu(M))*(GammaRhomu(M,0,1,0)+GammaPimu(M,0,1,0)+GammavmumuCC(M,0,1,0)+Gammavmumep(M,0,1,0))/((GammaKNmuX(M)*(Zus(FB)**2)+GammaDNmuX(M)*(Zcd(FB)**2)+GammaDsNmuX(M)*(Zcs(FB)**2))*(GammaPimuX(M)+GammaRhomuX(M))))/Zud(FB)
    else:
        z = 1
    return z

def RmuNuTeV4fXlow(M,FB):
    mDs = 1.97
    mPi = 0.13957
    mmu = 1.05e-1
    if M > (mPi+mmu):
        if M < (mDs-mmu):
            z = (GammaPimu(M,0,1,0)+GammavmumuCC(M,0,1,0)+Gammavmumep(M,0,1,0))/GammaPimuX(M)/(Zud(FB)**2)
    else:
        z = 1
    return z

# BEBC original: D to mu HNL, HNL to pi mu, leptons but CC only!

def RmuBEBC4fX(M,FB):
    mPi = 0.13957
    mmu = 1.05e-1 
    if M > (mPi+mmu):
        z = sqrt(GammaDNmu(M)*(GammaPimu(M,0,1,0)+GammavmumuCC(M,0,1,0)+Gammavmumep(M,0,1,0))/(GammaDNmuX(M)*(GammaPimuX(M))))/(Zcd(FB)*Zud(FB)) # re-scaling of U**2
    else:
        z = 1
    return z

def RmuBEBC4fXlow(M,FB):
    mPi = 0.13957
    mmu = 1.05e-1 
    if M > (mPi+mmu):
        z = (GammaPimu(M,0,1,0)+GammavmumuCC(M,0,1,0)+Gammavmumep(M,0,1,0))/GammaPimuX(M)/(Zud(FB)**2) # re-scaling of U**2
    else:
        z = 1
    return z

# def RmuNA34fYZ(M,FB):
#     mPi = 0.13957
#     mmu = 1.05e-1 
#     mD = 1.896
#     if M > (mPi+mmu):
#         if M >= (mD-mmu):
#             z = sqrt(GammaPimu(M,0,1,0)/(GammaPimuYZ(M)))/(Zcb(FB)*Zud(FB)) # re-scaling of U**2
#         if M < (mD-mmu):
#             z = sqrt(GammaPimu(M,0,1,0)/(GammaPimuYZ(M)))/(Zcs(FB)*Zud(FB)) 
#     else:
#         z = 1
#     return z

# def RmuNA34fYZlow(M,FB):
#     mPi = 0.13957
#     mmu = 1.05e-1
#     if M > (mPi+mmu):
#         z = GammaPimu(M,0,1,0)/GammaPimuYZ(M)/(Zud(FB)**2) # re-scaling of U**2
#     else:
#         z = 1
#     return z


# In[120]:


# ############### electron case ################

# def RePIENU4fX(M,FB):
#     z = (GammaPiNe(M)/GammaPiNeX(M))/(Zud(FB)**2) # re-scaling of U**2
#     return z

# def ReTRIUMF4fX(M,FB):
#     z = (GammaPiNe(M)/GammaPiNeX(M))/(Zud(FB)**2) # re-scaling of U**2
#     return z

# def ReNA624fX(M,FB):
#     z = (GammaKNe(M)/GammaKNeX(M))/(Zus(FB)**2) # re-scaling of U**2
#     return z

# def ReKENU4fX(M,FB):
#     z = (GammaKNe(M)/GammaKNeX(M))/(Zus(FB)**2) # re-scaling of U**2
#     return z

# # def ReKENU4fX2(M,FB):
# #     z = ((GammaKNe(M)/GammaKNmu(M))/(GammaKNeX(M)/GammaKNmuX(M)))/(Zus(FB)**2) # re-scaling of U**2
# #     return z

# def ReT2K4fX(M,FB):
#     mPi = 0.13957
#     me = 5.11e-4 
#     if M > (mPi+me):
#         z = sqrt(GammaKNe(M)*(GammaPie(M,1,0,0)+Gammavee(M,1,0,0)+Gammavmumu(M,1,0,0)+Gammavemmup(M,1,0,0))/(GammaKNeX(M)*(GammaPieX(M))))/(Zus(FB)*Zud(FB)) # re-scaling of U**2
#     else:
#         z = 1
#     return z

# # def ReT2K4fXlow(M,FB):
# #     z = sqrt((GammaPie(M,1,0,0)+Gammavee(M,1,0,0)+Gammavmumu(M,1,0,0)+Gammavemmup(M,1,0,0))/GammaPieX(M))/Zud(FB) # re-scaling of U**2
# #     return z

# # def ReBEBC4fX(M,FB):
# #     z = sqrt(GammaDsNe(M)*(GammaPi0(M,1,0,0)+GammaPie(M,1,0,0)+GammaPimu(M,1,0,0)+Gammavee(M,1,0,0)+Gammavmumu(M,1,0,0)+Gammavemmup(M,1,0,0))/(GammaDsNeX(M)*GammaPieX(M)))/(Zcs(FB)*Zud(FB)) # re-scaling of U**2
# #     return z

# # def ReBEBC4fXlow(M,FB):
# #     z = sqrt((GammaPi0(M,1,0,0)+GammaPie(M,1,0,0)+GammaPimu(M,1,0,0)+Gammavee(M,1,0,0)+Gammavmumu(M,1,0,0)+Gammavemmup(M,1,0,0))/GammaPieX(M))/Zud(FB) # re-scaling of U**2
# #     return z

# def ReBEBC4fX(M,FB):
#     mD = 1.896
#     mDs = 1.97
#     me = 5.11e-4
#     if M >= mD-me:    
#         z = sqrt(GammaDsNe(M)*(GammaPi0(M,1,0,0)+GammaPie(M,1,0,0)+GammaPimu(M,1,0,0)+Gammavee(M,1,0,0)+Gammavmumu(M,1,0,0)+Gammavemmup(M,1,0,0))/(GammaDsNeX(M)*GammaPieX(M)))/(Zcs(FB)*Zud(FB)) # re-scaling of U**2
#     if M < mD-me:
#         z = sqrt(GammaDNe(M)*(GammaPi0(M,1,0,0)+GammaPie(M,1,0,0)+GammaPimu(M,1,0,0)+Gammavee(M,1,0,0)+Gammavmumu(M,1,0,0)+Gammavemmup(M,1,0,0))/(GammaDNeX(M)*GammaPieX(M)))/(Zcd(FB)*Zud(FB)) + sqrt(GammaDsNe(M)*(GammaPi0(M,1,0,0)+GammaPie(M,1,0,0)+GammaPimu(M,1,0,0)+Gammavee(M,1,0,0)+Gammavmumu(M,1,0,0)+Gammavemmup(M,1,0,0))/(GammaDsNeX(M)*GammaPieX(M)))/(Zcs(FB)*Zud(FB)) # re-scaling of U**2
#     return z

# def ReBEBC4fXlow(M,FB):
#     z = (GammaPi0(M,1,0,0)+GammaPie(M,1,0,0)+GammaPimu(M,1,0,0)+Gammavee(M,1,0,0)+Gammavmumu(M,1,0,0)+Gammavemmup(M,1,0,0))/GammaPieX(M)/(Zud(FB)**2) # re-scaling of U**2
#     return z

# # def ReBelle4fX(M,FB):
# #     z = sqrt(GammaBNe(M)*(GammaPie(M,1,0,0)+GammaPimu(M,1,0,0))/(GammaBNeX(M)*GammaPieX(M)))/(Zub(FB)*Zud(FB)) # re-scaling of U**2
# #     return z

# # def ReBelle4fXlow(M,FB):
# #     z = sqrt((GammaPie(M,1,0,0)+GammaPimu(M,1,0,0))/(GammaPieX(M)))/Zud(FB) # re-scaling of U**2
# #     return z

# def ReBelle4fX(M,FB):
#     z = sqrt((GammaPie(M,1,0,0)+GammaPimu(M,1,0,0))/(GammaBDNeX(M)*GammaPieX(M)))/(Zcb(FB)*Zud(FB)) + sqrt(GammaBNe(M)*(GammaPie(M,1,0,0)+GammaPimu(M,1,0,0))/(GammaBNeX(M)*GammaPieX(M)))/(Zub(FB)*Zud(FB)) # re-scaling of U**2
#     return z

# def ReBelle4fXlow(M,FB):
#     z = (GammaPie(M,1,0,0)+GammaPimu(M,1,0,0))/(GammaPieX(M))/(Zud(FB)**2) # re-scaling of U**2
#     return z

# ############### muon case ################

# def RmuPSPion4fX(M,FB):
#     z = (GammaPiNmu(M)/GammaPiNmuX(M))/(Zud(FB)**2) # re-scaling of U**2
#     return z

# def RmuPSKaon4fX(M,FB):
#     z = (GammaKNmu(M)/GammaKNmuX(M))/(Zus(FB)**2) # re-scaling of U**2
#     return z

# def RmuT2K4fX(M,FB):
#     mPi = 0.13957
#     mmu = 1.05e-1 
#     if M > (mPi+mmu):
#         z = sqrt(GammaKNmu(M)*(GammaPimu(M,0,1,0)+Gammavee(M,0,1,0)+Gammavmumu(M,0,1,0)+Gammavemmup(M,0,1,0))/(GammaKNmuX(M)*(GammaPimuX(M))))/(Zus(FB)*Zud(FB)) # re-scaling of U**2
#     else:
#         z = 1
#     return z

# # def RmuNuTeV4fX(M,FB):
# #     mPi = 0.13957
# #     mmu = 1.05e-1 
# #     if M > (mPi+mmu):
# #         z = sqrt(GammaDNmu(M)*(GammaRhomu(M,0,1,0)+GammaPimu(M,0,1,0)+Gammavmumu(M,0,1,0)+Gammavmumep(M,0,1,0))/(GammaDNmuX(M)*(GammaPimuX(M)+GammaRhomuX(M))))/(Zcd(FB)*Zud(FB)) # re-scaling of U**2
# #     else:
# #         z = 1
# #     return z

# def RmuNuTeV4fX(M,FB):
#     mK = 0.493
#     mPi = 0.13957
#     mmu = 1.05e-1 
#     if M > (mPi+mmu):
#         if M >= (mK-mmu):
#             z = sqrt(GammaDNmu(M)*(GammaRhomu(M,0,1,0)+GammaPimu(M,0,1,0)+Gammavmumu(M,0,1,0)+Gammavmumep(M,0,1,0))/(GammaDNmuX(M)*(GammaPimuX(M)+GammaRhomuX(M))))/(Zcd(FB)*Zud(FB)) + sqrt(GammaDsNmu(M)*(GammaRhomu(M,0,1,0)+GammaPimu(M,0,1,0)+Gammavmumu(M,0,1,0)+Gammavmumep(M,0,1,0))/(GammaDsNmuX(M)*(GammaPimuX(M)+GammaRhomuX(M))))/(Zcs(FB)*Zud(FB)) # re-scaling of U**2
#         if M < (mK-mmu):
#             z = sqrt(GammaKNmu(M)*(GammaRhomu(M,0,1,0)+GammaPimu(M,0,1,0)+Gammavmumu(M,0,1,0)+Gammavmumep(M,0,1,0))/(GammaKNmuX(M)*(GammaPimuX(M)+GammaRhomuX(M))))/(Zus(FB)*Zud(FB)) + sqrt(GammaDNmu(M)*(GammaRhomu(M,0,1,0)+GammaPimu(M,0,1,0)+Gammavmumu(M,0,1,0)+Gammavmumep(M,0,1,0))/(GammaDNmuX(M)*(GammaPimuX(M)+GammaRhomuX(M))))/(Zcd(FB)*Zud(FB)) + sqrt(GammaDsNmu(M)*(GammaRhomu(M,0,1,0)+GammaPimu(M,0,1,0)+Gammavmumu(M,0,1,0)+Gammavmumep(M,0,1,0))/(GammaDsNmuX(M)*(GammaPimuX(M)+GammaRhomuX(M))))/(Zcs(FB)*Zud(FB)) # re-scaling of U**2
#     else:
#         z = 1
#     return z

# def RmuNuTeV4fXlow(M,FB):
#     mPi = 0.13957
#     mmu = 1.05e-1 
#     if M > (mPi+mmu):
#         z = (GammaRhomu(M,0,1,0)+GammaPimu(M,0,1,0)+Gammavmumu(M,0,1,0)+Gammavmumep(M,0,1,0))/(GammaPimuX(M)+GammaRhomuX(M))/(Zud(FB)**2) # re-scaling of U**2
#     else:
#         z = 1
#     return z

# # BEBC original: D to mu HNL, HNL to pi mu, leptons but CC only!

# def RmuBEBC4fX(M,FB):
#     mPi = 0.13957
#     mmu = 1.05e-1 
#     if M > (mPi+mmu):
#         z = sqrt(GammaDNmu(M)*(GammaPimu(M,0,1,0)+GammavmumuCC(M,0,1,0)+Gammavmumep(M,0,1,0))/(GammaDNmuX(M)*(GammaPimuX(M))))/(Zcd(FB)*Zud(FB)) # re-scaling of U**2
#     else:
#         z = 1
#     return z

# def RmuBEBC4fXlow(M,FB):
#     mPi = 0.13957
#     mmu = 1.05e-1 
#     if M > (mPi+mmu):
#         z = (GammaPimu(M,0,1,0)+GammavmumuCC(M,0,1,0)+Gammavmumep(M,0,1,0))/GammaPimuX(M)/(Zud(FB)**2) # re-scaling of U**2
#     else:
#         z = 1
#     return z

# # def RmuNA34fX(M,FB):
# #     mPi = 0.13957
# #     mmu = 1.05e-1 
# #     mB = 5.279
# #     mD = 1.896
# #     mDs = 1.97
# #     if M > (mPi+mmu):
# #         if M >= (mDs-mmu):
# #             z = 1
# #         #     z = sqrt(GammaBNmu(M)*GammaPimu(M,0,1,0)/(GammaBNmuX(M)*(GammaPimuX(M))))/(Zub(FB)*Zud(FB)) # re-scaling of U**2
# #         #sqrt(GammaBNmu(M)*GammaPimu(M,0,1,0)/(GammaBNmuX(M)*(GammaPimuX(M))))/(Zub(FB)*Zud(FB)) +
# #         if (M >= (mD-mmu)) & (M < (mDs-mmu)):
# #             z = sqrt(GammaDsNmu(M)*GammaPimu(M,0,1,0)/(GammaDsNmuX(M)*(GammaPimuX(M))))/(Zcs(FB)*Zud(FB)) # re-scaling of U**2
# #         if M < (mD-mmu):
# #             z = sqrt(GammaDsNmu(M)*GammaPimu(M,0,1,0)/(GammaDsNmuX(M)*(GammaPimuX(M))))/(Zcs(FB)*Zud(FB)) + sqrt(GammaDNmu(M)*GammaPimu(M,0,1,0)/(GammaDNmuX(M)*(GammaPimuX(M))))/(Zcd(FB)*Zud(FB)) 
# #     else:
# #         z = 1
# #     return z

# def RmuNA34fX(M,FB):
#     mPi = 0.13957
#     mmu = 1.05e-1 
#     #mB = 5.279
#     mD = 1.896
#     #mDs = 1.97
#     if M > (mPi+mmu):
#         if M >= (mD-mmu):
#             z = sqrt(GammaPimu(M,0,1,0)/(GammaBDNmuX(M)*(GammaPimuX(M))))/(Zcb(FB)*Zud(FB)) # re-scaling of U**2
#         if M < (mD-mmu):
#             z = sqrt(GammaPimu(M,0,1,0)/(GammaBDNmuX(M)*(GammaPimuX(M))))/(Zcb(FB)*Zud(FB)) + sqrt(GammaPimu(M,0,1,0)/(GammaDK0NmuX(M)*(GammaPimuX(M))))/(Zcs(FB)*Zud(FB)) 
#     else:
#         z = 1
#     return z

# def RmuNA34fXlow(M,FB):
#     mPi = 0.13957
#     mmu = 1.05e-1
#     if M > (mPi+mmu):
#         z = GammaPimu(M,0,1,0)/GammaPimuX(M)/(Zud(FB)**2) # re-scaling of U**2
#     else:
#         z = 1
#     return z

# # def ReBorexinoCC(M):
# #     me = 5.11e-4
# #     if M >= 2*me:
# #         z = sqrt(Gammavee(M, 1, 0, 0)/GammaveeCC(M, 1, 0, 0))
# #     else:
# #         z = 0
# #     return z

# # def ReT2KCC(M):
# #     me = 5.11e-4
# #     if M >= 2*me:
# #         z =  sqrt((GammaPie(M,1,0,0)+Gammavee(M,1,0,0)+Gammavemmup(M,1,0,0)+Gammavmumu(M,1,0,0))/(GammaveeCC(M,1,0,0)+GammaPie(M,1,0,0)+Gammavemmup(M,1,0,0)))
# #     else:
# #         z = 0
# #     return z

# # def ReCHARMCC(M):
# #     me = 5.11e-4
# #     if M >= 2*me:
# #         z = sqrt((Gammavee(M, 1, 0, 0) + Gammavemmup(M, 1, 0, 0) + Gammavmumu(M, 1, 0, 0))/
# #                  (GammaveeCC(M, 1, 0, 0) + Gammavemmup(M, 1, 0, 0)))
# #     else:
# #         z = 0
# #     return z

# # def ReLHCCC(M):
# #     me = 5.11e-4
# #     if M >= 2*me:
# #         z = sqrt((Gammavee(M, 1, 0, 0) + Gammavemmup(M, 1, 0, 0) + Gammavmumu(M, 1, 0, 0))/(GammaveeCC(M, 1, 0, 0) + Gammavemmup(M, 1, 0, 0)))
# #     else:
# #         z = 0
# #     return z


# In[132]:


############ Decay widths for the operators LNQd (label-y) and LNQu (label-z) ############

##### HNL decay into meson #####

def GammaPlYZ(f,mp,mN,ml,mq1,mq2):
    if mN >= (mp+ml):
        z = (f**2)*(mp**4)/(256*mN*(mq1 + mq2)**2)*(mN**2 + ml**2 - mp**2)*sqrt(l(1, xfrac(ml, mN)**2, xfrac(mp, mN)**2))
    else:
        z = 0
    return z

def GammaPieYZ(mN):
    mPi = 0.13957
    fPi = 0.13
    me = 5.11e-4
    mu = 216e-5
    md = 467e-5
    z = GammaPlYZ(fPi, mPi, mN, me, mu, md)
    return z

def GammaPimuYZ(mN):
    mPi = 0.13957
    fPi = 0.13
    mmu = 1.05e-1
    mu = 216e-5
    md = 467e-5
    z = GammaPlYZ(fPi, mPi, mN, mmu, mu, md)
    return z

def GammaKeYZ(mN):
    mK = 0.493
    fK = 0.156
    me = 5.11e-4
    mu = 216e-5
    ms = 9.3e-2
    z = GammaPlYZ(fK, mK, mN, me, mu, ms)
    return z

def GammaKmuYZ(mN):
    mK = 0.493
    fK = 0.156
    mmu = 1.05e-1
    mu = 216e-5
    ms = 9.3e-2
    z = GammaPlYZ(fK, mK, mN, mmu, mu, ms)
    return z

##### Meson decay into HNL #####

def GammaPNlYZ(f,mp,mN,ml,mq1,mq2):
    if mp >= (mN+ml):
        z = (f**2)*(mp**3)/(128*(mq1 + mq2)**2)*(mp**2 - mN**2 - ml**2)*sqrt(l(1, xfrac(ml, mp)**2, xfrac(mN, mp)**2))
    else:
        z = 0
    return z

def GammaPiNeYZ(mN):
    mPi = 0.13957
    fPi = 0.13
    me = 5.11e-4
    mu = 216e-5
    md = 467e-5
    z = GammaPNlYZ(fPi, mPi, mN, me, mu, md)
    return z

def GammaPiNmuYZ(mN):
    mPi = 0.13957
    fPi = 0.13
    mmu = 1.05e-1
    mu = 216e-5
    md = 467e-5
    z = GammaPNlYZ(fPi, mPi, mN, mmu, mu, md)
    return z

def GammaKNeYZ(mN):
    mK = 0.493
    fK = 0.156
    me = 5.11e-4
    mu = 216e-5
    ms = 9.3e-2
    z = GammaPNlYZ(fK, mK, mN, me, mu, ms)
    return z

def GammaKNmuYZ(mN):
    mK = 0.493
    fK = 0.156
    mmu = 1.05e-1
    mu = 216e-5
    ms = 9.3e-2
    z = GammaPNlYZ(fK, mK, mN, mmu, mu, ms)
    return z

# def GammaRhomuYZ(mN):
#     mRho = 0.775
#     fRho = 0.17
#     mmu = 1.05e-1
#     mu = 216e-5
#     md = 467e-5
#     z = GammaPNlYZ(fRho, mRho, mN, mmu, mu, md)
#     z = 0
#     return z

def GammaDNeYZ(mN):
    mD = 1.896
    fD = 0.212
    me = 5.11e-4
    md = 467e-5
    mc = 1.67
    z = GammaPNlYZ(fD, mD, mN, me, mc, md)
    return z

def GammaDsNeYZ(mN):
    mDs = 1.97
    fDs = 0.249
    me = 5.11e-4
    ms = 9.3e-2
    mc = 1.67
    z = GammaPNlYZ(fDs, mDs, mN, me, mc, ms)
    return z

def GammaDNmuYZ(mN):
    mD = 1.896
    fD = 0.212
    mmu = 1.05e-1
    md = 467e-5
    mc = 1.67
    z = GammaPNlYZ(fD, mD, mN, mmu, mc, md)
    return z

def GammaDsNmuYZ(mN):
    mDs = 1.97
    fDs = 0.249
    mmu = 1.05e-1
    ms = 9.3e-2
    mc = 1.67
    z = GammaPNlYZ(fDs, mDs, mN, mmu, mc, ms)
    return z

# def GammaBNeYZ(mN):
#     mB = 5.279
#     fB = 0.19
#     me = 5.11e-4
#     mu = 216e-5
#     mb = 4.68
#     z = GammaPNlYZ(fB, mB, mN, me, mu, mb)
#     return z


# In[141]:


#### electron case ####

# Peak searches. Kinematics of the meson decay products.

def RePSPion4fYZ(M,FB):
    z = (GammaPiNe(M)/GammaPiNeYZ(M))/(Zud(FB)**2) # re-scaling of U**2
    return z

def RePSKaon4fYZ(M,FB):
    z = (GammaKNe(M)/GammaKNeYZ(M))/(Zus(FB)**2) # re-scaling of U**2
    return z

# Decay on flight of the HNL into pion. The HNL come from K leptonic decay.

def ReT2K4fYZ(M,FB):
    mPi = 0.13957
    me = 5.11e-4
    if M > (mPi-me):
        z = sqrt(GammaKNe(M)*(GammaPie(M,1,0,0)+Gammavee(M,1,0,0)+Gammavmumu(M,1,0,0)+Gammavemmup(M,1,0,0))/(GammaKNeYZ(M)*GammaPieYZ(M)))/(Zus(FB)*Zud(FB))
    return z

# Re-analysis of BEBC. D and Ds (they do not say it, but the mass range of the HNL needs a Ds) leptonic decay. HNL decay to leptons, and pion channel. 

def ReBEBC4fYZ(M,FB):
    # mD = 1.896
    mDs = 1.97
    mPi = 0.13957
    me = 5.11e-4
    if M > (mPi+me):
        if M < (mDs-me):
            z = sqrt((GammaDNe(M)+GammaDsNe(M))*(GammaPie(M,1,0,0)+Gammavee(M,1,0,0)+Gammavmumu(M,1,0,0)+Gammavemmup(M,1,0,0))/((GammaDNeYZ(M)*Zcd(FB)+GammaDsNeYZ(M)*Zcs(FB))*GammaPieYZ(M)))/Zud(FB)
    # else:
    #     z = 1
    return z

def ReBEBC4fYZlow(M,FB):
    # mD = 1.896
    mDs = 1.97
    mPi = 0.13957
    me = 5.11e-4
    if M > (mPi+me):
        if M < (mDs-me):
            z = (GammaPie(M,1,0,0)+Gammavee(M,1,0,0)+Gammavmumu(M,1,0,0)+Gammavemmup(M,1,0,0))/GammaPieYZ(M)/(Zud(FB)**2) ## the upper line goes with the width
    # else:
    #     z = 1
    return z

# Belle. B to XlN, with X = D,D*,pi,rho..., nothing and N > pi e. We approximate the width of the semileptonic B decay to that of the standard mixing. Fix this!

# def ReBelle4fYZ(M,FB):
#     mB = 5.279
#     mPi = 0.13957
#     me = 5.11e-4 
#     if M > (mPi+me):
#         if M < (mB-me):
#             #z = sqrt(GammaPie(M,1,0,0)/(GammaPieYZ(M)))/(Zcb(FB)*Zud(FB)) #+ 
#             z = sqrt(GammaBNe(M)*GammaPie(M,1,0,0)/(GammaBNeYZ(M)*GammaPieYZ(M)))/(Zub(FB)*Zud(FB)) # re-scaling of U**2
#     else:
#         z = 1
#     return z

# def ReBelle4fYZlow(M,FB):
#     mB = 5.279
#     mPi = 0.13957
#     me = 5.11e-4
#     if M > (mPi+me):
#         if M < (mB-me):
#             z = GammaPie(M,1,0,0)/GammaPieYZ(M)/(Zud(FB)**2) # re-scaling of U**2
#     else:
#         z = 1
#     return z

#### muon case ####

# Peak searches. Kinematics of the meson decay products.

def RmuPSPion4fYZ(M,FB):
    z = (GammaPiNmu(M)/GammaPiNmuYZ(M))/(Zud(FB)**2) # re-scaling of U**2
    return z

def RmuPSKaon4fYZ(M,FB):
    z = (GammaKNmu(M)/GammaKNmuYZ(M))/(Zus(FB)**2) # re-scaling of U**2
    return z

# Decay on flight of the HNL into pion. The HNL come from K leptonic decay.

def RmuT2K4fYZ(M,FB):
    mPi = 0.13957
    mmu = 1.05e-1
    if M > (mPi+mmu):
        z = sqrt((GammaKNmu(M))*(GammaPimu(M,0,1,0)+Gammavee(M,0,1,0)+Gammavmumu(M,0,1,0)+Gammavmumep(M,0,1,0))/(GammaKNmuYZ(M)*GammaPimuYZ(M)))/(Zus(FB)*Zud(FB))
    else:
        z = 1
    return z

# NuTeV: Kaon, D and Ds to mu HNL, HNL to mu pi, mu rho, and leptonic channels. Only CC channels! Rho width = 0 with this operator!

def RmuNuTeV4fYZ(M,FB):
    mDs = 1.97
    mPi = 0.13957
    mmu = 1.05e-1
    if M > (mPi+mmu):
        if M < (mDs-mmu):
            z = sqrt((GammaKNmu(M)+GammaDNmu(M)+GammaDsNmu(M))*(GammaPimu(M,0,1,0)+GammavmumuCC(M,0,1,0)+Gammavmumep(M,0,1,0))/((GammaKNmuYZ(M)*(Zus(FB)**2)+GammaDNmuYZ(M)*(Zcd(FB)**2)+GammaDsNmuYZ(M)*(Zcs(FB)**2))*(GammaPimuYZ(M))))/Zud(FB)
    else:
        z = 1
    return z

def RmuNuTeV4fYZlow(M,FB):
    mDs = 1.97
    mPi = 0.13957
    mmu = 1.05e-1
    if M > (mPi+mmu):
        if M < (mDs-mmu):
            z = (GammaPimu(M,0,1,0)+GammavmumuCC(M,0,1,0)+Gammavmumep(M,0,1,0))/GammaPimuYZ(M)/(Zud(FB)**2)
    else:
        z = 1
    return z

# BEBC original: D to mu HNL, HNL to pi mu, leptons but CC only!

def RmuBEBC4fYZ(M,FB):
    mPi = 0.13957
    mmu = 1.05e-1 
    if M > (mPi+mmu):
        z = sqrt(GammaDNmu(M)*(GammaPimu(M,0,1,0)+GammavmumuCC(M,0,1,0)+Gammavmumep(M,0,1,0))/(GammaDNmuYZ(M)*(GammaPimuYZ(M))))/(Zcd(FB)*Zud(FB)) # re-scaling of U**2
    else:
        z = 1
    return z

def RmuBEBC4fYZlow(M,FB):
    mPi = 0.13957
    mmu = 1.05e-1 
    if M > (mPi+mmu):
        z = (GammaPimu(M,0,1,0)+GammavmumuCC(M,0,1,0)+Gammavmumep(M,0,1,0))/GammaPimuYZ(M)/(Zud(FB)**2) # re-scaling of U**2
    else:
        z = 1
    return z

# def RmuNA34fYZ(M,FB):
#     mPi = 0.13957
#     mmu = 1.05e-1 
#     mD = 1.896
#     if M > (mPi+mmu):
#         if M >= (mD-mmu):
#             z = sqrt(GammaPimu(M,0,1,0)/(GammaPimuYZ(M)))/(Zcb(FB)*Zud(FB)) # re-scaling of U**2
#         if M < (mD-mmu):
#             z = sqrt(GammaPimu(M,0,1,0)/(GammaPimuYZ(M)))/(Zcs(FB)*Zud(FB)) 
#     else:
#         z = 1
#     return z

# def RmuNA34fYZlow(M,FB):
#     mPi = 0.13957
#     mmu = 1.05e-1
#     if M > (mPi+mmu):
#         z = GammaPimu(M,0,1,0)/GammaPimuYZ(M)/(Zud(FB)**2) # re-scaling of U**2
#     else:
#         z = 1
#     return z

