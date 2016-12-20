#!/usr/bin/env python2

from tensor import tensor
from sym_var import sym_var
from lorentz_structures import mink_metric, gamma_0, Gamma, ProjP, ProjM, gamma_5, four_identity

from sympy import init_printing, simplify, pprint, mathematica_code, Symbol, sqrt, I
from sympy.functions import conjugate as cgt
from sympy.functions import Abs
from sympy import N as evaluate

from random import uniform

class tensor1d(tensor):

    def __init__(self, lst, key):
        components = [tensor([item], None) for item in lst]
        super(tensor1d, self).__init__(components, key)

def conjugate_tensor1d(tens):
    assert(len(tens.key_dim_dict())==1)
    arr = [cgt(el._array[0]) for el in tens.elements()]
    return tensor1d(arr, tens._toplevel_key)

# Normalization: uminus*uminusbar = 2M
# Normalization: uplus *uplusbar  = 2M

def uplus(p0,p1,p2,p3, key):
    pbar  = sqrt(p1**2+p2**2+p3**2)
    ppl = pbar + p3 
    pmi = pbar - p3
    ptr = p1 + I*p2

    return 1/sqrt(2*pbar*ppl)*tensor1d([sqrt(p0-pbar)*ppl,
                                        sqrt(p0-pbar)*ptr,
                                        sqrt(p0+pbar)*ppl,
                                        sqrt(p0+pbar)*ptr], key)

def vminus(p0,p1,p2,p3, key):
    pbar  = sqrt(p1**2+p2**2+p3**2)
    ppl = pbar + p3 
    pmi = pbar - p3
    ptr = p1 + I*p2

    return 1/sqrt(2*pbar*ppl)*tensor1d([-sqrt(p0-pbar)*ppl,
                                        -sqrt(p0-pbar)*ptr,
                                        +sqrt(p0+pbar)*ppl,
                                        +sqrt(p0+pbar)*ptr], key)
                                       
def uminus(p0,p1,p2,p3, key):
    pbar  = sqrt(p1**2+p2**2+p3**2)
    ppl = pbar + p3 
    pmi = pbar - p3
    ptr = p1 + I*p2

    return -1/sqrt(2*pbar*ppl)*tensor1d([sqrt(p0+pbar)*(-cgt(ptr)),
                                         sqrt(p0+pbar)*ppl,
                                         sqrt(p0-pbar)*(-cgt(ptr)),
                                         sqrt(p0-pbar)*ppl], key)

def vplus(p0,p1,p2,p3, key):
    pbar  = sqrt(p1**2+p2**2+p3**2)
    ppl = pbar + p3 
    pmi = pbar - p3
    ptr = p1 + I*p2

    return -1/sqrt(2*pbar*ppl)*tensor1d([+sqrt(p0+pbar)*(-cgt(ptr)),
                                         +sqrt(p0+pbar)*ppl,
                                         -sqrt(p0-pbar)*(-cgt(ptr)),
                                         -sqrt(p0-pbar)*ppl], key)

def uplusbar(k0,k1,k2,k3, key):
    return conjugate_tensor1d(uplus (k0,k1,k2,k3,'dummy'))*gamma_0('dummy', key)

def uminusbar(k0,k1,k2,k3, key):
    return conjugate_tensor1d(uminus(k0,k1,k2,k3,'dummy'))*gamma_0('dummy', key)

def vplusbar(k0,k1,k2,k3, key):
    return conjugate_tensor1d(vplus (k0,k1,k2,k3,'dummy'))*gamma_0('dummy', key)

def vminusbar(k0,k1,k2,k3, key):
    return conjugate_tensor1d(vminus(k0,k1,k2,k3,'dummy'))*gamma_0('dummy', key)

def dirac_op(p,
             M,
             key_a, key_b):
    [p0,p1,p2,p3] = p
    return tensor1d([+p0,p1,p2,p3], 'mu')*mink_metric('mu','nu')*Gamma('nu', key_a, key_b) - four_identity(key_a, key_b)*M

def epsilonplus(p0,p1,p2,p3, key, k0,k1,k2,k3):
    return 1/sqrt(2)*(uminusbar(k0,k1,k2,k3, 'a')*Gamma(key,'a','b')*uminus(p0,p1,p2,p2,'b'))/(uminusbar(k0,k1,k2,k3, 'c')*uplus(p0,p1,p2,p2,'c'))

def epsilonminus(p0,p1,p2,p3, key, k0,k1,k2,k3):
    return -1/sqrt(2)*(uplusbar(k0,k1,k2,k3, 'a')*Gamma(key,'a','b')*uplus(p0,p1,p2,p2,'b'))/(uplusbar(k0,k1,k2,k3, 'c')*uminus(p0,p1,p2,p2,'c'))

def test_spinors():

    p1 = Symbol('p1', real=True)
    p2 = Symbol('p2', real=True)
    p3 = Symbol('p3', real=True)
    M  = Symbol('M',  real=True)
    p0 = sqrt(p1**2+p2**2+p3**2+M**2)

    dct = {var:uniform(1.,100.) for var in [p1,p2,p3,M]}
    
    # Check normalization of spinors: u*ubar = +2M, v*vbar = -2M
    assert(evaluate((uminus(p0,p1,p2,p3,'a')*uminusbar(p0,p1,p2,p3,'a') )._array[0].subs(dct)) - 2.0*M.subs(dct) < 1.e-11)
    assert(evaluate((vminus(p0,p1,p2,p3,'a')*vminusbar(p0,p1,p2,p3,'a') )._array[0].subs(dct)) + 2.0*M.subs(dct) < 1.e-11)
    assert(evaluate((uplus (p0,p1,p2,p3,'a')*uplusbar (p0,p1,p2,p3,'a') )._array[0].subs(dct)) - 2.0*M.subs(dct) < 1.e-11)
    assert(evaluate((vplus (p0,p1,p2,p3,'a')*vplusbar (p0,p1,p2,p3,'a') )._array[0].subs(dct)) + 2.0*M.subs(dct) < 1.e-11)

    dop_u = dirac_op([p0,p1,p2,p3], +M, 'a', 'b')
    dop_v = dirac_op([p0,p1,p2,p3], -M, 'a', 'b')

    d1 = dop_u *(uplus    (p0,p1,p2,p3,'b'))
    d2 = dop_u *(uminus   (p0,p1,p2,p3,'b'))
    d3 = dop_v *(vplus    (p0,p1,p2,p3,'b'))
    d4 = dop_v *(vminus   (p0,p1,p2,p3,'b'))

    for i in range(4):
        assert(Abs(evaluate(d1._array[i]._array[0].subs(dct))) < 1.e-11)
        assert(Abs(evaluate(d2._array[i]._array[0].subs(dct))) < 1.e-11)
        assert(Abs(evaluate(d3._array[i]._array[0].subs(dct))) < 1.e-11)
        assert(Abs(evaluate(d4._array[i]._array[0].subs(dct))) < 1.e-11)

def test_polvecs():
    # Outgoing gluon mom
    k1 = Symbol('k1', real=True)
    k2 = Symbol('k2', real=True)
    k3 = Symbol('k3', real=True)
    k0 = sqrt(k1**2+k2**2+k3**2)

    # Outgoing gluon reference mom
    g1 = Symbol('g1', real=True)
    g2 = Symbol('g2', real=True)
    g3 = Symbol('g3', real=True)
    g0 = sqrt(k1**2+k2**2+k3**2)

    t = epsilonplus(k0,k1,k2,k3, 'mu', g0, g1, g2, g3)*conjugate_tensor1d(epsilonplus(k0,k1,k2,k3, 'nu', g0, g1, g2, g3)) + epsilonminus(k0,k1,k2,k3, 'mu', g0, g1, g2, g3)*conjugate_tensor1d(epsilonminus(k0,k1,k2,k3, 'nu', g0, g1, g2, g3))

    s = -mink_metric('mu','nu') + (tensor1d([k0,k1,k2,k3],'mu')*tensor1d([g0,g1,g2,g3],'nu')+tensor1d([k0,k1,k2,k3],'nu')*tensor1d([g0,g1,g2,g3],'mu'))/(tensor1d([k0,k1,k2,k3],'nu')*mink_metric('mu','nu')*tensor1d([g0,g1,g2,g3],'mu'))

    dct = {k:uniform(1.,100.) for k in [k1,k2,k3,g1,g2,g3]}
    print evaluate(t._array[0]._array[0]._array[0].subs(dct))
    print evaluate(s._array[0]._array[0]._array[0].subs(dct))


def t0():
    # Incoming quark mom
    t1 = Symbol('t1', real=True)
    t2 = Symbol('t2', real=True)
    t3 = Symbol('t3', real=True)
    t0 = sqrt(t1**2+t2**2+t3**2)
    # Outgoing quark mom
    p1 = Symbol('p1', real=True)
    p2 = Symbol('p2', real=True)
    p3 = Symbol('p3', real=True)
    p0 = sqrt(p1**2+p2**2+p3**2)
    
    # Outgoing gluon mom
    k1 = Symbol('k1', real=True)
    k2 = Symbol('k2', real=True)
    k3 = Symbol('k3', real=True)
    k0 = sqrt(k1**2+k2**2+k3**2)

    # Outgoing gluon reference mom
    g1 = Symbol('g1', real=True)
    g2 = Symbol('g2', real=True)
    g3 = Symbol('g3', real=True)
    g0 = sqrt(k1**2+k2**2+k3**2)


    m2 = uplusbar(p0,p1,p2,p3,'a')*Gamma('mu','a','b')*uplus(k0,k1,k2,k3,'b')*mink_metric('mu','nu')*conjugate_tensor1d(epsilonplus(k0,k1,k2,k3, 'nu', g0,g1,g2,g3))

    print m2._array[0].simplify()

    
test_spinors()
#t0()

#test_polvecs()

