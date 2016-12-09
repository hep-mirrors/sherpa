#!/usr/bin/env python2

from tensor import tensor
from sym_var import sym_var
from lorentz_structures import mink_metric, gamma_0, Gamma, ProjP, ProjM, gamma_5, four_identity

from sympy import init_printing, simplify, pprint, mathematica_code, Symbol, sqrt, I
from sympy.functions import conjugate as cgt

class tensor1d(tensor):

    def __init__(self, lst, key):
        components = [tensor([item], None) for item in lst]
        super(tensor1d, self).__init__(components, key)

def conjugate_tensor1d(tens):
    assert(len(tens.key_dim_dict())==1)
    arr = [cgt(el._array[0]) for el in tens.elements()]
    return tensor1d(arr, tens._toplevel_key)

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

def test_bar():
    p1 = Symbol('p1', real=True)
    p2 = Symbol('p2', real=True)
    p3 = Symbol('p3', real=True)
    M  = 0 # Symbol('M',  real=True)
    p0 = sqrt(p1**2+p2**2+p3**2+M**2)

    u  = uplus   (p0,p1,p2,p3,'a')
    ub = uplusbar(p0,p1,p2,p3,'a')

    assert(u._array[0]._array[0] - cgt(ub._array[2]._array[0])==0)
    assert(u._array[1]._array[0] - cgt(ub._array[3]._array[0])==0)
    assert(u._array[2]._array[0] - cgt(ub._array[0]._array[0])==0)
    assert(u._array[3]._array[0] - cgt(ub._array[1]._array[0])==0)

    u  = uminus   (p0,p1,p2,p3,'a')
    ub = uminusbar(p0,p1,p2,p3,'a')

    assert(u._array[0]._array[0] - cgt(ub._array[2]._array[0])==0)
    assert(u._array[1]._array[0] - cgt(ub._array[3]._array[0])==0)
    assert(u._array[2]._array[0] - cgt(ub._array[0]._array[0])==0)
    assert(u._array[3]._array[0] - cgt(ub._array[1]._array[0])==0)

def test_spinors(massless=True):

    p1 = Symbol('p1', real=True)
    p2 = Symbol('p2', real=True)
    p3 = Symbol('p3', real=True)
    M  = 0 if massless else Symbol('M',  real=True)
    p0 = sqrt(p1**2+p2**2+p3**2+M**2)
    
    assert((uminus(p0,p1,p2,p3,'a')*uminusbar(p0,p1,p2,p3,'a') )._array[0] == 0)
    assert((vminus(p0,p1,p2,p3,'a')*vminusbar(p0,p1,p2,p3,'a') )._array[0] == 0)
    assert((uminus(p0,p1,p2,p3,'a')*vplusbar (p0,p1,p2,p3,'a') )._array[0] == 0)
    assert((uplus (p0,p1,p2,p3,'a')*vminusbar(p0,p1,p2,p3,'a') )._array[0] == 0)

    dop = dirac_op([p0,p1,p2,p3], M, 'a', 'b')

    d1 = dop *(uplus    (p0,p1,p2,p3,'b'))
    d2 = dop *(uminus   (p0,p1,p2,p3,'b'))
    d3 = dop *(vplus    (p0,p1,p2,p3,'b'))
    d4 = dop *(vminus   (p0,p1,p2,p3,'b'))

    for i in range(4):
        assert(d1._array[i]._array[0].simplify() == 0 )
        assert(d2._array[i]._array[0].simplify() == 0 )
        assert(d3._array[i]._array[0].simplify() == 0 )
        assert(d4._array[i]._array[0].simplify() == 0 )

test_bar()
test_spinors()


