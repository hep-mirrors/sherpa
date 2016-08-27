#!/bin/env python2 

def require(test, message):
    if not test:
        raise RuntimeError("Test \"{0}\" failed".format(message))
    print "Test \"{0}\" passed".format(message)

from math import sqrt
from color_structures import *

def fuzz_eq(a,b):
    if isinstance(a,complex):
        return abs(a.real - b)<1e-14 and b.imag<1e-14
    else:
        return abs(a-b)<1e-14

def test_f():

    a = f('a','b','c')

    val = complex(a[{color_key('a', 'ad'):0,color_key('b', 'ad'):1,color_key('c', 'ad'):2}]._array[0])
    assert(fuzz_eq(val, 1.0))

    val = complex(a[{color_key('a', 'ad'):0,color_key('b', 'ad'):3,color_key('c', 'ad'):6}]._array[0])
    assert(fuzz_eq(val, 0.5))

    val = complex(a[{color_key('a', 'ad'):0,color_key('b', 'ad'):3,color_key('c', 'ad'):6}]._array[0])
    assert(fuzz_eq(val, 0.5))

    val = complex(a[{color_key('a', 'ad'):0,color_key('b', 'ad'):4,color_key('c', 'ad'):5}]._array[0])
    assert(fuzz_eq(val,-0.5))

    val = complex(a[{color_key('a', 'ad'):1,color_key('b', 'ad'):3,color_key('c', 'ad'):5}]._array[0])
    assert(fuzz_eq(val, 0.5))

    val = complex(a[{color_key('a', 'ad'):1,color_key('b', 'ad'):4,color_key('c', 'ad'):6}]._array[0])
    assert(fuzz_eq(val, 0.5))

    val = complex(a[{color_key('a', 'ad'):2,color_key('b', 'ad'):3,color_key('c', 'ad'):4}]._array[0])
    assert(fuzz_eq(val, 0.5))

    val = complex(a[{color_key('a', 'ad'):2,color_key('b', 'ad'):5,color_key('c', 'ad'):6}]._array[0])
    assert(fuzz_eq(val,-0.5))

    val = complex(a[{color_key('a', 'ad'):3,color_key('b', 'ad'):4,color_key('c', 'ad'):7}]._array[0])
    assert(fuzz_eq(val,sqrt(3.)/2.0))

    val = complex(a[{color_key('a', 'ad'):5,color_key('b', 'ad'):6,color_key('c', 'ad'):7}]._array[0])
    assert(fuzz_eq(val,sqrt(3.)/2.0))

def test_trace():
    for a in range(8):
        for b in range(8):
            # check trace( Ta*Tb ) = delta_{ab}
            p =  T_a(a,'i','j')*T_a(b,'j','i')
            p = complex(p._array[0])
            if a!=b:
                assert(p==0)
            else:
                assert(abs(p.real-pfac) < 1e-14)
                assert(abs(p.imag     ) < 1e-14)


if __name__ == "__main__":

    

    #
    
    #from time import time
    #t0 = time()
    #print [el for el in (f(-2,1,2)).elements()]

    #f(-2,1,2)*f(-1,-2,3)*f(4,5,-1)
    #print time() - t0

    #     #f(1,2,3)*f(2,4,5)                # 14.2726690769

    #     #f(1,2,3)*f(6,4,5)                # 15.1672339439

    #     #f(1,2,3)*f(1,2,3)                # 5.99734210968
        
    #     #f(1,2,3) # 6.03008389473
    #     #f(2,4,5) #
        
    #     a=f(-2,1,2)
    #     b=f(-1,-2,3)
    #     c=f(4,5,-1)
    #     print "Constructors ", time() - t0
    #     d = a*b
    #     print "a*b=d", time() - t0
    #     e = d*c
    #     print "e=d*c ", time() - t0

    test_trace()
    test_f()

