#!/bin/env python2 

def require(test, message):
    if not test:
        raise RuntimeError("Test \"{0}\" failed".format(message))
    print "Test \"{0}\" passed".format(message)

from tensor import tensor, multiply, lorentz_key

metric = tensor(
    [tensor([tensor([-1], None), tensor([0], None), tensor([0], None), tensor([0], None)], 'mu_2'),
     tensor([tensor([ 0], None), tensor([1], None), tensor([0], None), tensor([0], None)], 'mu_2'), 
     tensor([tensor([ 0], None), tensor([0], None), tensor([1], None), tensor([0], None)], 'mu_2'),
     tensor([tensor([ 0], None), tensor([0], None), tensor([0], None), tensor([1], None)], 'mu_2')],'mu_1') 

metric_prime = tensor(
    [tensor([tensor([-1], None), tensor([0], None), tensor([0], None), tensor([0], None)], 'mu_1'),
     tensor([tensor([ 0], None), tensor([1], None), tensor([0], None), tensor([0], None)], 'mu_1'), 
     tensor([tensor([ 0], None), tensor([0], None), tensor([1], None), tensor([0], None)], 'mu_1'),
     tensor([tensor([ 0], None), tensor([0], None), tensor([0], None), tensor([1], None)], 'mu_1')],'mu_2') 

identity_4 = tensor(
    [tensor([tensor([1], None), tensor([0], None), tensor([0], None), tensor([0], None)], 'mu_2'),
     tensor([tensor([0], None), tensor([1], None), tensor([0], None), tensor([0], None)], 'mu_2'), 
     tensor([tensor([0], None), tensor([0], None), tensor([1], None), tensor([0], None)], 'mu_2'),
     tensor([tensor([0], None), tensor([0], None), tensor([0], None), tensor([1], None)], 'mu_2')],'mu_1') 

vector_a = tensor([tensor([12.345], None), tensor([3], None), tensor([-12], None), tensor([1.2], None)], 'mu_1')
vector_b = tensor([tensor([-12.345], None), tensor([3], None), tensor([-12], None), tensor([1.2], None)], 'mu_2')
vector_c = tensor([tensor([12.345], None), tensor([3], None), tensor([-12], None), tensor([1.2], None)], 'mu_3')
vector_d = tensor([tensor([12.345], None), tensor([3], None), tensor([-12], None), tensor([1.2], None), tensor([1], None)], 'mu_1')

require(vector_a!=vector_b, "Four-vector inequality (value)")
require(vector_a!=vector_c, "Four-vector inequality (index key)")
require(vector_a!=vector_d, "Four-vector inequality (toplevel dim)")
require(metric!=identity_4, "4x4 matrix inequality (values)")
require(metric==metric_prime, "4x4 matrix key ordering")
require( all([(metric[{'mu_1':i}]._array == metric[{'mu_2':i}]._array) for i in range(4)]) , "Lorentz metric symmetry")
require(metric*metric == tensor([4], None), "Lorentz metric trace")
require((vector_a*metric==vector_b), "Four-vector/metric multiplication")
require((vector_b*metric==vector_a), "Four-vector/metric  multiplication")
require((vector_b*metric==metric*vector_b), "Four-vector/metric multiplication commutativity")

vector = tensor([tensor([12.345], None), tensor([3], None), tensor([-12], None), tensor([1.2], None)], 'mu_5')
vector[{'mu_5':0}] = tensor([0], None)
vector2 = tensor([tensor([5], None), tensor([5], None), tensor([5], None), tensor([5], None)], 'mu_2')
metric[{'mu_1':1}] = vector2
vector3 = tensor([tensor([5], None), tensor([5], None), tensor([5], None), tensor([5], None)], 'mu_1')
metric[{'mu_2':1}] = vector3
summetric = metric+metric
metric = tensor(
    [tensor([tensor([-1], None), tensor([0], None), tensor([0], None), tensor([0], None)], 'mu_2'),
     tensor([tensor([0],  None), tensor([1], None), tensor([0], None), tensor([0], None)], 'mu_2'), 
     tensor([ tensor([0], None), tensor([0], None), tensor([1], None), tensor([0], None)], 'mu_2'),
     tensor([ tensor([0], None), tensor([0], None), tensor([0], None), tensor([1], None)], 'mu_2')],'mu_1') 
vector = tensor([tensor([5], None), tensor([5], None), tensor([5], None), tensor([5], None)], 'mu_1')


from lorentz_structures import *

# check gamma_5 = i gamma_0*gamma_1*gamma_2*gamma_3
test  = gamma_0(1,'a')*gamma_1('a','b')*gamma_2('b','c')*gamma_3('c',2)*tensor([complex(0,1)], None)
test_ = gamma_5(1,2)
require ( test_ == test, "Fifth gamma matrix consistency")

# square of the charge conjugation matrix 
# in our convention is -identity
require(C('mu_1',3)*C(3,'mu_2')*tensor([-1], None) == identity_4, "Charge conjugation square")

# test anticommutator of dirac matrices, test entries on very low level (elementary tensors)
anticomm = Gamma('mu','i','k') *  Gamma('nu','k','j') + Gamma('nu','i','k') *  Gamma('mu','k','j')
test = 1
for ind in range(1,4):
    for i in range (4):
        test *= int(anticomm[{'mu':ind, 'nu':ind, 'i':i, 'j':i }]._array[0] == -2)
for i in range (4):
    test *= int(anticomm[{'mu':0, 'nu':0, 'i':i, 'j':i }]._array[0] == 2)
require(test, "Dirac matrix anticommutator")

vector  = tensor([tensor([2], None), tensor([0], None), tensor([0], None), tensor([2], None)], lorentz_key('mu_5'))
vector2 = tensor([tensor([2], None), tensor([0], None), tensor([0], None), tensor([2], None)], lorentz_key('mu_5'))
require( ((vector*vector2)._array[0] == 0), "Mass")

vector  = tensor([tensor([2], None), tensor([0], None), tensor([0], None), tensor([2], None)], lorentz_key('mu_5'))
vector2 = tensor([tensor([2], None), tensor([0], None), tensor([2], None), tensor([0], None)], lorentz_key('mu_5'))
require(((vector*vector2)._array[0] == 4), "Mass")

vector  = tensor([tensor([0], None), tensor([0], None), tensor([0], None), tensor([2], None)], lorentz_key('mu_5'))
vector2 = tensor([tensor([2], None), tensor([0], None), tensor([0], None), tensor([1], None)], lorentz_key('mu_5'))
require(((vector*vector2)._array[0] == -2), "Mass")

even = [0,1,2,3]
require((Epsilon(1,2,3,4))[{1:even[0], 2:even[1], 3:even[2], 4:even[3]}]._array[0] == 1, "Epsilon 1")

even = [1,0,2,3]
require((Epsilon(1,2,3,4))[{1:even[0], 2:even[1], 3:even[2], 4:even[3]}]._array[0] == -1, "Epsilon 2")

even = [1,0,3,2]
require((Epsilon(1,2,3,4))[{1:even[0], 2:even[1], 3:even[2], 4:even[3]}]._array[0] == 1, "Epsilon 3")

vector  = tensor([tensor([2], None), tensor([3.21], None), tensor([12.2], None), tensor([2], None)], lorentz_key('mu_5'))
vec2 = 2.0*2.0 - 3.21*3.21 - 12.2*12.2 - 2.0*2.0
require(
    (vector**2 == vector*vector) and
    (vector**2 == tensor([vec2], None) and
     (vector**4 == tensor([vec2*vec2], None))), "Tensor power operator")

