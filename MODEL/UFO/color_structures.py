from math import sqrt
from tensor import tensor, multiply, new, color_key

# Generators of SU(3), in following normalisation
# pfac = 1./sqrt(2.0) : (T^a*T^b) = 1.0 * \delta_{ab}
# pfac = 1./(2.0)     : (T^a*T^b) = 1/2 * \delta_{ab}
pfac = 1./2.0

def T1(i,j):
    return pfac*tensor([tensor([tensor([ 0.0], None), tensor([ 1.0], None), tensor([ 0.0], None)], j),
                        tensor([tensor([ 1.0], None), tensor([ 0.0], None), tensor([ 0.0], None)], j),
                        tensor([tensor([ 0.0], None), tensor([ 0.0], None), tensor([ 0.0], None)], j)], i)

def T2(i,j):
    I = complex(0,1)
    return pfac*tensor([tensor([tensor([ 0.0], None), tensor([ -I ], None), tensor([ 0.0], None)], j),
                        tensor([tensor([  I ], None), tensor([ 0.0], None), tensor([ 0.0], None)], j),
                        tensor([tensor([ 0.0], None), tensor([ 0.0], None), tensor([ 0.0], None)], j)], i)
    
def T3(i,j):
    return pfac*tensor([tensor([tensor([ 1.0], None), tensor([ 0.0], None), tensor([ 0.0], None)], j),
                        tensor([tensor([ 0.0], None), tensor([-1.0], None), tensor([ 0.0], None)], j),
                        tensor([tensor([ 0.0], None), tensor([ 0.0], None), tensor([ 0.0], None)], j)], i)

def T4(i,j):
    return pfac*tensor([tensor([tensor([ 0.0], None), tensor([ 0.0], None), tensor([ 1.0], None)], j),
                        tensor([tensor([ 0.0], None), tensor([ 0.0], None), tensor([ 0.0], None)], j),
                        tensor([tensor([ 1.0], None), tensor([ 0.0], None), tensor([ 0.0], None)], j)], i)

def T5(i,j):
    I = complex(0,1)
    return pfac*tensor([tensor([tensor([ 0.0], None), tensor([ 0.0], None), tensor([ -I ], None)], j),
                        tensor([tensor([ 0.0], None), tensor([ 0.0], None), tensor([ 0.0], None)], j),
                        tensor([tensor([  I ], None), tensor([ 0.0], None), tensor([ 0.0], None)], j)], i)

def T6(i,j):
    return pfac*tensor([tensor([tensor([ 0.0], None), tensor([ 0.0], None), tensor([ 0.0], None)], j),
                        tensor([tensor([ 0.0], None), tensor([ 0.0], None), tensor([ 1.0], None)], j),
                        tensor([tensor([ 0.0], None), tensor([ 1.0], None), tensor([ 0.0], None)], j)], i)

def T7(i,j):
    I = complex(0,1)
    return pfac*tensor([tensor([tensor([ 0.0], None), tensor([ 0.0], None), tensor([ 0.0], None)], j),
                        tensor([tensor([ 0.0], None), tensor([ 0.0], None), tensor([ -I ], None)], j),
                        tensor([tensor([ 0.0], None), tensor([  I ], None), tensor([ 0.0], None)], j)], i)

def T8(i,j):
    return pfac/sqrt(3)*tensor([tensor([tensor([ 1.0], None), tensor([ 0.0], None), tensor([ 0.0], None)], j),
                                tensor([tensor([ 0.0], None), tensor([ 1.0], None), tensor([ 0.0], None)], j),
                                tensor([tensor([ 0.0], None), tensor([ 0.0], None), tensor([-2.0], None)], j)], i)

def T_a(a, i, j):
    if a  == 0: return T1(i,j)
    if a  == 1: return T2(i,j)
    if a  == 2: return T3(i,j)
    if a  == 3: return T4(i,j)
    if a  == 4: return T5(i,j)
    if a  == 5: return T6(i,j)
    if a  == 6: return T7(i,j)
    if a  == 7: return T8(i,j)
    raise ValueError('Invalid argument: a={0}'.format(a))

def ef(a,b,c):
    """Return numerical value of f^{abc} with a,b,c integers"""
    # Calculate commutator, [T^a,T^b] = if^{abc}T^c
    comm = (T_a(a,'i','k')*T_a(b,'k','j') - T_a(b,'i','k')*T_a(a,'k','j'))
    # Project onto basis element T^c and return coefficient,
    # by using Tr(T^a*T^b) = 0.5 * \delta_{ab}
    return -complex(0,1.0)/(pfac)*comm*T_a(c,'j','i')

def replacer_T(a,i,j):
    return sqrt(2.0)*tensor([T1(i,j), T2(i,j), T3(i,j), T4(i,j), 
                             T5(i,j), T6(i,j), T7(i,j), T8(i,j)], a)

#######################################
# elementary color structures from    # 
# table 5, arXiv:1108.2040v2 [hep-ph] #
#######################################

def f(a,b,c):
    """Commutator structure constants"""
    ret = new(color_key(a, 'ad'), 8, {color_key(b, 'ad'):8,color_key(c, 'ad'):8})
    for i in range(8):
        for j in range(8):
            for k in range(8):
                if(i==j or i==k or j==k): continue
                ret[{color_key(a, 'ad'):i,
                     color_key(b, 'ad'):j,
                     color_key(c, 'ad'):k}] = ef(i,j,k)
    return ret

def T(a,i,j):
    """Fundamental representation matrices"""
    return tensor([T1(color_key(i, 'fu'),color_key(j, 'af')),
                   T2(color_key(i, 'fu'),color_key(j, 'af')),
                   T3(color_key(i, 'fu'),color_key(j, 'af')),
                   T4(color_key(i, 'fu'),color_key(j, 'af')),
                   T5(color_key(i, 'fu'),color_key(j, 'af')),
                   T6(color_key(i, 'fu'),color_key(j, 'af')),
                   T7(color_key(i, 'fu'),color_key(j, 'af')),
                   T8(color_key(i, 'fu'),color_key(j, 'af'))], color_key(a, 'ad'))

def Identity(i,j):
    """Kronecker deltas of fundamental representation"""
    return tensor([tensor([tensor([1.0], None), tensor([0.0], None), tensor([0.0], None)], color_key(i, 'fu')),
                   tensor([tensor([0.0], None), tensor([1.0], None), tensor([0.0], None)], color_key(i, 'fu')),
                   tensor([tensor([0.0], None), tensor([0.0], None), tensor([1.0], None)], color_key(i, 'fu'))], color_key(j, 'af'))

def IdentityG(i,j):
    """Kronecker deltas of adjoint representation"""
    ret = new(color_key(i, 'ad'), 8, {color_key(j, 'ad'):8})
    for c in range(8):
        ret[{color_key(i, 'ad'):c, 
             color_key(j, 'ad'):c}]._array[0] = 1.0
    return ret
