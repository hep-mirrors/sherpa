from tensor import tensor, new, lorentz_key
from sympy.parsing.sympy_parser import parse_expr
from itertools import permutations

#################################
# Dirac matrices as defined in  #
# master's thesis (A.3) - (A.5) #
#################################
 
I = complex(0,1)

def four_identity(i,j):

    return tensor([tensor([tensor([ 1.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None)], j),
                   tensor([tensor([ 0.0 ], None), tensor([ 1.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None)], j),
                   tensor([tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ 1.0 ], None), tensor([ 0.0 ], None)], j),
                   tensor([tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ 1.0 ], None)], j)], i)

def mink_metric(i,j):

    return tensor([tensor([tensor([ 1.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None)], j),
                   tensor([tensor([ 0.0 ], None), tensor([-1.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None)], j),
                   tensor([tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([-1.0 ], None), tensor([ 0.0 ], None)], j),
                   tensor([tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([-1.0 ], None)], j)], i)
    
def gamma_0(i,j):

    return tensor([tensor([tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ 1.0 ], None), tensor([ 0.0 ], None)], j),
                   tensor([tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ 1.0 ], None)], j),
                   tensor([tensor([ 1.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None)], j),
                   tensor([tensor([ 0.0 ], None), tensor([ 1.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None)], j)], i)

def gamma_1(i,j):

    return tensor([tensor([tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ 1.0 ], None)], j),
                   tensor([tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ 1.0 ], None), tensor([ 0.0 ], None)], j),
                   tensor([tensor([ 0.0 ], None), tensor([-1.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None)], j),
                   tensor([tensor([-1.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None)], j)], i)

def gamma_2(i,j):

    return tensor([tensor([tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ -I  ], None)], j),
                   tensor([tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([  I  ], None), tensor([ 0.0 ], None)], j),
                   tensor([tensor([ 0.0 ], None), tensor([  I  ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None)], j),
                   tensor([tensor([ -I  ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None)], j)], i)

def gamma_3(i,j):

    return tensor([tensor([tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ 1.0 ], None), tensor([ 0.0 ], None)], j),
                   tensor([tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([-1.0 ], None)], j),
                   tensor([tensor([-1.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None)], j),
                   tensor([tensor([ 0.0 ], None), tensor([ 1.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None)], j)], i)

def gamma_5(i,j):
    return tensor([tensor([tensor([-1.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None)], j),
                   tensor([tensor([ 0.0 ], None), tensor([-1.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None)], j),
                   tensor([tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ 1.0 ], None), tensor([ 0.0 ], None)], j),
                   tensor([tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ 0.0 ], None), tensor([ 1.0 ], None)], j)], i)


#######################################
# elementary lorentz structures from  # 
# table 6, arXiv:1108.2040v2 [hep-ph] #
#######################################

class C(tensor):
    
    def __init__(self,i,j):
        array = [tensor([tensor([ 0.0], None), tensor([ 1.0], None), tensor([ 0.0], None), tensor([ 0.0], None)], j),
                 tensor([tensor([-1.0], None), tensor([ 0.0], None), tensor([ 0.0], None), tensor([ 0.0], None)], j),
                 tensor([tensor([ 0.0], None), tensor([ 0.0], None), tensor([ 0.0], None), tensor([-1.0], None)], j),
                 tensor([tensor([ 0.0], None), tensor([ 0.0], None), tensor([ 1.0], None), tensor([ 0.0], None)], j)]

        super(C,self).__init__(array, i)


class Gamma(tensor):
    
    def __init__(self, mu, i, j):
        mul = lorentz_key(mu) if (mu<0) else mu
        array = [gamma_0(i,j),
                 gamma_1(i,j),
                 gamma_2(i,j),
                 gamma_3(i,j),]
        
        super(Gamma,self).__init__(array, mul)

class Gamma5(tensor):
    
    def __init__(self,i, j):
        tmp = gamma_5(i,j)
        super(Gamma5,self).__init__(tmp._array, tmp._toplevel_key)


class Metric(tensor):
    
    def __init__(self,i,j):
        il = lorentz_key(i) if (i<0) else i
        jl = lorentz_key(j) if (j<0) else j
        array = [tensor([tensor([ 1.0], None), tensor([ 0.0], None), tensor([ 0.0], None), tensor([ 0.0], None)], jl),
                 tensor([tensor([ 0.0], None), tensor([-1.0], None), tensor([ 0.0], None), tensor([ 0.0], None)], jl),
                 tensor([tensor([ 0.0], None), tensor([ 0.0], None), tensor([-1.0], None), tensor([ 0.0], None)], jl),
                 tensor([tensor([ 0.0], None), tensor([ 0.0], None), tensor([ 0.0], None), tensor([-1.0], None)], jl)]

        super(Metric,self).__init__(array,il)

class P(tensor):
    
    def __init__(self,i,n):
        # our inde scheme starts at 0, so subtract 1 from n
        # when naming the momenta
        k = n-1
        il = lorentz_key(i) if (i<0) else i
        p_str = "p{0}{1}"
        array = [tensor([parse_expr(p_str.format(k,0))], None), 
                 tensor([parse_expr(p_str.format(k,1))], None), 
                 tensor([parse_expr(p_str.format(k,2))], None), 
                 tensor([parse_expr(p_str.format(k,3))], None)]

        super(P,self).__init__(array,il)


class ProjM(tensor):

    def __init__(self,i,j):
        temp = (four_identity(i,j) - gamma_5(i,j))*(0.5)
        
        super(ProjM,self).__init__(temp._array,temp._toplevel_key)


class ProjP(tensor):

    def __init__(self,i,j):
        temp = (four_identity(i,j) + gamma_5(i,j))*(0.5)
        
        super(ProjP,self).__init__(temp._array,temp._toplevel_key)
        

class Identity(tensor):

    def __init__(self,i,j):
        array = [tensor([tensor([1.0], None), tensor([0.0], None), tensor([0.0], None), tensor([0.0], None)], j),
                 tensor([tensor([0.0], None), tensor([1.0], None), tensor([0.0], None), tensor([0.0], None)], j),
                 tensor([tensor([0.0], None), tensor([0.0], None), tensor([1.0], None), tensor([0.0], None)], j),
                 tensor([tensor([0.0], None), tensor([0.0], None), tensor([0.0], None), tensor([1.0], None)], j)]

        super(Identity,self).__init__(array,i)

class Epsilon(tensor):

    @staticmethod
    def parity(lst):
        parity = 1
        for i in range(0,len(lst)-1):
            if lst[i] != i:
                parity *= -1
                mn = min(range(i,len(lst)), key=lst.__getitem__)
                lst[i],lst[mn] = lst[mn],lst[i]
        return parity    

    def __init__(self,i,j,k,l):

        perms = list(permutations([0,1,2,3]))
        
        il = lorentz_key(i) if (i<0) else i
        jl = lorentz_key(j) if (j<0) else j
        kl = lorentz_key(k) if (k<0) else k
        ll = lorentz_key(l) if (l<0) else l
        tmp = new({il:4,jl:4,kl:4,ll:4})

        for perm in perms:
            tmp.__setitem__({il:perm[0], jl:perm[1], kl:perm[2], ll:perm[3]}, 
                            tensor([self.parity(list(perm))], None))
        
        super(Epsilon,self).__init__(tmp._array, tmp._toplevel_key)
        

##########################
# Some helper functions  # 
##########################

# Check if tensor is usual
# gamma_mu structure. Need
# to stupidly check all per-
# mutations of key assignments
def is_ffv(tns):
    if not isinstance(tns, tensor): return False
    keys = tns.key_dim_dict().keys()
    if len(keys)!=3: return False
    dims = tns.key_dim_dict().values()
    if not all([dim==4 for dim in dims]): return False
    if (tns==Gamma(keys[0],keys[1],keys[2])
        or tns==Gamma(keys[0],keys[2],keys[1])
        or tns==Gamma(keys[1],keys[0],keys[2])
        or tns==Gamma(keys[1],keys[2],keys[0])
        or tns==Gamma(keys[2],keys[1],keys[0])
        or tns==Gamma(keys[2],keys[0],keys[1])):
        return True
    return False
    
    
# Check if tensor is usual
# triple-gluon-like structure
def is_vvv(tns):
    if not isinstance(tns, tensor): return False
    keys = tns.key_dim_dict().keys()
    if len(keys)!=3: return False
    dims = tns.key_dim_dict().values()
    if not all([dim==4 for dim in dims]): return False
    vvv = (P(3,1)*Metric(1,2) - P(3,2)*Metric(1,2) - P(2,1)*Metric(1,3) + P(2,3)*Metric(1,3) + P(1,2)*Metric(2,3) - P(1,3)*Metric(2,3))
    if tns==vvv:
        return True
    return False
