from tensor import tensor, new, lorentz_key
from c_variable import c_variable
from itertools import permutations

#################################
# Dirac matrices as defined in  #
# master's thesis (A.3) - (A.5) #
#################################
 
I = complex(0,1)

def four_identity(i,j):

    return tensor([tensor([tensor([ 1 ], None), tensor([ 0 ], None), tensor([ 0 ], None), tensor([ 0 ], None)], j),
                   tensor([tensor([ 0 ], None), tensor([ 1 ], None), tensor([ 0 ], None), tensor([ 0 ], None)], j),
                   tensor([tensor([ 0 ], None), tensor([ 0 ], None), tensor([ 1 ], None), tensor([ 0 ], None)], j),
                   tensor([tensor([ 0 ], None), tensor([ 0 ], None), tensor([ 0 ], None), tensor([ 1 ], None)], j)], i)

def mink_metric(i,j):

    return tensor([tensor([tensor([ 1 ], None), tensor([ 0 ], None), tensor([ 0 ], None), tensor([ 0 ], None)], j),
                   tensor([tensor([ 0 ], None), tensor([-1 ], None), tensor([ 0 ], None), tensor([ 0 ], None)], j),
                   tensor([tensor([ 0 ], None), tensor([ 0 ], None), tensor([-1 ], None), tensor([ 0 ], None)], j),
                   tensor([tensor([ 0 ], None), tensor([ 0 ], None), tensor([ 0 ], None), tensor([-1 ], None)], j)], i)
    
def gamma_0(i,j):

    return tensor([tensor([tensor([ 0 ], None), tensor([ 0 ], None), tensor([ 1 ], None), tensor([ 0 ], None)], j),
                   tensor([tensor([ 0 ], None), tensor([ 0 ], None), tensor([ 0 ], None), tensor([ 1 ], None)], j),
                   tensor([tensor([ 1 ], None), tensor([ 0 ], None), tensor([ 0 ], None), tensor([ 0 ], None)], j),
                   tensor([tensor([ 0 ], None), tensor([ 1 ], None), tensor([ 0 ], None), tensor([ 0 ], None)], j)], i)

def gamma_1(i,j):

    return tensor([tensor([tensor([ 0 ], None), tensor([ 0 ], None), tensor([ 0 ], None), tensor([ 1 ], None)], j),
                   tensor([tensor([ 0 ], None), tensor([ 0 ], None), tensor([ 1 ], None), tensor([ 0 ], None)], j),
                   tensor([tensor([ 0 ], None), tensor([-1 ], None), tensor([ 0 ], None), tensor([ 0 ], None)], j),
                   tensor([tensor([-1 ], None), tensor([ 0 ], None), tensor([ 0 ], None), tensor([ 0 ], None)], j)], i)

def gamma_2(i,j):

    return tensor([tensor([tensor([ 0 ], None), tensor([ 0 ], None), tensor([ 0 ], None), tensor([-I ], None)], j),
                   tensor([tensor([ 0 ], None), tensor([ 0 ], None), tensor([ I ], None), tensor([ 0 ], None)], j),
                   tensor([tensor([ 0 ], None), tensor([ I ], None), tensor([ 0 ], None), tensor([ 0 ], None)], j),
                   tensor([tensor([-I ], None), tensor([ 0 ], None), tensor([ 0 ], None), tensor([ 0 ], None)], j)], i)

def gamma_3(i,j):

    return tensor([tensor([tensor([ 0 ], None), tensor([ 0 ], None), tensor([ 1 ], None), tensor([ 0 ], None)], j),
                   tensor([tensor([ 0 ], None), tensor([ 0 ], None), tensor([ 0 ], None), tensor([-1 ], None)], j),
                   tensor([tensor([-1 ], None), tensor([ 0 ], None), tensor([ 0 ], None), tensor([ 0 ], None)], j),
                   tensor([tensor([ 0 ], None), tensor([ 1 ], None), tensor([ 0 ], None), tensor([ 0 ], None)], j)], i)

def gamma_5(i,j):
    return tensor([tensor([tensor([-1 ], None), tensor([ 0 ], None), tensor([ 0 ], None), tensor([ 0 ], None)], j),
                   tensor([tensor([ 0 ], None), tensor([-1 ], None), tensor([ 0 ], None), tensor([ 0 ], None)], j),
                   tensor([tensor([ 0 ], None), tensor([ 0 ], None), tensor([ 1 ], None), tensor([ 0 ], None)], j),
                   tensor([tensor([ 0 ], None), tensor([ 0 ], None), tensor([ 0 ], None), tensor([ 1 ], None)], j)], i)


#######################################
# elementary lorentz structures from  # 
# table 6, arXiv:1108.2040v2 [hep-ph] #
#######################################

class C(tensor):
    
    def __init__(self,i,j):
        array = [tensor([tensor([ 0], None), tensor([1], None), tensor([0], None), tensor([ 0], None)], j),
                 tensor([tensor([-1], None), tensor([0], None), tensor([0], None), tensor([ 0], None)], j),
                 tensor([tensor([ 0], None), tensor([0], None), tensor([0], None), tensor([-1], None)], j),
                 tensor([tensor([ 0], None), tensor([0], None), tensor([1], None), tensor([ 0], None)], j)]

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
        array = [tensor([tensor([ 1], None), tensor([ 0], None), tensor([ 0], None), tensor([ 0], None)], jl),
                 tensor([tensor([ 0], None), tensor([-1], None), tensor([ 0], None), tensor([ 0], None)], jl),
                 tensor([tensor([ 0], None), tensor([ 0], None), tensor([-1], None), tensor([ 0], None)], jl),
                 tensor([tensor([ 0], None), tensor([ 0], None), tensor([ 0], None), tensor([-1], None)], jl)]

        super(Metric,self).__init__(array,il)

class P(tensor):
    
    def __init__(self,i,n):
        # our inde scheme starts at 0, so subtract 1 from n
        # when naming the momenta
        k = n-1
        il = lorentz_key(i) if (i<0) else i
        p_str = "p{0}[{1}]"
        array = [tensor([c_variable(p_str.format(k,0))], None), 
                 tensor([c_variable(p_str.format(k,1))], None), 
                 tensor([c_variable(p_str.format(k,2))], None), 
                 tensor([c_variable(p_str.format(k,3))], None)]

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
        array = [tensor([tensor([ 1], None), tensor([ 0], None), tensor([ 0], None), tensor([ 0], None)], j),
                 tensor([tensor([ 0], None), tensor([ 1], None), tensor([ 0], None), tensor([ 0], None)], j),
                 tensor([tensor([ 0], None), tensor([ 0], None), tensor([ 1], None), tensor([ 0], None)], j),
                 tensor([tensor([ 0], None), tensor([ 0], None), tensor([ 0], None), tensor([ 1], None)], j)]

        super(Identity,self).__init__(array,i)

class Epsilon(tensor):

    def __init__(self,i,j,k,l):

        perms = list(permutations([0,1,2,3]))
        
        il = lorentz_key(i) if (i<0) else i
        jl = lorentz_key(j) if (j<0) else j
        kl = lorentz_key(k) if (k<0) else k
        ll = lorentz_key(l) if (l<0) else l
        tmp = new(il, 4,{jl:4,kl:4,ll:4})

        for perm in perms:
            tmp.__setitem__({il:perm[0], jl:perm[1], kl:perm[2], ll:perm[3]}, 
                            tensor([parity(list(perm))], None))
        
        super(Epsilon,self).__init__(tmp._array, tmp._toplevel_key)
        

def parity(lst):
    parity = 1
    for i in range(0,len(lst)-1):
        if lst[i] != i:
            parity *= -1
            mn = min(range(i,len(lst)), key=lst.__getitem__)
            lst[i],lst[mn] = lst[mn],lst[i]
    return parity    
