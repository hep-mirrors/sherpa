from copy import deepcopy, copy
from ufo_exception import ufo_exception
from py_to_cpp import c_string_from_num

class tensor(object):
    
    def __init__(self, array, toplevel_key = None):
        self._toplevel_key = toplevel_key
        self._array = array
        self._toplevel_dim = len(self._array)
        self._elementary = True if (self._toplevel_key is None) else False
        if self._elementary:
            assert(len(self._array)) == 1
        #self.check()

    def __repr__(self):
        return "{0} {1} {2}".format(self.__class__ , self._toplevel_key, (self._array).__repr__())

    def __str__(self, indent=""):

        if self._elementary:
            return self._array[0].__str__()

        ret = ""
        if not self._array[0]._elementary:
            for tens,i in zip(self._array, range(self._toplevel_dim)):
                title = indent+"{0}:{1}".format(self._toplevel_key, i)
                sub_indent = " "*len(title)
                ret += title+"\n"+tens.__str__(sub_indent)
        else:
            title = indent+"{0}: ".format(self._toplevel_key)
            ret += title
            for tens in self._array:
                ret += tens.__str__()+" || "
        ret += "\n"
        return ret

    def __eq__(self, rhs):

        # end of recursion: elementary tensor
        if self._elementary:
            return (rhs._elementary and (self._array == rhs._array))

        # toplevel key must be key of rhs
        # and corresp. dim. must be identical
        rhs_kdd = rhs.key_dim_dict()
        if not self._toplevel_key in rhs_kdd:
            return False
        if not self._toplevel_dim == rhs_kdd[self._toplevel_key]:
            return False

        # now check recursively
        for i in range(self._toplevel_dim):
            if (self.__getitem__({self._toplevel_key:i}) != rhs.__getitem__({self._toplevel_key:i})):
                return False
        return True

    # operator '!=' must be implemented separately
    def __ne__(self, rhs):
        return not self.__eq__(rhs)

    # indices is a dict of 'type index_key:index_value'
    def __getitem__(self, indices):

        ret = self

        # Avoid recursion as far as possible by dereferencing toplevel
        # until toplevel key is not in indices anymore
        while(ret._toplevel_key in indices):
            ret = ret._array[indices[ret._toplevel_key]]

        # Trivial case
        if ret._elementary:
            return ret

        # Now have to rely on recursion
        return tensor([tens[indices] for tens in ret._array], ret._toplevel_key)

    def __setitem__(self, indices, value):
        
        # which element to alter?
        tens = self.__getitem__(indices)

        # easy case: elementary tensor, reinit with value
        if tens._elementary:
            tens.__init__(value._array, value._toplevel_key)

        else:
            # more complicated: tens=eta_{mu,...}, value=gamma_{mu,theta,...}
            # all indices in eta must also appear in gamma with same 
            # dimensionality, implicitly done by common_key_dict()
            ccd = common_key_dict(tens, value)
            key, dim = ccd.popitem()
            for i in range(dim):
                tens.__getitem__({key:i}).__setitem__({},value.__getitem__({key:i}))

    def copy(self):
        """Return a deep copy"""
        if self._elementary:
            return tensor([deepcopy(self._array[0])], None)
        return tensor([el.copy() for el in self._array], self._toplevel_key)

    def keys(self):
        """Return list of all keys in the tensor"""
        if self._elementary:
            return []
        ret = [self._toplevel_key]
        ret.extend(self._array[0].keys())
        return ret

    def key_dim_dict(self):
        """Return dictionary mapping the keys
        to the corresponding dimensionality"""
        if self._elementary:
            return {}
        ret = {self._toplevel_key:self._toplevel_dim}
        ret.update(self._array[0].key_dim_dict())
        return ret

    def key_index_dict(self,n=0):
        """Return dictionary mapping the keys
        to the depth in the array where it occurs"""
        if self._elementary:
            return {}
        ret = {self._toplevel_key:n}
        ret.update(self._array[0].key_index_dict(n+1))
        return ret

    def index_key_dict(self,n=0):
        """Return dictionary mapping the depth
        in the array to the corresponding key"""
        if self._elementary:
            return {}
        ret = {n:self._toplevel_key}
        ret.update(self._array[0].index_key_dict(n+1))
        return ret

    def __add__(self, rhs):
        # so far, support/define
        # only sum of identical type
        # i.e. demand equality of key_dim_dict()
        kdd = self.key_dim_dict()
        if (cmp(kdd, rhs.key_dim_dict())!=0):
            raise ufo_exception("Inconsistent tensor addition")
        
        if self._elementary:
            ret = tensor([0.0], None)
            if isinstance(self._array[0], str) or isinstance(rhs._array[0], str):
                ret._array[0] = str(self._array[0]) + " + " + str(rhs._array[0])
            else:
                ret._array[0] = self._array[0] + rhs._array[0]
            return ret

        # Build a return value
        ret = new(kdd)

        # Avoid recursion and loop over all elements of the return
        # tensor explicitly
        for el, ind in self.elements_indices():
            #assert(ret.__getitem__(ind)._elementary)
            ret.__getitem__(ind)._array[0] = el._array[0] + rhs.__getitem__(ind)._array[0]
        return ret

    def __iadd__(self, rhs):
        self = self+rhs
        return self

    def __sub__(self, rhs):
        return (self + tensor([-1], None)*rhs)

    def __mul__(self, rhs):
        
        if not isinstance(rhs, tensor):
            return self.__mul__(tensor([rhs], None))

        # If no indices to be summed over: return simple product
        if len(common_keys(self,rhs)) == 0:
            return multiply(self, rhs)

        # Else: perform contraction of repeated indices
        return contract(self, rhs)

    def __rmul__(self, lhs):
        return self.__mul__(lhs)

    def __div__(self, rhs):
        if isinstance(rhs, tensor):
            if rhs._elementary:
                return self.__mul__(tensor([1.0/rhs._array[0]], None))
        return self.__mul__(1.0/rhs)
    
    def __neg__(self):
        return tensor([-1], None)*self

    def __pow__(self, exp):
        assert(isinstance(exp, int))
        assert(exp>0)
        ret = self.copy()
        for i in range(exp-1):
            ret *= self
        return ret

    def multipliers(self):
        if self._elementary:
            yield self.in_place_elementary_multiply
            return
        for i in range(self._toplevel_dim):
            for mul in self._array[i].multipliers():
                yield mul

    def elements(self):
        """Generator providing all elementary elements of the tensor"""
        if self._elementary:
            yield self
            return
        for i in range(self._toplevel_dim):
            for j in self._array[i].elements():
                yield j

    def elements_indices(self):
        """Generator providing all elementary elements of the tensor along
        with the corresponding index assignment for each element
        """
        if self._elementary:
            yield self, dict()
            return
        for ind in all_indices(self.key_dim_dict()):
            yield self.__getitem__(ind), ind

    def in_place_elementary_multiply(self, other):
        """In-place multiplication of an elementary tensor by a number or by
        another tensor. This is used indirectly by the __mul__ method.
        """
        if not self._elementary:
            raise ValueError("Can only call this method on elementary tensors")

        # If other is just a number, scale up
        if not isinstance(other, tensor):
            self._array[0] = other*self._array[0]

        else:
            rt  = self*other
            self.__init__(rt._array, rt._toplevel_key)

    def check(self):
        # elementary tensors don't have tensors as members of _array
        if not self._elementary:
            # check if _toplevel_dim of all members of _array match
            if len(self._array)>0:
                dim = self._array[0]._toplevel_dim
                key = self._array[0]._toplevel_key
                kdd = self._array[0].key_dim_dict()
            for tens in self._array[1:]:
                if tens._toplevel_dim != dim:
                    raise RuntimeError("Inconsistent tensor")
                if tens._toplevel_key != key:
                    raise RuntimeError("Inconsistent tensor")
                if (cmp(tens.key_dim_dict(), kdd)!=0):
                    raise RuntimeError("Inconsistent tensor")

###################
# unbound functions
###################

def all_indices(key_dim_dict):
    """Argument of this function is a dictionary assigning positive
    integers to keys. Function returns all possible dictionaries
    mapping the keys in the argument dictionary to all possible
    combination of integers smaller than the integer in the argument
    dictionary.
    """

    if len(key_dim_dict)==0:
        yield dict()
        return

    key, dim = key_dim_dict.popitem()
    for rec in all_indices(key_dim_dict):
        for i in range(dim):
            ret = rec.copy()
            ret[key] = i
            yield ret


def multiply(tens_a,tens_b):
    """This method is just a helper function
    for the __mul__ method, no common keys
    should appear, when this is called.
    This implementation is faster by a factor of 10
    compared to the one commented out above."""

    assert(len(common_keys(tens_a,tens_b))==0)
    if (tens_a._elementary and  tens_b._elementary):
        return tensor([tens_a._array[0]*tens_b._array[0]], None)
    ret=(tens_b).copy()
    for mul in ret.multipliers():
        mul(tens_a)
    return ret

def new(key_dim_dict):
    """Create a new tensor with 'tensor([0], None)' as elementary entries.
    Implemented non-recursively for performance reasons

    """

    ret = tensor([0.0], None)
    while key_dim_dict:
        key, dim = key_dim_dict.popitem()
        arr = [ret.copy() for i in range(dim)]
        ret.__init__(arr, key)
    return ret

def contract(tens_a, tens_b):
    """Perform contraction of repeated keys in the two tensors. Avoid
    recursion for the sake of performance.
    """

    ckd   = common_key_dict(tens_a, tens_b)
    kdd   = tens_a.key_dim_dict()
    kdd.update(tens_b.key_dim_dict())
    kdd   = {key:dim for key,dim in kdd.iteritems() if key not in ckd}
    
    ret = new(kdd)

    for el,indices in ret.elements_indices():
        for sum_index in all_indices(ckd.copy()):
            inds = indices.copy()
            inds.update(sum_index)

            # Fuckin ugly hack to fix implicit Minkowski metric
            # insertions in UFO
            pf = 1.0
            for k,i in sum_index.iteritems():
                if isinstance(k, lorentz_key) and i!=0:
                    pf *= -1.0
                    
            el._array[0] += pf * tens_b.__getitem__(inds)._array[0] * tens_a.__getitem__(inds)._array[0]
            
    return ret

def common_keys(tens_a, tens_b):
    return [key for key in tens_a.keys() if key in tens_b.keys()]

def common_key_dict(tens_a, tens_b):
    comm_keys = common_keys(tens_a,tens_b)
    new_dict = dict()
    dict_a = tens_a.key_dim_dict()
    dict_b = tens_b.key_dim_dict()
    for key in comm_keys:
        dim_a = dict_a[key]
        if (dim_a != dict_b[key]):
            raise ufo_exception("Cannot create common tensor key dictionary")
        new_dict.update( {key:dim_a} )
    return new_dict

class lorentz_key(object):
    """A class to acommodate implicit raising/lowering
    of lorentz indices when contracting, as required
    by UFO"""
    
    def __init__(self, key):
        self._key = key

    def __eq__(self, rhs):
        if isinstance(rhs, lorentz_key):
            return (self._key == rhs._key)
        return False

    def __hash__(self):
        return self._key.__hash__()

    def __ne__(self, rhs):
        return (not self.__eq__(rhs))

    def __repr__(self):
        return self._key.__repr__()

    def __str__(self):
        return self._key.__str__()

class color_key(object):
    """A class to acommodate distinction between adjoing, fundamental, and
    anti-fundamental indices, as required by Comix
    """
    
    def __init__(self, key, rep, mapped=None):
        # fu: fundamental
        # af: anti-fundamental
        # ad: adjoint
        if not rep in ['fu','af','ad']:
            raise ValueError("Unknown representation of color_key {0}".format(rep))
        self.rep        = rep
        self.key        = key
        self.mapped_key = mapped if (mapped is not None) else self

    def __eq__(self, rhs):
        if isinstance(rhs, color_key):
            return (self.key == rhs.key)
        return False

    def __hash__(self):
        return self.key.__hash__()

    def __ne__(self, rhs):
        return (not self.__eq__(rhs))

    def __repr__(self):
        return self.key.__repr__()

    def __str__(self):
        return self.key.__str__()

    def mapped_key(self):
        return self.mapped_key
