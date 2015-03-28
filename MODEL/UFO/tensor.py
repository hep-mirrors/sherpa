from copy import deepcopy
from ufo_exception import ufo_exception
from py_to_cpp import c_string_from_num

class tensor(object):
    
    def __init__(self, array, toplevel_key):
        self._toplevel_key = toplevel_key
        self._array = array
        self._toplevel_dim = len(self._array)
        self._elementary = True if (self._toplevel_key is None) else False
        if self._elementary:
            assert(len(self._array)) == 1
        self.check()

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
            title = indent+"{0} ".format(self._toplevel_key)
            ret += title
            for tens in self._array:
                ret += tens.__str__()+" "
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

        # end of recursion: elementary tensor
        if self._elementary:
            return self

        # if indices contains value for toplevel index
        if self._toplevel_key in indices:
            # return corresponding element and pass indices
            return self._array[indices[self._toplevel_key]].__getitem__(indices)
        else:
            return tensor([tens[indices] for tens in self._array], self._toplevel_key)

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

        assert(self[indices]==value)
        self.check()

    # return list of all keys
    def keys(self):
        if self._elementary:
            return []
        else:
            ret = [self._toplevel_key]
            ret.extend(self._array[0].keys())
            return ret

    def key_dim_dict(self):
        if self._elementary:
            return {}
        else:
            ret = {self._toplevel_key:self._toplevel_dim}
            ret.update(self._array[0].key_dim_dict())
            return ret

    def __add__(self, rhs):
        # so far, support/define
        # only sum of identical type
        # i.e. demand equality of key_dim_dict()
        if (cmp(self.key_dim_dict(), rhs.key_dim_dict())!=0):
            raise ufo_exception("Inconsistent tensor addition")
        
        # get a deepcopy so we don't need
        # to build a return value from scratch
        ret = deepcopy(self)

        # end of recursion: elementary tensor
        if ret._elementary:
            if isinstance(ret._array[0], str) or isinstance(rhs._array[0], str):
                ret._array[0] = str(ret._array[0]) + " + " + str(rhs._array[0])
            else:
                ret._array[0] = ret._array[0] + rhs._array[0]
            return ret

        # apply addition recursively to elements in array
        for i in range(self._toplevel_dim):
            ret[{self._toplevel_key:i}] = ret[{self._toplevel_key:i}].__add__(rhs[{self._toplevel_key:i}])

        return ret

    def __iadd__(self, rhs):
        self = self+rhs
        return self

    def __sub__(self, rhs):
        return (self + tensor([-1], None)*rhs)

    def __mul__(self, rhs):
        
        if isinstance(rhs, int) or isinstance(rhs, float) or isinstance(rhs, complex) or isinstance(rhs,str):
            return self.__mul__(tensor([rhs], None))

        if not isinstance(rhs, tensor):
            raise ufo_exception("Tensor multiplication for type {0} not supported".format(type(rhs)))

        # if no indices to be summed over:
        # return simple product
        if len(common_keys(self,rhs)) == 0:
            return multiply(self, rhs)

        return contract(common_key_dict(self, rhs), self, rhs)

    def __rmul__(self, lhs):
        return self.__mul__(lhs)

    def __div__(self, rhs):
        if isinstance(rhs, int):
            return self.__mul__(1/rhs)
        if isinstance(rhs, float) or isinstance(rhs, complex):
            return self.__mul__(1.0/rhs)
        else:
            raise ufo_exception("Tensor division for this type not supported")

    def __neg__(self):
        return tensor([-1], None)*self

    def __pow__(self, exp):
        assert(isinstance(exp, int))
        assert(exp>0)
        ret = deepcopy(self)
        for i in range(exp-1):
            ret *= self
        return ret

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
        
def multiply(tens_a,tens_b):
    # this method is just a helper function
    # for the __mul__ method, no common keys
    # should appear, when this is called
    assert(len(common_keys(tens_a,tens_b))==0)

    if (tens_a._elementary and  tens_b._elementary):
        return tensor([tens_a._array[0]*tens_b._array[0]], None)
        
    new_dict = tens_a.key_dim_dict()
    new_dict.update(tens_b.key_dim_dict())

    dict_a = tens_a.key_dim_dict()
    dict_b = tens_b.key_dim_dict()

    key, dim = new_dict.popitem()
    ret = new(key, dim, new_dict)

    for i in range(dim):
        ret.__setitem__({key:i}, multiply(tens_a[{key:i}],tens_b) if (key in dict_a) else multiply(tens_a,tens_b[{key:i}]))

    return ret

# create a new tensor with 'tensor([0], None)' as
# elementary entries. Since this works revursively, 
# need to provide first toplevel_key and toplevel_dim
# as well as lower level key:dim pairs in dict
def new(top_key, top_dim, key_dim_dict):
        
    # termination of recursion: elementary tensor
    if top_key == None:
        return tensor([0], None)

    # if lower level tensors are elementary, dict will
    # be empty
    if len(key_dim_dict) != 0:
        new_key, new_dim = key_dim_dict.popitem()
    else:
        new_key = None
        new_dim = None

    # now fill array of appropriate dimension
    # recursively
    array = []
    for i in range(top_dim):
        array.append(new(new_key, new_dim, deepcopy(key_dim_dict)))

    return tensor(array, top_key)

# this method performs no checks ans assumes
# all keys in the dict are valid keys of the tensors
def contract(key_dim_dict, tens_a, tens_b):

    if len(key_dim_dict)==0:
        assert(len(common_keys(tens_a,tens_b))==0)
        return multiply(tens_a,tens_b)

    key, dim = key_dim_dict.popitem()

    if dim==0:
        raise ufo_exception("Cannot contract empty tensors")

    # need to make a deepcopy here since
    # 'contract' keeps popping items from dict
    ret = [contract(deepcopy(key_dim_dict), tens_a[{key:i}],tens_b[{key:i}]) for i in range(dim)]

    # built-in support for implicit metric
    # multiplication when contracting lorentz indices
    if hasattr(key, '_lorentz'):
        if key._lorentz:
            return -sum(ret[1:], -ret[0])
    return sum(ret[1:], ret[0])

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


# a class to acommodate implicit raising/lowering
# of lorentz indices when contracting, as required
# by UFO
class lorentz_key(object):
    
    def __init__(self, key):
        self._key = key
        self._lorentz = True

    def __eq__(self, rhs):

        if isinstance(rhs, lorentz_key):
            if not (self._lorentz == rhs._lorentz):
                raise ufo_exception("Internal error, ambiguous lorentz key comparison")
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
