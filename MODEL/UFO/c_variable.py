from copy import deepcopy
from py_to_cpp import c_string_from_num
from ufo_exception import ufo_exception

class c_variable(object):
    
    def __init__(self, string, pre=1.0):
        self._prefac = pre
        self._string = string
        
    def __mul__(self, other):

        if isinstance(other, c_variable):
            # string can potentially be a sum, insert parantheses
            return c_variable("("+self._string + ") * (" + other._string + ")", self._prefac * other._prefac)

        if isinstance(other, complex) or isinstance(other, float):
            return c_variable(self._string, self._prefac * other)            

        if isinstance(other, int):
            return c_variable(self._string, self._prefac * float(other))
            
        else:
            raise ufo_exception("Cannot perform multiplication of c_variable with {0}".format(type(other)))

    def __str__(self):
        
        if self._prefac == 1.0:
            return self._string
        if self.is_zero():
            return "0.0"
        if self._prefac == -1.0:
            # string can potentially be a sum, insert parantheses
            return "-(" + self._string + ")"
        return c_string_from_num(self._prefac) + " * (" + self._string + ")"

    def is_zero(self):
        return self._prefac == 0.0

    def __rmul__(self, other):
        return self.__mul__(other)

    def __add__(self, other):
        
        if isinstance(other, c_variable):
            # avoid explicit additions of zeros
            if other.is_zero():
                return c_variable(self._string, self._prefac)
            if self.is_zero():
                return c_variable(other._string, other._prefac)
            # make use of associativity if both have same prefac
            if (self._prefac == other._prefac):
                string = self._string + " + " + other._string
                return c_variable(string, self._prefac)
            # make use of associativity if both have same prefac
            if (self._prefac == -other._prefac):
                string = self._string + " - (" + other._string +  ")"
                return c_variable(string, self._prefac)
            string = self.__str__() + " + " + other.__str__()
            return c_variable(string, 1.0)
        else:
            if other == 0:
                return c_variable(self._string, self._prefac)
            string = self.__str__() + " + " + c_string_from_num(other) 
            return c_variable(string, 1.0)

    def __radd__(self, other):
        return self.__add__(other)

    def __iadd__(self, other):
        tmp = self.__add__(other)
        self.__init__(tmp._string, tmp._prefac)


def numerical_prefactors(var_a,var_b):
    common_facs = []
    facs_a = []
    facs_b = deepcopy(var_b._num_facs)
    for fac in var_a._num_facs:
        if fac in facs_b:
            common_facs.append(fac)
            facs_b.remove(fac)
        else:
            facs_a.append(fac)

    return facs_a, facs_b, common_facs
            
    
def symbolic_prefactors(var_a, var_b):
    common_facs = []
    facs_a = []
    facs_b = deepcopy(var_b._sym_facs)
    for fac in var_a._sym_facs:
        if fac in facs_b:
            common_facs.append(fac)
            facs_b.remove(fac)
        else:
            facs_a.append(fac)

    return facs_a, facs_b, common_facs
