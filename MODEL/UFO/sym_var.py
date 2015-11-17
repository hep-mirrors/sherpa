from ufo_exception import ufo_exception
from sympy import ccode, expand, Eq
from sympy.parsing.sympy_parser import parse_expr
from sympy.core.basic import Basic

class sym_var(object):

    def __init__(self, init):
        if isinstance(init, str):
            self._expr = parse_expr(init)
        elif isinstance(init, Basic):
            self._expr = init
        else:
            raise ufo_exception("Expect string or sympy expression in constructor")

    def __mul__(self, other):
        if isinstance(other, sym_var):
            return sym_var(self._expr*other._expr)
        return sym_var(self._expr*other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __add__(self, other):
        if isinstance(other, sym_var):
            return sym_var(self._expr+other._expr)
        return sym_var(self._expr+other)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, sym_var):
            return sym_var(self._expr-other._expr)
        return sym_var(self._expr-other)

    def __rsub__(self, other):
        return self.__sub__(other)

    def __div__(self, other):
        if isinstance(other, sym_var):
            return sym_var(self._expr/other._expr)
        return sym_var(self._expr/other)

    def __rdiv__(self, other):
        return sym_var(other/self._expr)

    def __str__(self):
        return ccode(expand(self._expr))

    def __eq__(self, other):
        return Eq((self-other)._expr, 0.0)
