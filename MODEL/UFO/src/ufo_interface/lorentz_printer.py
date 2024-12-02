from sympy.printing.cxx import CXX11CodePrinter
from .lorentz_simplify import simplify_tensor
import numpy as np
import sympy


VECT_GAUGE_DICT = {
    0: "0",
    1: "ATOOLS::Spinor<SType>::R1()",
    2: "ATOOLS::Spinor<SType>::R2()",
    3: "ATOOLS::Spinor<SType>::R3()",
}


class LorentzPrinter(CXX11CodePrinter):
    def __init__(self, settings=None):
        super().__init__(settings)
        self.index = 0
        self.spin = 0
        self.form_factor = None
        self.form_factor_name = ""
        self.numeric = None
        self.handling_ff = False

    def set_form_factor(self, form_factor):
        if not isinstance(form_factor, (int, float)):
            self.handling_ff = True
            numeric = self.numeric
            self.numeric = None
            self.form_factor_name = form_factor.name
            args = [self._print(arg) for arg in form_factor.args]
            self.form_factor = f'p_ff->FF({", ".join([str(arg) for arg in args])})'
            self.form_factor = self._print(self.form_factor)
            self.form_factor = self.form_factor.replace('P', 'p')
            self.form_factor = self.form_factor.replace('[', '').replace(']', '')
            self.numeric = numeric
            self.handling_ff = False
        else:
            self.form_factor = None

    def set_numeric(self, numeric):
        self.numeric = numeric

    def set_index_spin(self, index, spin):
        self.index = index
        self.spin = spin

    def _sanitize(self, result):
        result = self._print(result)
        result = str(result).replace('[', '').replace(']', '')
        return result.replace('P', 'p')

    def _get_result(self, result):
        if not hasattr(result, 'shape') or result.shape == ():
            return f'(*j{self.index}) = {self._sanitize(result)};\n'
        elif len(result.shape) == 1:
            string = ''
            for i, val in enumerate(result):
                if self.spin == 3:
                    index = VECT_GAUGE_DICT[i]
                else:
                    index = i
                val = self._sanitize(val)
                string += f'(*j{self.index})[{index}] = {val};\n'
            return string
        elif len(result.shape) == 2:
            string = ''
            for i, val in enumerate(result):
                for j, subval in enumerate(val):
                    subval = self._sanitize(subval)
                    index = f"{VECT_GAUGE_DICT[i]}*{result.shape[0]}+{j}"
                    string += f'(*j{self.index})[{index}] = {val[j]};\n'
        else:
            raise ValueError(f'Unknown shape {result.shape}')
        return string.replace('P', 'p')

    def _print(self, expr, **kwargs):
        if isinstance(expr, np.ndarray):
            return self._print_numpy_ndarray(expr)
        return super()._print(expr, **kwargs)

    def _print_ImmutableDenseNDimArray(self, expr):
        rank = expr.rank()
        if rank == 0:
            return self._print(expr._array)
        else:
            return self._get_result(expr._array)
        # TODO: Handle rank > 1

    def _print_ArrayContraction(self, expr):
        if self.numeric is None:
            result = simplify_tensor(expr)
        else:
            result = simplify_tensor(expr, **self.numeric)
        return self._print(result)

    def _print_ArrayAdd(self, expr):
        if self.numeric is None:
            result = sympy.sum([simplify_tensor(arg)
                               for arg in expr.args])
        else:
            result = sympy.sum([simplify_tensor(arg, **self.numeric)
                               for arg in expr.args])
        return self._print(result)

    def _print_ArrayTensorProduct(self, expr):
        if self.numeric is None:
            result = sympy.prod([simplify_tensor(arg)
                                for arg in expr.args])
        else:
            result = sympy.prod([simplify_tensor(arg, **self.numeric)
                                for arg in expr.args])
        if isinstance(result, (sympy.Mul, sympy.ImmutableDenseNDimArray)):
            return self._get_result(result)
        return self._print(result)

    def _print_numpy_ndarray(self, expr):
        # Form factors should always have scalar arguments
        if self.handling_ff:
            return self._sanitize(expr.item())
        if expr.shape == ():
            return self._get_result(expr.item())
        return self._get_result(expr)
