from sympy.printing.cxx import CXX11CodePrinter
from .lorentz_simplify import simplify_tensor
import math
import re


VECT_GAUGE_DICT = {
    0: "0",
    1: "ATOOLS::Spinor<SType>::R1()",
    2: "ATOOLS::Spinor<SType>::R2()",
    3: "ATOOLS::Spinor<SType>::R3()",
}


def squarerepl(match):
    return match.group(1) + '*' + match.group(1)


class LorentzPrinter(CXX11CodePrinter):
    def __init__(self, settings=None):
        super().__init__(settings)
        self.index = 0
        self.spin = 0
        self.form_factor = None
        self.form_factor_name = ""
        self.regex = re.compile(r'([^\*]+)\*\*2')
        self.numeric = None

    def set_form_factor(self, form_factor):
        if not isinstance(form_factor, (int, float)):
            self.form_factor_name = form_factor.name
            args = [simplify_tensor(arg) for arg in form_factor.args]
            self.form_factor = f'p_ff->FF({", ".join([str(arg) for arg in args])})'
            self.form_factor = self._print(self.form_factor)
            self.form_factor = self.form_factor.replace('P', 'p')
            self.form_factor = self.form_factor.replace('[', '').replace(']', '')
            self.form_factor = self.regex.sub(squarerepl, self.form_factor)
        else:
            self.form_factor = None

    def set_numeric(self, numeric):
        self.numeric = numeric

    def set_index_spin(self, index, spin):
        self.index = index
        self.spin = spin

    def _sanitize(self, result):
        result = str(result).replace('[', '').replace(']', '')
        result = result.replace('P', 'p')
        return self.regex.sub(squarerepl, result)

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

    def _print_ArrayContraction(self, expr):
        if self.numeric is None:
            return self._get_result(simplify_tensor(expr))
        return self._get_result(simplify_tensor(expr, **self.numeric))

    def _print_ArrayAdd(self, expr):
        if self.numeric is None:
            return self._get_result(sum([simplify_tensor(arg)
                                         for arg in expr.args]))
        return self._get_result(sum([simplify_tensor(arg, **self.numeric)
                                     for arg in expr.args]))

    def _print_ArrayTensorProduct(self, expr):
        if self.numeric is None:
            return self._get_result(math.prod([simplify_tensor(arg)
                                               for arg in expr.args]))
        return self._get_result(math.prod([simplify_tensor(arg, **self.numeric)
                                           for arg in expr.args]))
