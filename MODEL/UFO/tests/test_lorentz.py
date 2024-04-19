""" Testing Lorentz parsing. """

import pytest
from sympy.tensor.array.expressions import ArraySymbol
import sympy
import numpy as np

from ufo_interface import lorentz_algebra
from ufo_interface import lorentz_simplify

# ArraySymbols
P1 = ArraySymbol('P1', (4,))
P2 = ArraySymbol('P2', (4,))
Metric = lorentz_algebra.LORENTZMAPPING["Metric"]
Gamma = ArraySymbol('Gamma', (4, 4, 4))
ProjM = ArraySymbol('ProjM', (4, 4))

# Symbols
L_1 = sympy.Symbol('L_1')
L_3 = sympy.Symbol('L_3')
S_1 = sympy.Symbol('S_1')
S_2 = sympy.Symbol('S_2')
S_3 = sympy.Symbol('S_3')
L_m1 = sympy.Symbol('L_m1')
S_m1 = sympy.Symbol('S_m1')
S_m2 = sympy.Symbol('S_m2')
dummy_L_m1 = sympy.Symbol('dummy_L_m1')


def test_momentum():
    """ Test Lorentz visitor momentum. """
    value, _ = lorentz_algebra.calc_lorentz('P(1, 3)')
    assert str(value) == 'P2[L_1]'

    value, _ = lorentz_algebra.calc_lorentz('P(-1, 3)')
    assert str(value) == 'P2[L_m1]'


def test_lorentz_mapping():
    """ Test Lorentz visitor mapping. """
    value, _ = lorentz_algebra.calc_lorentz('Identity(1, 3)')
    assert str(value) == 'Identity[S_1, S_3]'

    value, _ = lorentz_algebra.calc_lorentz('IdentityL(1, 3)')
    assert str(value) == 'IdentityL[L_1, L_3]'

    value, _ = lorentz_algebra.calc_lorentz('Gamma(3, 2, 1)')
    assert str(value) == 'Gamma[L_3, S_2, S_1]'

    value, _ = lorentz_algebra.calc_lorentz('Gamma5(2, 1)')
    assert str(value) == 'Gamma5[S_2, S_1]'

    value, _ = lorentz_algebra.calc_lorentz('ProjM(2, 1)')
    assert str(value) == 'ProjM[S_2, S_1]'

    value, _ = lorentz_algebra.calc_lorentz('ProjP(2, 1)')
    assert str(value) == 'ProjP[S_2, S_1]'

    value, _ = lorentz_algebra.calc_lorentz('Sigma(4, 3, 2, 1)')
    assert str(value) == 'Sigma[L_4, L_3, S_2, S_1]'

    value, _ = lorentz_algebra.calc_lorentz('C(2, 1)')
    assert str(value) == 'C[S_2, S_1]'

    value, _ = lorentz_algebra.calc_lorentz('Metric(2, 1)')
    assert str(value) == 'Metric[L_2, L_1]'

    value, _ = lorentz_algebra.calc_lorentz('Epsilon(4, 3, 2, 1)')
    assert str(value) == 'Epsilon[L_4, L_3, L_2, L_1]'


def test_form_factor():
    """ Test Lorentz visitor form factor. """
    _, value = lorentz_algebra.calc_lorentz('ff(x)')
    assert str(value) == 'ff(x)'

    with pytest.raises(ValueError) as excinfo:
        lorentz_algebra.calc_lorentz('ff()')

    assert 'No arguments for form factor' in str(excinfo.value)


def test_single_binop():
    """ Test Lorentz visitor for single binary operators. """

    value, _ = lorentz_algebra.calc_lorentz('P(1, 2) + P(1, 3)')
    assert value == P1[L_1] + P2[L_1]

    with pytest.raises(ValueError) as excinfo:
        lorentz_algebra.calc_lorentz('P(2, 2) + P(1, 3)')
    assert 'Different indices' in str(excinfo.value)

    value, _ = lorentz_algebra.calc_lorentz('P(1, 2) - P(1, 3)')
    assert value == P1[L_1] - P2[L_1]

    with pytest.raises(ValueError) as excinfo:
        lorentz_algebra.calc_lorentz('P(2, 2) - P(1, 3)')
    assert 'Different indices' in str(excinfo.value)

    value, _ = lorentz_algebra.calc_lorentz('P(-1, 2) * P(-1, 3)')
    assert value == sympy.Sum(P1[L_m1]*Metric[L_m1, dummy_L_m1]*P2[dummy_L_m1],
                              (dummy_L_m1, 0, 3), (L_m1, 0, 3))

    expr = 'Gamma(-1, 1, -2) * Gamma(-1, -2, 3)'
    value, _ = lorentz_algebra.calc_lorentz(expr)
    expected = Gamma[L_m1, S_1, S_m2]*Gamma[dummy_L_m1, S_m2, S_3]
    expected *= Metric[L_m1, dummy_L_m1]
    expected = sympy.Sum(expected,
                         (dummy_L_m1, 0, 3), (L_m1, 0, 3), (S_m2, 0, 3))
    assert value == expected

    value, _ = lorentz_algebra.calc_lorentz('3 * P(1, 3)')
    assert value == 3.0*P2[L_1]

    value, _ = lorentz_algebra.calc_lorentz('Gamma(3,2,-1)*ProjM(-1,1)')
    assert value == sympy.Sum(Gamma[L_3, S_2, S_m1]*ProjM[S_m1, S_1],
                              (S_m1, 0, 3))

    value, _ = lorentz_algebra.calc_lorentz('P(1, 3) / 3')
    assert value == P2[L_1] / 3.0

    value, _ = lorentz_algebra.calc_lorentz('P(-1, 3)**2')
    assert value == sympy.Sum(P2[L_m1]*Metric[L_m1, dummy_L_m1]*P2[dummy_L_m1],
                              (dummy_L_m1, 0, 3), (L_m1, 0, 3))


def test_multiple_binop():
    expr = 'Gamma(-1, 1, 2) * P(-1, 3) + ProjM(1, 2)'
    value, _ = lorentz_algebra.calc_lorentz(expr)
    sum_term = Gamma[L_m1, S_1, S_2]*Metric[L_m1, dummy_L_m1]*P2[dummy_L_m1]
    sum_term = sympy.Sum(sum_term, (dummy_L_m1, 0, 3), (L_m1, 0, 3))
    expected = sum_term + ProjM[S_1, S_2]
    assert value == expected

    expr = '(Gamma(-1, 1, 2) + Gamma(-1, 1, -2)*ProjM(-2, 2)) * P(-1, 3)'
    value, _ = lorentz_algebra.calc_lorentz(expr)
    expected = sympy.Sum(Gamma[L_m1, S_1, S_m2]*ProjM[S_m2, S_2], (S_m2, 0, 3))
    expected += Gamma[L_m1, S_1, S_2]
    expected *= Metric[L_m1, dummy_L_m1]*P2[dummy_L_m1]
    expected = sympy.Sum(expected, (dummy_L_m1, 0, 3), (L_m1, 0, 3))
    assert value == expected


def test_simplify():
    from sympy.tensor.array.expressions import convert_indexed_to_array

    expr = 'Gamma(-1, 1, -2)*Gamma(-1, -2, 2)'
    value, _ = lorentz_algebra.calc_lorentz(expr)
    value = convert_indexed_to_array(value)
    simplify = sympy.Array(lorentz_simplify.simplify_tensor(value))
    expected = sympy.Array(4*np.eye(4))
    assert simplify == expected

    expr = 'Gamma(-1, 1, -2)*Gamma(-1, -2, 2)'
    value, _ = lorentz_algebra.calc_lorentz(expr)
    simplify, expr, indices = lorentz_simplify.simplify_symbolic(value)
    assert simplify == expected
    assert expr == ArraySymbol('Tensor_0', (4, 4))
    assert indices == [S_1, S_2]

    _, expr, _ = lorentz_simplify.simplify_symbolic(value)
    assert expr == ArraySymbol('Tensor_1', (4, 4))
