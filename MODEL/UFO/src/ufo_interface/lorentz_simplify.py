from sympy.tensor.array.expressions import ArrayTensorProduct, ArraySymbol, \
    ArrayAdd, PermuteDims, ArrayContraction
import sympy
from sympy.tensor.array.expressions import convert_indexed_to_array
import numpy as np
from opt_einsum import contract
from .lorentz_algebra import all_indices

Identity = np.diag([1, 1, 1, 1]).astype(np.float32)
IdentityL = np.diag([1, 1, 1, 1]).astype(np.float32)
Gamma = np.array([[[0, 0, 1, 0],
                   [0, 0, 0, 1],
                   [1, 0, 0, 0],
                   [0, 1, 0, 0]],
                  [[0, 0, 0, 1],
                   [0, 0, 1, 0],
                   [0, -1, 0, 0],
                   [-1, 0, 0, 0]],
                  [[0, 0, 0, -1j],
                   [0, 0, 1j, 0],
                   [0, 1j, 0, 0],
                   [-1j, 0, 0, 0]],
                  [[0, 0, 1, 0],
                   [0, 0, 0, -1],
                   [-1, 0, 0, 0],
                   [0, 1, 0, 0]]])
Gamma5 = np.array([[-1, 0, 0, 0],
                   [0, -1, 0, 0],
                   [0, 0, 1, 0],
                   [0, 0, 0, 1]], dtype=np.float32)
ProjM = np.array([[1, 0, 0, 0],
                  [0, 1, 0, 0],
                  [0, 0, 0, 0],
                  [0, 0, 0, 0]], dtype=np.float32)
ProjP = np.array([[0, 0, 0, 0],
                  [0, 0, 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]], dtype=np.float32)
Sigma = 0.5j*(np.einsum('abc,dce->adbe', Gamma, Gamma)
              - np.einsum('abc,dce->dabe', Gamma, Gamma))
C = np.array([[0, 0, 1, 0],
              [0, 0, 0, -1],
              [1, 0, 0, 0],
              [0, -1, 0, 0]], dtype=np.float32)
Metric = np.diag([1, -1, -1, -1]).astype(np.float32)
Epsilon = -np.array([(a-b)*(a-c)*(a-d)*(b-c)*(b-d)*(c-d)/12
                     for a in range(1, 5)
                     for b in range(1, 5)
                     for c in range(1, 5)
                     for d in range(1, 5)],
                    dtype=np.float32).reshape((4, 4, 4, 4))

# Create mapping from string to numpy array
LORENTZNUMERIC = {'Identity': Identity,
                  'IdentityL': IdentityL,
                  'Gamma': Gamma,
                  'Gamma5': Gamma5,
                  'ProjM': ProjM,
                  'ProjP': ProjP,
                  'Sigma': Sigma,
                  'C': C,
                  'Metric': Metric,
                  'Epsilon': Epsilon}


def create_symbolic_array(name, shape):
    result = np.empty(shape, dtype=object)
    for idx, _ in np.ndenumerate(result):
        result[idx] = sympy.Symbol(f'{name}[{", ".join([str(i) for i in idx])}]')
    return result


def simplify_tensor(arg, **kwargs):
    if isinstance(arg, ArrayContraction):
        base = arg.expr
        indices = arg.contraction_indices

        if isinstance(base, ArrayTensorProduct):
            tensors = [simplify_tensor(arg, **kwargs) for arg in base.args]
            ranks = base.subranks
        else:
            tensors = base
            ranks = [len(base.shape)]

        contraction, letters_free = _get_einsum_string(ranks, indices)
        letters_free = "".join(sorted(letters_free))
        contraction = contraction + "->" + letters_free
        result = contract(contraction, *tensors, backend='object')
        return result
    elif isinstance(arg, ArrayTensorProduct):
        tensors = [simplify_tensor(arg, **kwargs) for arg in arg.args]
        result = tensors[0]
        for tensor in tensors[1:]:
            if isinstance(result, np.ndarray):
                result = np.tensordot(result, tensor, axes=0)
            else:
                result *= tensor
        return result
    elif isinstance(arg, PermuteDims):
        return np.transpose(simplify_tensor(arg.expr, **kwargs),
                            arg.permutation.array_form)
    elif isinstance(arg, ArraySymbol):
        name = str(arg.name)
        shape = arg.shape
        if name in LORENTZNUMERIC:
            return LORENTZNUMERIC[name]
        elif name in kwargs:
            return kwargs[name]
        else:
            return np.array([sympy.Symbol(f'{name}[{i}]')
                             for i in range(np.prod(shape))],
                            dtype=object).reshape(shape)
    elif isinstance(arg, ArrayAdd):
        return sum([simplify_tensor(arg, **kwargs) for arg in arg.args])
    elif isinstance(arg, sympy.Number):
        return float(arg)
    elif isinstance(arg, sympy.core.numbers.ImaginaryUnit):
        return 1j
    elif isinstance(arg, (float, complex, sympy.Symbol)):
        return arg
    else:
        print(arg)
        raise ValueError(f'Unknown type {type(arg)}')


def simplify_symbolic(input):
    tensor = convert_indexed_to_array(input)
    tensor = simplify_tensor(tensor)
    indices = all_indices(input)
    name = NameGenerator.generate()
    if not isinstance(tensor, np.ndarray):
        shape = (1,)
    else:
        shape = tensor.shape
    name = ArraySymbol(name, shape)
    return sympy.Array(tensor), name, indices


def _generate_name():
    i = 0
    while True:
        yield f'Tensor_{i}'
        i += 1


class NameGenerator:
    _generator = _generate_name()

    @classmethod
    def generate(cls):
        return next(cls._generator)


def _get_letter_generator_for_einsum():
    for i in range(97, 123):
        yield chr(i)
    for i in range(65, 91):
        yield chr(i)
    raise ValueError('Ran out of letters')


# Implementation taken from sympy.printing.pycode.ArrayPrinter
def _get_einsum_string(ranks, contraction_indices):
    letters = _get_letter_generator_for_einsum()
    contraction_string = ''
    counter = 0
    d = {j: min(i) for i in contraction_indices for j in i}
    indices = []
    for rank in ranks:
        lindices = []
        for i in range(rank):
            if counter in d:
                lindices.append(d[counter])
            else:
                lindices.append(counter)
            counter += 1
        indices.append(lindices)
    mapping = {}
    letters_free = []
    for i in indices:
        for j in i:
            if j not in mapping:
                letter = next(letters)
                mapping[j] = letter
            else:
                letter = mapping[j]
            contraction_string += letter
            if j not in d:
                letters_free.append(letter)
        contraction_string += ","
    contraction_string = contraction_string[:-1]
    return contraction_string, letters_free
