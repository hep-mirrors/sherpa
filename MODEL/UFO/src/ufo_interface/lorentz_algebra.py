import ast
from sympy.tensor.array.expressions import ArraySymbol
from sympy.tensor.array.expressions.array_expressions import ArrayElement
from sympy.tensor.array.expressions import convert_indexed_to_array
import sympy
from collections import OrderedDict

METRIC = ArraySymbol('Metric', (4, 4))

LORENTZMAPPING = {'Identity': ArraySymbol('Identity', (4, 4)),
                  'IdentityL': ArraySymbol('IdentityL', (4, 4)),
                  'Gamma': ArraySymbol('Gamma', (4, 4, 4)),
                  'Gamma5': ArraySymbol('Gamma5', (4, 4)),
                  'ProjM': ArraySymbol('ProjM', (4, 4)),
                  'ProjP': ArraySymbol('ProjP', (4, 4)),
                  'Sigma': ArraySymbol('Sigma', (4, 4, 4, 4)),
                  'C': ArraySymbol('C', (4, 4)),
                  'P': ArraySymbol('P', (4,)),
                  'Metric': METRIC,
                  'Epsilon': ArraySymbol('Epsilon', (4, 4, 4, 4))}

INDEX_MAPPING = {'Identity': ['S', 'S'],
                 'IdentityL': ['L', 'L'],
                 'Gamma': ['L', 'S', 'S'],
                 'Gamma5': ['S', 'S'],
                 'ProjM': ['S', 'S'],
                 'ProjP': ['S', 'S'],
                 'Sigma': ['L', 'L', 'S', 'S'],
                 'C': ['S', 'S'],
                 'P': ['L', 'P'],
                 'Metric': ['L', 'L'],
                 'Epsilon': ['L', 'L', 'L', 'L']}


def _get_repeated_indices(indices):
    unique = OrderedDict()
    for idx in indices:
        if idx in unique:
            unique[idx] = 0
        else:
            unique[idx] = 1
    repeated = [i for i, v in unique.items() if v == 0]
    return repeated


def _flatten(lst):
    flattened = []
    for item in lst:
        if isinstance(item, list):
            flattened.extend(_flatten(item))
        else:
            flattened.append(item)
    return flattened


def combine_lists(indices_list):
    """
    Combines a list of lists of strings, keeping only duplicates that exist entirely
    within a single list and excluding repeated elements if they exist in both lists.

    Args:
        indices_list (list): A list of lists of strings.

    Returns:
        list: The combined list of strings.
    """
    combined_list = indices_list[0]
    for indices in indices_list[1:]:
        sub_list = []
        for index in indices:
            if index in combined_list:
                continue
            sub_list.append(index)
        combined_list.extend(sub_list)

    return combined_list


def remove_repeated_indices(indices):
    """
    Removes repeated indices from a list of strings.

    Args:
        indices (list): A list of strings.

    Returns:
        list: The list of strings with repeated indices removed.
    """
    unique = OrderedDict()
    for idx in indices:
        if idx in unique:
            unique[idx] = 0
        else:
            unique[idx] = 1
    return [i for i, v in unique.items() if v == 1]


def all_indices(expr):
    if isinstance(expr, ArrayElement):
        return _flatten(expr.indices)
    elif isinstance(expr, sympy.Mul):
        return _flatten([all_indices(arg) for arg in expr.args])
    elif isinstance(expr, sympy.Add):
        indices = [_flatten(all_indices(arg)) for arg in expr.args]
        return _flatten(combine_lists(indices))
    elif isinstance(expr, sympy.Sum):
        return remove_repeated_indices(all_indices(expr.function))
    else:
        return []


def get_repeated_indices(expr):
    indices = all_indices(expr)
    indices = _get_repeated_indices(indices)
    return indices


def MakeMomentum(args):
    name = f'P{int(args[1])-1}'
    if args[0] > 0:
        index = sympy.Symbol(f'L_{int(args[0])}')
    else:
        index = sympy.Symbol(f'L_m{-int(args[0])}')
    return ArraySymbol(name, (4,))[index]


def _check_indices(left, right):
    lidx = set(all_indices(left))
    ridx = set(all_indices(right))
    if lidx != ridx:
        return False
    return True


def handle_add(left, right, op):
    if not _check_indices(left, right):
        raise ValueError('Different indices')
    if isinstance(op, ast.Add):
        return left+right
    else:
        return left-right


def handle_mult(left, right):
    # Handle scalar multiplication
    if isinstance(left, float) or not isinstance(right, ArrayElement):
        return left*right

    # Collect indices on left and right
    lidx = all_indices(left)
    ridx = all_indices(right)
    repeated = _get_repeated_indices(lidx+ridx)
    if len(repeated) == 0:
        return left*right

    # Add dummy indices for repeated Lorentz indices
    dummy_indices = {index: sympy.Symbol(f'dummy_{index}')
                     for index in repeated
                     if 'L' in str(index)}

    # Handle case where there are no repeated Lorentz indices
    if len(dummy_indices) == 0:
        result = left*right
        for index in repeated:
            result = sympy.Sum(result, (index, 0, 3))
        return result

    # Add metric for repeated Lorentz indices
    result = left
    ridx = [dummy_indices.get(index, index) for index in ridx]
    for index, dummy_index in dummy_indices.items():
        result *= METRIC[index, dummy_index]
        repeated.remove(index)
    result *= right.name[*ridx]
    for index, dummy_index in dummy_indices.items():
        result = sympy.Sum(result, (dummy_index, 0, 3), (index, 0, 3))

    # Sum over repeated indices
    for index in repeated:
        result = sympy.Sum(result, (index, 0, 3))
    return result


class lorentz_visitor(ast.NodeVisitor):
    def __init__(self):
        self.value = None
        self.form_factor = 1

    def to_symbol(self, args, name):
        if isinstance(args, list):
            if name not in INDEX_MAPPING:
                return args
            indices = INDEX_MAPPING[name]
            args = [f'{indices[i]}_{int(arg)}' if arg > 0
                    else f'{indices[i]}_m{-int(arg)}'
                    for i, arg in enumerate(args)]
            args = [sympy.Symbol(arg) for arg in args]
            return args
        else:
            return args

    def __str__(self):
        return self.value

    def visit_Module(self, node):
        self.generic_visit(node)

    def visit_Expr(self, node):
        self.value = self.visit(node.value)

    def visit_Attribute(self, node):
        return node.attr

    def visit_Call(self, node):
        func_name = self.visit(node.func)
        args = [self.visit(arg) for arg in node.args]
        if func_name == 'P':
            return MakeMomentum(args)
        elif func_name == 'complex':
            return args[0]+sympy.I*args[1]
        elif func_name in LORENTZMAPPING:
            args = self.to_symbol(args, func_name)
            func = LORENTZMAPPING[func_name]
            return func[*args]
        else:
            if len(args) == 0:
                raise ValueError(f'No arguments for form factor: {func_name}')
            args = self.to_symbol(args, func_name)
            repeated = get_repeated_indices(args[0])
            for index in repeated:
                args[0] = sympy.Sum(args[0], (index, 0, 3))
            args[0] = convert_indexed_to_array(args[0])
            self.form_factor *= sympy.Function(func_name)(*args)
            return 1

    def visit_Name(self, node):
        return str(node.id)

    def visit_Constant(self, node):
        return float(node.value)

    def visit_UnaryOp(self, node):
        return -self.visit(node.operand)

    def visit_BinOp(self, node):
        left = self.visit(node.left)
        right = self.visit(node.right)
        if isinstance(node.op, (ast.Add, ast.Sub)):
            return handle_add(left, right, node.op)
        elif isinstance(node.op, ast.Mult):
            return handle_mult(left, right)
        elif isinstance(node.op, ast.Pow):
            if int(right) != 2:
                raise ValueError('Only square of a element is supported')
            dummy_index = sympy.Symbol(f'dummy_{left.indices[0]}')
            metric = METRIC[left.indices[0], dummy_index]
            result = left*metric*left.name[dummy_index]
            return sympy.Sum(result, (dummy_index, 0, 3), (left.indices[0], 0, 3))
        elif isinstance(node.op, ast.Div):
            return left/right


def calc_lorentz(lorentz):
    visitor = lorentz_visitor()
    visitor.visit(ast.parse(lorentz))
    return visitor.value, visitor.form_factor
