import ast
import regex

cmath_dict = {
    "cos": "cos",
    "sin": "sin",
    "tan": "tan",
    "acos": "acos",
    "asin": "asin",
    "atan": "atan",
    "sqrt": "sqrt",
    "pi": "M_PI",
    "log": "log",
    "exp": "exp"
}

color_dict = {
    "1": "Color_Function(cf::None",
    "Identity": "Color_Function(cf::D",
    "IdentityG": "Color_Function(cf::G",
    "T": "Color_Function(cf::T",
    "f": "Color_Function(cf::F",
    "d": "Color_Function(cf::d",
    # Color Sextets are not supported
    #     "Epsilon": "Color_Function(cf::Eps",
    #     "EpsilonBar": "Color_Function(cf::EpsBar",
    #     "T6": "Color_Function(cf::T6",
    #     "K6": "Color_Function(cf::K6",
    #     "K6Bar": "Color_Function(cf::K6Bar",
}


def ensuring_matching_parens(string):
    while string.count('(') != string.count(')'):
        string += ')'
    return string


class param_visitor(ast.NodeVisitor):
    def __init__(self):
        self.value = ""

    def __str__(self):
        return self.value

    def visit_Module(self, node):
        self.generic_visit(node)

    def visit_Expr(self, node):
        self.value = self.visit(node.value)

    def visit_Call(self, node):
        call = self.visit(node.func)
        call += "("
        if len(node.args) > 0:
            call += ",".join([self.visit(arg) for arg in node.args])
        call += ")"
        return call

    def visit_Name(self, node):
        return str(node.id)

    def visit_Attribute(self, node):
        return cmath_dict[node.attr]

    def visit_Constant(self, node):
        return str(float(node.value))

    def visit_UnaryOp(self, node):
        return f"({self.visit(node.op)}{self.visit(node.operand)})"

    def visit_BinOp(self, node):
        if isinstance(node.op, ast.Pow):
            return f"pow({self.visit(node.left)}, {self.visit(node.right)})"
        else:
            left = self.visit(node.left)
            op = self.visit(node.op)
            right = self.visit(node.right)
            return f"({left} {op} {right})"

    def visit_IfExp(self, node):
        test = self.visit(node.test)
        body = self.visit(node.body)
        orelse = self.visit(node.orelse)
        return f"(({test}) ? ({body}) : ({orelse}))"

    def visit_Mult(self, _):
        return '*'

    def visit_Add(self, _):
        return '+'

    def visit_UAdd(self, _):
        return '+'

    def visit_Sub(self, _):
        return '-'

    def visit_USub(self, _):
        return '-'

    def visit_Div(self, _):
        return '/'


class color_visitor(ast.NodeVisitor):
    def __init__(self):
        self.value = ""

    def __str__(self):
        return self.value

    def visit_Module(self, node):
        self.generic_visit(node)

    def visit_Expr(self, node):
        self.value = self.visit(node.value)

    def visit_Call(self, node):
        call = self.visit(node.func)
        call += ","
        if len(node.args) > 0:
            call += ",".join([self.visit(arg) for arg in node.args])
        call += ")"
        return call

    def visit_Name(self, node):
        return color_dict[node.id]

    def visit_Attribute(self, node):
        return color_dict[node.attr]

    def visit_Constant(self, node):
        return str(node.value)

    def visit_UnaryOp(self, node):
        return f"{self.visit(node.op)}{self.visit(node.operand)}"

    def visit_BinOp(self, node):
        if not isinstance(node.op, ast.Mult):
            raise
        left = self.visit(node.left)
        right = self.visit(node.right)
        return f"{left[:-1]}, new {right}"

    def visit_USub(self, _):
        return '-'


def calc_parameter(parameter):
    visitor = param_visitor()
    visitor.visit(ast.parse(parameter))
    return ensuring_matching_parens(visitor.value)


def calc_color(color):
    if color == '1':
        return 'Color_Function(cf::None)'
    visitor = color_visitor()
    visitor.visit(ast.parse(color))
    return ensuring_matching_parens(visitor.value)
