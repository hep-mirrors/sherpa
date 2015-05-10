"""
Conversion of python expressions to C++ compatible strings. 
This code is based on 'converter.py', part of the HERWIG++ UFO interface.
Many thanks to David Grellscheid for the permission to use this code.

Does not support custom functions other than the ones in the cmath_dictionary defined below.
Does not support complex atan function (C++11 feature)
"""

import ast
from ufo_exception import ufo_exception

# how to translate cmath python names 
# into names appropriate for C++
cmath_dictionary = {
    "cos": "cos",
    "sin": "sin",
    "tan": "tan",
    "acos": "acos",
    "asin": "asin",
    "atan": "atan",
    "sqrt": "sqrt",
    "pi": "M_PI",
    "log":"log"
}

def py_to_cpp(expr):
    return cpp_visitor().cpp_string(expr)

def c_string_from_num(num):
    # where this is used, we have 'complex' typedef'd
    if isinstance(num, complex):
        if num == 0:
            return "(0.0)"
        return "(complex({0},{1}))".format(num.real,num.imag)
    # do not want integers floating around in generated c code
    if isinstance(num, int):
        return "({0})".format(float(num))
    if isinstance(num, float):
        return "({0})".format(num)
    raise ufo_exception("Can't convert {0}".format(num))


class cpp_visitor(ast.NodeVisitor):

    def __init__(self):
        pass

    def cpp_string(self, expr):
        self.string = ""
        self.vars = set()
        self.visit(ast.parse(expr))
        return self.string

    def generic_visit(self, node):
        raise NotImplementedError("Node of type \"{0}\" is not implemented".format(type(node).__name__))

    def pass_super(self,node):
        super(type(self),self).generic_visit(node)
        
    def visit_Module(self, node):
        self.pass_super(node)

    def visit_Expr(self, node):
        self.pass_super(node)

    def visit_Attribute(self,node):
        if node.value.id != "cmath":
            raise NotImplementedError("Attribute \"{0}\" is not implemented".format(node.value.id))
        self.string += cmath_dictionary[node.attr]
                
    def visit_UnaryOp(self,node):
        self.string += "("
        self.visit(node.op)
        self.visit(node.operand)
        self.string += ")"

    def visit_BinOp(self, node):
        if type(node.op) == type(ast.Pow()):
            self.handle_power(node)
        else:
            self.string += "("
            self.visit(node.left)
            self.visit(node.op)
            self.visit(node.right)
            self.string += ")"

    def handle_power(self, node):
        self.string += "pow("
        self.visit(node.left)
        self.string += ","
        self.visit(node.right)
        self.string += ")"

    def visit_Call(self,node):
        self.visit(node.func)
        self.string += "("
        if len(node.args) > 0:
            self.visit(node.args[0])
        for a in node.args[1:]:
            self.string += ","
            self.visit(a)
        self.string += ")"

    def visit_Name(self,node):
        text = str(node.id)
        self.vars.add(text)
        self.string += text

    def visit_Num(self, node):
        # some zeros are encoded as 0j
        self.string += "0.0" if node.n == 0 else str(float(node.n))

    def visit_Mult(self, node):
        self.string += "*"

    def visit_Add(self, node):
        self.string += "+"

    def visit_Sub(self, node):
        self.string += "-"

    def visit_USub(self, node):
        self.string += "-"

    def visit_UAdd(self, node):
        self.string += "+"

    def visit_Div(self, node):
        self.string += "/"
