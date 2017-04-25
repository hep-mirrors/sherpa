"""
Conversion of UFO colour structure strings into
C++ constructor expressions that create corresponding
sherpa Colour_Function instances
"""

cf_dict = {
    # '1' is not a 'name' in the node sense,
    # therefore it has to be replaced by the
    # placeholder string 'None'. 
    # This is done elsewhere.
    'None'     : 'UFO::UFO_CF_1',
    'Identity' : 'UFO::UFO_CF_Identity',
    'T'        : 'UFO::UFO_CF_T',
    'f'        : 'UFO::UFO_CF_f'
}

import ast

def colour_translate(expr):
    return colour_visitor().cpp_string(expr)

class colour_visitor(ast.NodeVisitor):

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

    def visit_Call(self,node):
        self.visit(node.func)
        self.string += "("
        # UFO syntax relies on integers as 
        # arguments of colour structures
        for arg in node.args:
            assert(type(arg).__name__ == 'Num')
        if len(node.args) > 0:
            self.visit(node.args[0])
        for a in node.args[1:]:
            self.string += ","
            self.visit(a)
        self.string += ")"

    def visit_Name(self,node):
        name = str(node.id)
        self.string += cf_dict[name]

    def visit_Num(self, node):
        self.string += str(node.n)

    def visit_Mult(self, node):
        self.string += "*"

    # def visit_Add(self, node):
    #     self.string += "+"

    # def visit_Sub(self, node):
    #     self.string += "-"

    # def visit_USub(self, node):
    #     self.string += "-"

    # def visit_UAdd(self, node):
    #     self.string += "+"

    # def visit_Div(self, node):
    #     self.string += "/"
