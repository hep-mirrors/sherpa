from py_to_cpp import py_to_cpp

class s_coupling():

    def __init__(self, ufo_coupling):
        self.ufo_coupling = ufo_coupling
    
    def cpp_value(self):
        return py_to_cpp(self.ufo_coupling.value)

    def name(self):
        return self.ufo_coupling.name

    def order(self, string):
        if string in self.ufo_coupling.order.keys():
            return self.ufo_coupling.order[string]
        else:
            return 0
    
