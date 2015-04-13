from py_to_cpp import py_to_cpp

class s_parameter():

    def __init__(self, ufo_parameter):
        assert (ufo_parameter.nature in ["external", "internal"])
        if ufo_parameter.nature == "external":
            assert ufo_parameter.lhacode is not None, "No lhacode for external parameter present"
            assert len(ufo_parameter.lhacode) in [1,2], "Unknown lhacode format \"{0}\"".format(ufo_parameter.lhacode)
        assert ufo_parameter.type in ["real","complex"], "Unknown parameter type \"{0}\"".format(ufo_parameter.type)
        self.ufo_parameter = ufo_parameter
        
    def is_external(self):
        return self.ufo_parameter.nature=="external"

    def is_internal(self):
        return self.ufo_parameter.nature=="internal"

    def is_complex(self):
        return self.ufo_parameter.type=="complex"

    def is_real(self):
        return self.ufo_parameter.type=="real"

    def name(self):
        return self.ufo_parameter.name

    def lha_indices(self):
        return self.ufo_parameter.lhacode

    def lha_block(self):
        return self.ufo_parameter.lhablock
            
    def cpp_value(self):
        return py_to_cpp(str(self.ufo_parameter.value))

    def raw_value(self):
        return self.ufo_parameter.value
    
