from ufo_exception import ufo_exception

class color_function:
    
    def __init__(self, ufo_string):

        # products of elementary color functions can occur: f(-1,1,2)*f(3,4,-1)
        self.multiplicands = [elementary_color_function(string) for string in ufo_string.split("*")]
        if not len(self.multiplicands) in [1,2]:
            raise not_implemented("Cannot handle UFO color function \"{0}\"".format(ufo_string))

        if len(self.multiplicands) == 2:

            uncontracted = ([i for i in self.multiplicands[0].arguments.values() if i > 0] + [i for i in self.multiplicands[1].arguments.values() if i > 0])
            uncontracted.sort()
            assert range(1,len(uncontracted)+1) == uncontracted
            assert self.multiplicands[0].contraction_positions.keys() == self.multiplicands[1].contraction_positions.keys()

            # now replace the summation indices (negative in UFO format) in both color functions with positive integers
            newindex = len(uncontracted)+1
            for key,position1 in self.multiplicands[0].contraction_positions.items():
                position2 = self.multiplicands[1].contraction_positions[key]
                self.multiplicands[0].arguments[position1] = newindex
                self.multiplicands[1].arguments[position2] = newindex
                newindex += 1

    def sherpa_constructor(self):
        if len(self.multiplicands) == 1:
            return self.multiplicands[0].sherpa_constructor()

        return self.multiplicands[0].sherpa_constructor().rstrip(")")+",new "+self.multiplicands[1].sherpa_constructor() + ")"
        
class elementary_color_function:
    
    # expected format:
    # ----------------
    # 1
    # Identity(1,2)
    # T(1,2,3)
    # ...

    def __init__(self, ufo_string):

        self.elementary_dict = {
            "1"          : "cf::None",
            "Identity"   : "cf::D",
            "T"          : "cf::T",
            "f"          : "cf::F",
            #    "d"          : "",
            #    "Epsilon"    : "",
            #    "EpsilonBar" : "", 
            #    "T6"         : "",
            #    "K6"         : "",
            #    "K6Bar"      : ""
            }
        
        self.arguments = {}
        self.contraction_positions = {}
        
        splitted = ufo_string.split("(")
        self.key = splitted[0]
        if not self.key in self.elementary_dict.keys():
            raise not_implemented("UFO color function \"{0}\" not implemented".format(self.key))

        if self.key == "1": 
            assert len(splitted) == 1
            return
        assert len(splitted) == 2

        ints = [int(string) for string in (splitted[1]).rstrip(")").split(",")]
        self.arguments = { i:j for i,j in zip(range(len(ints)),ints)}
        self.contraction_positions = { j:i for i,j in zip(range(len(ints)),ints) if j<0}

    def sherpa_constructor(self):
        ret = "Color_Function("+self.elementary_dict[self.key]
        nargs = len(self.arguments.keys())
        missing = 3 - nargs
        for ind in range(nargs):
            ret += ","+str(self.arguments[ind]-1)
        for i in range(missing):
            ret += ",-1"
        for ind in range(nargs):
            ret += ",\'{0}\'".format(self.arguments[ind]-1)
        for i in range(missing):
            ret += ",\'?\'"
        ret += ")"
        return ret
