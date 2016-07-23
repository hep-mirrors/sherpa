from templates import color_calc_template

class s_color(object):

    # Here we store color representations of the coupling structures.
    # Keys are the strings representing the color structure in UFO. Do
    # this, because computing the tensor representations is expensive.
    tensor_cache = dict()
    
    def __init__(self,ufo_color):
        self.ufo_color = ufo_color
        self.unique_id = self.ufo_color.replace('(','_').replace(')','_').replace('-','m').replace(',','_').replace('*','').rstrip('_')

    def name(self):
        return self.unique_id

    def c_name(self):
        return self.unique_id+'.C'

    # get a tensor representation
    def get_cpl_tensor(self):
        if not self.ufo_color in s_color.tensor_cache:
            s_color.tensor_cache[self.ufo_color] = eval(self.ufo_color)
        return s_color.tensor_cache[self.ufo_color]

    def write(self, path):
        with open(path+"/"+self.c_name(), "w") as outfile:
            outfile.write(color_calc_template.substitute(color_name = self.name()))
