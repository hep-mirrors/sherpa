from s_particle import s_particle
from s_coupling import s_coupling
from s_lorentz import s_lorentz
from s_color import s_color
from ufo_exception import ufo_exception
from colour_converter import colour_translate

spin_dict = {0 : "S",
             1 : "F",
             2 : "V",
             4 : "T"}

translate = {"SVV" : "VVS",
             "SFF" : "FFS"}

# UFO vertices may have multiple couplings which
# can might have incompatible coupling orders.
# Split the vertex up into subvertices such that
# each of them has well defined coupling orders.
# Return a list of s_vertex instances.
# Hierarchy is a list of strings, identifying
# the hierarchy of coupling orders, e.g.
# hierarchy = ['QCD','QED',...]
def split_by_orders(ufo_vertex, hierarchy):

    # we want couplings, lorentz, and colour
    # structures in flat list, such that
    # complete cpl structure of vertex 
    # is \sum_i cpl_list[i]*lor_list[i]*col_list[i]
    cpl_list = []
    col_list = []
    lor_list = []
    for col_ind in range(len(ufo_vertex.color)):
        for lor_ind in range(len(ufo_vertex.lorentz)):
            if (col_ind, lor_ind) in ufo_vertex.couplings:
                cpl_list.append(s_coupling(ufo_vertex.couplings[(col_ind, lor_ind)]))
                lor_list.append(s_lorentz (ufo_vertex.lorentz  [lor_ind]) )
                col_list.append(           ufo_vertex.color    [col_ind]  )
                #col_list.append(s_color   (ufo_vertex.color    [col_ind]) )

    ret = []
    assert(len(cpl_list)>0)
    while (len(cpl_list)>0):
        cur_vertex   = s_vertex(ufo_vertex,hierarchy)
        cur_cpl      =  cpl_list.pop(-1)
        cur_cpl_list = [cur_cpl         ]
        cur_lor_list = [lor_list.pop(-1)]
        cur_col_list = [col_list.pop(-1)]
        for i in xrange(len(cpl_list)-1,-1,-1):
            if cur_cpl.ufo_coupling.order == cpl_list[i].ufo_coupling.order:
                cur_cpl_list.append(cpl_list.pop(i))
                cur_lor_list.append(lor_list.pop(i))
                cur_col_list.append(col_list.pop(i))
        cur_vertex.set_coupling_list(cur_cpl_list)
        cur_vertex.set_colour_list(cur_col_list)
        cur_vertex.set_lorentz_list(cur_lor_list)
        ret.append(cur_vertex)

    return ret
        
class s_vertex():
    
    def __init__(self, _ufo_vertex, hierarchy):
        self._ufo_vertex = _ufo_vertex
        self._coupling_list = []
        self._lorentz_list = []
        self._colour_list = []
        # Hierarchy is a list of strings, identifying
        # the hierarchy of coupling orders, e.g.
        # hierarchy = ['QCD','QED',...]
        self._hierarchy = hierarchy

    def set_coupling_list(self, cpl_list):
        self._coupling_list = cpl_list

    def set_lorentz_list(self, lor_list):
        self._lorentz_list = lor_list

    def set_colour_list(self, col_list):
        self._colour_list = col_list

    def lorentz_type(self):
        ret = [s_particle(part).spin_times_two() for part in self._ufo_vertex.particles]
        ret.sort()
        ret = "".join([ spin_dict[val] for val in ret ])
        return ret if not ret in translate.keys() else translate[ret]

    def name(self):
        return self._ufo_vertex.name

    def particles(self):
        return [s_particle(part) for part in self._ufo_vertex.particles]

    # Return a list of coupling orders
    # corresponding to the strings in 'order_ids'
    # that identify the orders relevant in the model,
    # e.g. order_ids = ['QCD','QED',...]
    def orders(self):
        # take an arbitrary element, since coupling
        # orders should be equal for all couplings
        cpl = self.coupling_list()[0]
        return [cpl.order(string) for string in self._hierarchy]

    def has_ghosts(self):
        return any([ part.spin_times_two()<0 for part in self.particles()])

    def has_goldstones(self):
        return any([ part.is_goldstone() for part in self.particles()])

    def coupling_list(self):
        return self._coupling_list
            
    def colour_list(self):
        return self._colour_list
    
    def lorentz_list(self):
        return self._lorentz_list
    

class vertex_collection(object):

    def __init__(self, vertex_list, id_string):
        self._vertex_list = vertex_list
        self._id_string  = id_string


    def implementation_string(self):

        indent =  "\n    "
        string =  indent+"void {0} () {{".format(self._id_string)
        indent += "  "

        for vert in self._vertex_list:
            string += (indent + 
                       "m_v.push_back(Single_Vertex());")
            for part in vert.particles():
                string += (indent + 
                           "m_v.back().AddParticle( ATOOLS::Flavour((kf_code){0},{1}) );"
                           .format(abs(part.kf_code()),0 if part.kf_code()>0 else 1))
            for cpl in vert.coupling_list():
                string += (
                    indent + 
                    "m_v.back().cpl.push_back( ATOOLS::Kabbala(\"{0}\",ComplexConstant(string(\"{0}\"))) );" 
                    .format(cpl.name()))
            for col in  vert.colour_list():
                # stupid '1' needs to be replaced by some string placeholder
                col = 'None()' if col == '1' else col
                string += (indent +
                           "m_v.back().Color.push_back({0});"
                           .format(colour_translate(col)))
                           # "m_v.back().Color.push_back(\"{0}\");"
                           # .format(col.name()))
                
            for lor in  vert.lorentz_list():
                string += (indent +
                           "m_v.back().Lorentz.push_back(\"{0}\");"
                           .format(lor.name()))
            orders = vert.orders()
            string += indent + "m_v.back().order.resize({0});".format(len(orders))
            for i in range(len(orders)):
                string += indent + "m_v.back().order[{0}]    = {1};".format(i,orders[i])
        
        indent = indent[:-2]
        string += indent + "}"
        return string


    def call_string(self):
        return "\n    {0}();".format(self._id_string)

