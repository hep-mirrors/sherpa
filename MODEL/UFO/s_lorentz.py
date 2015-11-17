from copy import deepcopy
from tensor import tensor
from templates import lorentz_calc_template
from code_snippets import *
from lorentz_structures import C,Gamma,Gamma5,Metric,P,ProjM,ProjP,Epsilon,Identity,mink_metric,is_ffv,is_vvv
from ufo_exception import ufo_exception
from sympy import Eq

class s_lorentz():

    # Here we store tensor representations
    # of the coupling structures. Keys are
    # the 'structure' strings in the 
    # corresponding UFO lorentz class.
    # Do this, because computing the
    # tensor representations is expensive.
    tensor_cache = dict()
    
    def __init__(self,ufo_lorentz):
        self.ufo_lorentz = ufo_lorentz
        self._key_spin_dict     = dict()
        self._key_tens_dict     = dict()
        self._ferm_partner_dict = dict()
        for key,spin in enumerate(self.ufo_lorentz.spins):
            self._key_spin_dict[key] = spin
            self._key_tens_dict[key] = get_in_current_tens(key,spin)  if not self.has_ghosts() else None

            # for fermions, find the fermion flow partner, i.e. the next fermion in the list
            if (spin == 2) and not (key in self._ferm_partner_dict.values()):
                for key2 in range(key+1, len(self.ufo_lorentz.spins)):
                    if self.ufo_lorentz.spins[key2] == 2:
                        self._ferm_partner_dict[key] = key2
                        break

        # now assume that UFO vertices are arranged such that _ferm_partner_dict.keys() are 
        # all 'bar' type spinors and _ferm_partner_dict.values() are 'non-bar', i.e.
        # fermion flow always continuous through 'adjacent' fermions in the 'spins' list
        for key in self._ferm_partner_dict.keys():
            assert(key not in self._ferm_partner_dict.values())

        # store fermionic keys for convenience
        self._ferm_keys  = [key for key,spin in self._key_spin_dict.iteritems() if spin==2]
        self._n_ferms    = len(self._ferm_keys)

    def n_ext(self):
        return len(self.ufo_lorentz.spins)

    def write(self, path, ferm_optimize=False):
        with open(path+"/"+self.c_name(), "w") as outfile:
            outfile.write(lorentz_calc_template.substitute(vertex_name = self.calc_class_name(),
                                                           implementation = self.get_all_implementations(ferm_optimize)))

    # map to old hard-wired calculators if possible
    def mapped_name(self):
        if hasattr(self, '_mapped_name'):
            return self._mapped_name
        # is_ffv cannot differentiate between ffv and
        # similar structures with arbitrary scalars
        # attached, so check here for correct spin structure
        if (sorted(self.ufo_lorentz.spins) == sorted([2,2,3])) and is_ffv(self.get_cpl_tensor()):
            self._mapped_name = "FFV"
        # similarly, is_vvv can't tell if we have scalars attached
        elif (self.ufo_lorentz.spins == [3,3,3]) and is_vvv(self.get_cpl_tensor()):
            self._mapped_name = "VVV"
        else:
            self._mapped_name =None
        return self.mapped_name()
        
    # name of the C++ file for a corresponding
    # calculator
    def c_name(self):
        return self.ufo_lorentz.name+".C"

    def name(self):
        return self.ufo_lorentz.name
        
    # name of the C++ calculator class in the
    def calc_class_name(self):
        return self.ufo_lorentz.name

    def has_ghosts(self):
        return any([spin<0 for spin in self.ufo_lorentz.spins])

    def spins(self):
        return self.ufo_lorentz.spins

    # get a tensor representation of the lorentz
    # coupling structure
    def get_cpl_tensor(self):
        if not self.ufo_lorentz.structure in s_lorentz.tensor_cache:
            s_lorentz.tensor_cache[self.ufo_lorentz.structure] = eval(self.ufo_lorentz.structure)
        return s_lorentz.tensor_cache[self.ufo_lorentz.structure]

    def get_all_implementations(self,ferm_optimize):
        imp = ""
        # for each key i, create one calculator
        # corresponding to key i outgoing
        for i in range(self.n_ext()):
            imp += "\n// if outgoind UFO-index is {0}\n".format(i)
            imp += "if (p_v->V()->id.back()=={0}){{\n".format(i)
            imp += self.get_implementation(i,ferm_optimize)
            imp += "\n}\n"
        return imp

    def get_case_dicts(self, keys, dict_list=None):
        if dict_list==None: dict_list = []
        if (len(keys) == 0):
            return dict_list
        key = keys.pop()
        assert(key in self._ferm_keys)
        if (len(dict_list)==0):
            return self.get_case_dicts(keys,[{key:1}, {key:2}, {key:3}])
        else:
            tmp = []
            for dct in dict_list:
                dct.update({key:1})
                tmp.append(deepcopy(dct))
                dct.update({key:2})
                tmp.append(deepcopy(dct))
                dct.update({key:3})
                tmp.append(deepcopy(dct))
            return self.get_case_dicts(keys,tmp)

    def get_key_index_dict(self, out_key):

        # ufo--vertex    rotated version
        # 1  2   3  4      3  4   5  1
        #  \ \   / /        \ \   / /
        #   \ \ / /          \ \ / /
        #     XXX              XXX
        #      |                |
        #      |                |
        #      |                |
        #      5                2
        #                    out_key
        #
        # Via the index keys of the ufo vertex, it is
        # rotiationally invariant. Just need to make sure
        # we give the tensor representation of the
        # external currents the appropriate keys, according
        # to the rotation we're looking at.
        # In the above example, the 0th external current at the rotated
        # vertex must get key '3'. Need a dict {key:index}.
        # Since starting at 0 is more convenient in python and c,
        # rather use a dict {key-1:index} and make sure we add
        # a 1 to each key when constructing a tensor representation.

        num_ext = self.n_ext()
        assert(out_key in range(num_ext))
        # Python 2.7 or later:
        # return { i:(i+(num_ext-1-out_key))%(num_ext) for i in range(num_ext)}
        ret = dict()
        for i in range(num_ext):
            ret[i] = (i+(num_ext-1-out_key))%(num_ext)
        return ret

    # are external momenta required for writing out
    # lorentz coupling structure: determine from UFO string
    def needs_external_momenta(self):
        return "P(" in self.ufo_lorentz.structure

    def get_implementation(self, out_key, ferm_optimize):
        ferm_opt = (len(self._ferm_keys) > 0) if ferm_optimize else False
        assert((not self.has_ghosts()))

        imp = ""

        key_index_dict = self.get_key_index_dict(out_key)
        out_spin       = self._key_spin_dict[out_key]
        out_tens       = self._key_tens_dict[out_key]
        in_keys        = [ key for key in self._key_spin_dict.keys() if key!=out_key ]

        # declare incoming currents and momenta
        for key in in_keys:
            imp += get_in_current_declaration(self._key_spin_dict[key],
                                              key, key_index_dict[key],
                                              1 if (key in self._ferm_partner_dict.keys()) else -1)
            if self.needs_external_momenta():
                imp += get_in_mom_declaration(key, key_index_dict[key])

        # declare and initialize outgoing momentum
        if self.needs_external_momenta():
            imp += get_out_mom_declaration(out_key, key_index_dict)

        # declare the return value, i.e. outgoing current
        imp += get_out_current_declaration(out_spin, out_key)

        # optimized version for massless spinors:
        # can ask Comix via the 'C_Spinor::On()' 
        # method, which components are nonzero:
        # (x,x,0,0)::On() = 1
        # (0,0,x,x)::On() = 2
        # (x,x,x,x)::On() = 3
        # Now create one calculator for each of
        # the possible combinations of incoming 
        # external 'On' values. Use binary rep.
        # of an integer to identify a specific
        # cobination. E.g. three incoming spinors:
        # case_id = spinor_1.On()+spinor_2.On()*4+spinor_3.On()*16

        in_ferm_keys = [key for key in self._ferm_keys if key!=out_key]
        case_dicts = self.get_case_dicts(deepcopy(in_ferm_keys)) if ferm_opt else [{}]
                
        if ferm_opt:
            imp += "switch({0}){{\n".format(self.construct_case_id(in_ferm_keys))

        for case_dict in case_dicts:
            if (len(case_dicts)>1):
                imp += "case {0}:\n".format(self.get_case_id(case_dict))

            # tensor representing coupling struct
            # to be contracted with ext. wavefct.
            cpl  =  deepcopy(self.get_cpl_tensor())

            # contract incoming wavefct's with the cpl tensor
            for key in in_keys:
                tens = deepcopy(self._key_tens_dict[key])

                # for massless spinor optimization: assume
                # only two components are nonzero according 
                # to the value stored in the dict
                if key in case_dict:
                    if case_dict[key] == 1 :
                        tens._array[2]  = tensor([0], None)
                        tens._array[3]  = tensor([0], None)
                    elif case_dict[key] == 2 :
                        tens._array[0]  = tensor([0], None)
                        tens._array[1]  = tensor([0], None)
                
                cpl =  cpl * (tens)
              
            # check if contraction of external tensors with
            # coupling tensor yields desired outgoing tensor type
            assert(cmp(cpl.key_dim_dict(),out_tens.key_dim_dict())==0)

            # is the return tensor arithmetically zero?
            return_zero = self.is_zero(out_spin, cpl)
            
            # initialize the return value
            if not return_zero:
                out_bar_type = 1 if out_key not in self._ferm_partner_dict.keys() else -1
                imp += get_out_current_initialization(out_spin, out_key, cpl, out_bar_type)

            # set the 'S' property
            if not return_zero:
                imp += "j{0}->SetS(j{1}.S()".format(out_key,in_keys[0])
                for k in in_keys[1:]:
                    imp += "|j{0}.S()".format(k)
                imp += ");\n"

            imp += "return j{0};\n".format(out_key)

        # close the switch statement
        if ferm_opt:
            imp += "default:\n THROW(fatal_error, \"Massless spinor optimization error in Lorentz calculator\");\n"
            imp += "}\n"

        return imp

    def construct_case_id(self,in_ferm_keys):
        ret = ""
        for key in in_ferm_keys:
            ret += "+(j{0}.On()<<({1}))".format(key,key*2)
        return ret

    def get_case_id(self,case_dict):
        ret = 0
        for key, on in case_dict.iteritems():
            ret += on << (key*2)
        return ret

    # check if a return value tensor is zero for perfomance optimizations
    def is_zero(self, out_spin, out_tensor):
        if out_spin == 1:
            return out_tensor._array[0]==0.0
        return all([(out_tensor._array[i]._array[0]==0.0) for i in range(len(out_tensor._array))])
