from templates import color_calc_template
from color_structures import T, replacer_T, f, Identity, IdentityG
from tensor import color_key

class s_color(object):

    def __init__(self,ufo_color):

        self.ufo_color = ufo_color
        self.unique_id = self.ufo_color.replace('(','_').replace(')','_').replace('-','m').replace(',','_').replace('*','').rstrip('_')
        self._tens = None

    def name(self):
        return self.unique_id

    def c_name(self):
        return self.unique_id+'.C'

    def original_keys(self):
        """Get the original keys of the tensor representation
        before replacing gluon keys by pairs of quark keys"""
        # This attribute is set in the get_cpl_tensor method
        if not hasattr(self, "_original_keys"):
            self.get_cpl_tensor()
        return self._original_keys

    def get_cpl_tensor(self):
        """Get a tensor representation of the color coupling structure."""

        # Instead of setting up a class-wide cache, just store the
        # tensor as a member variable. Tensor can become large and
        # memory intensive, so we want to get rid of it along with its
        # s_color instance as soon as it gets destroyed.

        if self._tens is None:
            self._tens = eval(self.ufo_color)
            self._original_keys = self._tens.keys()
            # In order to translate to color-flow representation: For
            # each gluon-key, multiply with T of matching key. Thereby
            # swap gluon key for two fundamental keys 
            gluon_keys = [key for key,dim in
                          self._tens.key_dim_dict().iteritems() if dim==8]

            for gk in gluon_keys:
                # 'gk' will be replaced by str(gk)+'0' and str(gk)+'1',
                # i.e. gluon with key 2 replaced by keys "20" and "21"
                self._tens = self._tens*replacer_T(gk,
                                                   color_key(str(gk)+'0', 'fu', mapped=gk),
                                                   color_key(str(gk)+'1', 'af', mapped=gk))
        return self._tens

    def entries_string(self):
        return self.print_entries('m_cfacs',self.get_cpl_tensor())

    def print_entries(self,tag, tens):
        """Generate code lines 'm_cfacs[i0][i1][i2]...[ix]=value' that set
        the color factors of an array representing the color structure"""

        # End of recursion: just append value to tag or skip if zero
        if tens._elementary:
            value = complex(tens._array[0])
            # If the color factor is zero, the assigment can be
            # omitted as all values are assigned to zero initially
            if (abs(value.real)<1.0e-12) and (abs(value.imag)<1.0e-12):
                return ''
            value = 'std::complex<double>({0},{1})'.format(value.real, value.imag)
            return tag + '=' + value + ';\n'

        # Dereference toplevel and continue recursively
        ret = ''
        for i in range(tens._toplevel_dim):
            ret += self.print_entries(tag+'[{0}]'.format(i),
                                      tens[{tens._toplevel_key:i}])
        return ret

    def array_declaration_string(self):
        """Generate a string that initialises the color factor array
        representing the color structure"""

        kdd = self.get_cpl_tensor().key_dim_dict()
        ikd = self.get_cpl_tensor().index_key_dict()
        assert(len(ikd)==len(kdd))
        ret = ''.join(['[{0}]'.format(kdd[ikd[i]]) for i in range(len(ikd))])
        return ret

    def array_init_string(self):
        kdd = self.get_cpl_tensor().key_dim_dict()
        ikd = self.get_cpl_tensor().index_key_dict()
        assert(len(ikd)==len(kdd))
        ret = "{"*len(ikd)+" std::complex<double>(0.0,0.0) "+"}"*len(ikd)
        return ret

    def cfs_string(self):
        ret       = ''
        out_keys = self.original_keys()
        kdd      = self.get_cpl_tensor().key_dim_dict()
        ikd      = self.get_cpl_tensor().index_key_dict()
        for k in out_keys:
            ret += "\ncase {0}:\n".format(k.key-1)

            # Find all keys in the actual coupling tensor that belong
            # to this original key (can be two if the original key
            # belongs to gluon)
            mapped_keys = [k0 for k0 in self.get_cpl_tensor().keys() if k0.mapped_key==k]
            assert(len(mapped_keys) in [1,2])

            # If outgoing quark or antiquark
            if len(mapped_keys)==1:
                mk = mapped_keys[0]
                ret += "for(size_t i(0); i<{0}; i++)\n".format(kdd[mk])
                aid = ""
                for i in range(len(ikd)):
                    if ikd[i] == mk:
                        aid += '[i]'
                    else:
                        aid += '[(*j[m_inds[{0}]])({1})-1]'.format(ikd[i].mapped_key.key-1,
                                                                   0 if ikd[i].rep=='fu' else 1)
                ret += 'if(m_cfacs{2}!=0.0)\nm_c.push_back(CInfo({0},{1},m_cfacs{2}));\n'.format("i+1"if mk.rep=='af' else "0",
                                                                                                 "i+1"if mk.rep=='fu' else "0", aid)

            # If outgoing gluon
            elif len(mapped_keys)==2:
                mk0 = mapped_keys[0]
                mk1 = mapped_keys[1]
                ret += "for(size_t i(0); i<{0}; i++)\n".format(kdd[mk0])
                ret += "for(size_t m(0); m<{0}; m++)\n".format(kdd[mk1])
                aid = ""
                for i in range(len(ikd)):
                    if ikd[i] == mk0:
                        aid += '[i]'
                    elif ikd[i] == mk1:
                        aid += '[m]'
                    else:
                        aid += '[(*j[m_inds[{0}]])({1})-1]'.format(ikd[i].mapped_key.key-1,
                                                                   0 if ikd[i].rep=='fu' else 1)
                assert(mk0.rep=='af' or mk1.rep=='af')
                assert(mk0.rep=='fu' or mk1.rep=='fu')
                ret += 'if(m_cfacs{2}!=0.0)\nm_c.push_back(CInfo({0},{1},m_cfacs{2}));\n'.format("i+1" if mk0.rep=='af' else "m+1",
                                                                                                 "m+1" if mk1.rep=='fu' else "i+1", aid)
            ret +='break;\n'
        ret += "\ndefault:\n"
        aid = ''.join(['[(*j[m_inds[{0}]])({1})-1]'.format(ikd[i].mapped_key.key-1,
                                                           0 if ikd[i].rep=='fu' else 1) for i in range(len(ikd))])
        ret += 'if(m_cfacs{0}!=0.0)\nm_c.push_back(CInfo(0,0,m_cfacs{0}));\n'.format(aid)

        return ret

    def write(self, path):
        with open(path+"/"+self.c_name(), "w") as outfile:
            outfile.write(color_calc_template.substitute(color_name        = self.name(),
                                                         array_declaration = self.array_declaration_string(),
                                                         array_init        = self.array_init_string(),
                                                         array_vals        = self.entries_string(),
                                                         get_cfs           = self.cfs_string()))
