# Some ugly functions required for generation of code.
# Functions return strings to be dumped directly
# into the lorentz calculator source code or tensors
# with symbolic string-like expressions.

from tensor import tensor
from lorentz_structures import mink_metric
from sym_var import sym_var

type_dict = {
    1 : "CScalar",
    2 : "CSpinor",
    3 : "CVec4",
    5 : "CTensor"
}

vect_gauge_dict = {
    0 : "0",
    1 : "ATOOLS::Spinor<SType>::R1()",
    2 : "ATOOLS::Spinor<SType>::R2()",
    3 : "ATOOLS::Spinor<SType>::R3()"
}

def get_in_current_declaration(spin, key, index, ferm_partner):
    ret = ""
    if spin != 2:
        ret += "const {0} <SType> & j{1} = *(jj[{2}]->Get< {0} <SType> >());\n".format(type_dict[spin], key, index)
    else:
        ret += ("const {0} <SType> & j{1} = ((jj[{2}]->Get< {0} <SType> >())->B() == {3}) ? " +
                "(*(jj[{2}]->Get< {0} <SType> >())) : " +
                "(*(jj[{2}]->Get< {0} <SType> >())).CConj() ;\n"
            ).format(type_dict[spin], key, index, ferm_partner)
    for i in range(1 if spin==1 else 4):
        ret += "const SComplex & j{0}{1} = j{0}[{2}];\n".format(key, i, vect_gauge_dict[i] if spin==3 else i)
    return ret
        

def get_in_mom_declaration(key, index):
    ret = "const ATOOLS::Vec4D & p{0} = p_v->J({1})->P();\n".format(key, index)
    for i in range(4):
        ret += "const double& p{0}{1} = p{0}[{2}];\n".format(key, i, vect_gauge_dict[i])
    return ret

def get_out_mom_declaration(out_key, key_index_dict):
    keys = key_index_dict.keys()
    keys.remove(out_key)
    assert(len(keys)>0)
    ret =  "ATOOLS::Vec4D p{0} = -p{1}".format(out_key, keys[0])
    for key in keys[1:]:
        ret += "-p{0}".format(key)
    ret += ';\n'
    for i in range(4):
        ret += "const double& p{0}{1} = p{0}[{2}];\n".format(out_key, i, vect_gauge_dict[i])
    return ret
            
def get_out_current_declaration(out_spin, out_key):
    if (out_spin == 1):
        return "CScalar<SType>* j{0} = NULL;\n".format(out_key)
    elif (out_spin == 2):
        return "CSpinor<SType>* j{0} = NULL;\n".format(out_key)
    elif (out_spin == 3):
        return "CVec4<SType>* j{0} = NULL;\n".format(out_key)
    else:
        raise ufo_exception("Cannot handle spin {0}".format(out_spin))

# get a tensor representation of current
# with an index-key 'key'
def get_in_current_tens(key, spin):
    # scalar
    if spin == 1:
        return tensor([sym_var("j{0}0".format(key))], None)
    # fermion
    if (spin == 2):
        # put key+1 as key to account for UFO key conv.
        return tensor([tensor([sym_var("j{0}{1}".format(key,0))] , None), 
                       tensor([sym_var("j{0}{1}".format(key,1))] , None), 
                       tensor([sym_var("j{0}{1}".format(key,2))] , None),
                       tensor([sym_var("j{0}{1}".format(key,3))] , None)], key+1)


    if (spin == 3):
        dummy =  tensor([tensor([sym_var("j{0}{1}".format(key,0))] , None), 
                         tensor([sym_var("j{0}{1}".format(key,1))] , None), 
                         tensor([sym_var("j{0}{1}".format(key,2))] , None),
                         tensor([sym_var("j{0}{1}".format(key,3))] , None)], 'dummy_key')
        # Incoming currents are alway contravariant in sherpa, 
        # so we need to multiply by metric.
        # Put key+1 as key to account for UFO key conv.
        return dummy * mink_metric(key+1, 'dummy_key')

    raise ufo_exception("External wavefunction for spin {0} not implemented".format(spin))

def get_out_current_initialization(out_spin, out_key, out_tensor, out_bar_type):
    # scalar
    if (out_spin == 1):
        assert(out_tensor._toplevel_dim == 1)
        return "j{0} = CScalar<SType>::New({1});\n".format(out_key, out_tensor._array[0])
    # fermion
    elif (out_spin == 2):
        assert(out_tensor._toplevel_dim == 4)
        components = [out_tensor._array[i]._array[0] for i in range(4)]
        if (components[2]==0.0) and (components[3]==0.0):
            on_type = 1 
        elif (components[0]==0.0) and (components[1]==0.0):
            on_type = 2
        else:
            on_type = 3
        string = ""
        string += "j{0} = CSpinor<SType>::New(m_r[{0}],{1},0,0,0,0,{2});\n".format(out_key,
                                                                                       out_bar_type,
                                                                                       on_type)
        string += "(*j{0})[0] = {1};\n".format(out_key,components[0])
        string += "(*j{0})[1] = {1};\n".format(out_key,components[1])
        string += "(*j{0})[2] = {1};\n".format(out_key,components[2])
        string += "(*j{0})[3] = {1};\n".format(out_key,components[3])
        return string
    # vector
    elif (out_spin == 3):
        assert(out_tensor._toplevel_dim == 4)
        string = ""
        string += "j{0} = CVec4<SType>::New();\n".format(out_key)
        string += "(*j{0})[{1}] = {2};\n".format(out_key, vect_gauge_dict[0], out_tensor._array[0]._array[0])
        string += "(*j{0})[{1}] = {2};\n".format(out_key, vect_gauge_dict[1], out_tensor._array[1]._array[0])
        string += "(*j{0})[{1}] = {2};\n".format(out_key, vect_gauge_dict[2], out_tensor._array[2]._array[0])
        string += "(*j{0})[{1}] = {2};\n".format(out_key, vect_gauge_dict[3], out_tensor._array[3]._array[0])
        return string

    raise ufo_exception("External wavefunction for spin {0} not implemented".format(out_spin))
