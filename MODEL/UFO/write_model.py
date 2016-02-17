from ufo_interface import s_vertex, s_parameter, s_particle, s_coupling, colour_translate, split_by_orders, vertex_collection
from ufo_interface.templates import model_template
from operator import attrgetter
from copy import deepcopy

def sort_by_hierarchy(orders):
    ret = deepcopy(orders)
    ret.sort(key=attrgetter('hierarchy'))
    return ret

def write_model(model, lorentzes, model_name, model_file_name):

    para_init = ""
    part_init = ""

    external_parameters = [s_parameter(param) for param in  model.all_parameters if (s_parameter(param).is_external())]
    internal_parameters = [s_parameter(param) for param in  model.all_parameters if (s_parameter(param).is_internal())]

    # external parameter initialization
    for param in external_parameters:
        statement = "    double "+param.name()+" = p_dataread->GetEntry<double>(\""+str(param.lha_block())+"\", "
        if len(param.lha_indices()) == 1:
            statement += str(param.lha_indices()[0])
        if len(param.lha_indices()) == 2:
            statement += str(param.lha_indices()[0])+", "+str(param.lha_indices()[1])
        statement+=");\n"
        statement+='    p_constants->insert(make_pair(string("{0}"),{0}));'.format(param.name())
        para_init += "\n"+statement

    # internal parameter initialization and calculation
    for param in internal_parameters:
        if param.is_complex():
            statement = "    Complex "+param.name()+" = "+param.cpp_value()+";"
        else:
            statement = "    double "+param.name()+" = ToDouble("+param.cpp_value()+");"
        para_init += "\n"+statement
        para_init += "\n    DEBUG_VAR({0});".format(param.name())

    # fill particle list
    for s_part in [ s_particle(p) for p in model.all_particles ]:
        kfcode = s_part.kf_code()
        # don't explicitly need to add antiparticles
        if kfcode < 0:
            continue
        massive = 0 if (s_part.ufo_particle.mass is model.parameters.ZERO) else 1
        part_init += ("\n    ATOOLS::s_kftable["+str(s_part.kf_code())+"] = new ATOOLS::Particle_Info("+
                      str(s_part.kf_code())+", "+                                 # kf_code
                      str(1000.0)+", "+                                           # mass
                      str(0.0)+", "+                                              # width
                      str(s_part.charge_times_three())+", "+                      # 3*(electrical_charge)
                      str(s_part.color())+", "+                                   # strong charge
                      str(s_part.spin_times_two())+", "+                          # 2*spin
                      str(s_part.self_conjugate())+", "+                          # self_conjugate
                      str(1)+", "+                                                # is active
                      str(0)+", "+                                                # stable
                      str(massive)+", "+                                          # massive
                      "\""+str(s_part.name())+"\", "+                             # name
                      "\""+str(s_part.antiname())+"\", "+                         # antiname
                      "\""+str(s_part.texname())+"\", "+                          # texname
                      "\""+str(s_part.antitexname())+"\");")                      # antitexname

        wstring = s_part.width().name() if s_part.width().is_external() else s_part.width().cpp_value()
        mstring = s_part.mass().name()  if s_part.mass().is_external()  else s_part.mass().cpp_value()
        para_init += "\n    ATOOLS::Flavour({0}).SetWidth(ToDouble({1}));".format(kfcode,wstring)
        para_init += "\n    ATOOLS::Flavour({0}).SetMass(ToDouble({1}));".format(kfcode,mstring)
        para_init += "\n    ATOOLS::Flavour({0}).SetHadMass(ToDouble({1}));".format(kfcode,mstring)

    # coupling initialization and calculation
    for coup in model.all_couplings:
        s_coup = s_coupling(coup)
        para_init += "\n    p_complexconstants->insert(make_pair(string(\""+s_coup.name()+"\"),"+s_coup.cpp_value()+"));"
        para_init += "\n    DEBUG_VAR((*p_complexconstants)[\"{0}\"]);".format(s_coup.name())

    sorted_orders = sort_by_hierarchy(model.all_orders)
    hierarchy     = [order.name for order in sorted_orders]
    declarations  = ""
    calls         = ""
    vertices      = list(sum([split_by_orders(vert, hierarchy) for vert in model.all_vertices],[]))
    vertices      = [vert for vert in vertices if not (vert.has_ghosts() or vert.has_goldstones()) ]
    i             = 0

    # initialization of vertices
    while (len(vertices) > 0):
        sub_vert_list = []
        id_string = "vertices_{0}".format(i)
        while ((len(sub_vert_list)<10) and len(vertices) > 0):
            sub_vert_list.append(vertices.pop())

        v_collection  = vertex_collection(sub_vert_list, id_string)
        declarations += v_collection.implementation_string()
        calls        += v_collection.call_string()
        i            += 1

    # fill the lorentz name map
    lorentz_map = ""
    for lor in lorentzes:
        name, mapped_name = lor.name(), lor.mapped_name()
        if  mapped_name is not None:
            lorentz_map +='\n      m_lorentz_map.insert(make_pair("{0}","{1}"));'.format(name,mapped_name)

    # write out model
    with open(model_file_name, "w") as outfile:
        outfile.write(model_template.substitute(model_name=model_name, particle_init=part_init, param_init=para_init,
                                                declarations = declarations, calls = calls, fill_lorentz_map=lorentz_map))
