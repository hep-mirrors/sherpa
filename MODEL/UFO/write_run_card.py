from ufo_interface import s_parameter, s_particle
from ufo_interface.templates import run_card_template

# For all parameters except the "decay" parameters
def table_format(nci, lha_indices, ncv, value, name):
    formatter = r"{0: <"+str(nci)+r"}"
    ret  = "\t"
    ret += " ".join([formatter.format(index) for index in lha_indices])
    ret += (r" {0: <"+str(ncv)+r"}").format(value)
    ret += " # {0}\n".format(name)
    return ret

def write_run_card(model, model_name, run_card_path):

    ext_params = [s_parameter(param) for param in  model.all_parameters if (s_parameter(param).is_external())]
    
    # length of the 'index'-fields in output
    nci        = max([max([len(str(index)) for index in param.lha_indices()]) for param in ext_params])
    # length of the 'value'-fields in output
    ncv        = max([len(str(param.raw_value())) for param in ext_params])
    
    blocks     = dict()
    for param in ext_params:
        cur_block = param.lha_block().lower()
        if not cur_block in blocks:
            blocks[cur_block]=[par for par in ext_params if par.lha_block().lower()==cur_block]

    yaml_indent = "  "

    ufo_params = ""

    for block,param_list in blocks.iteritems():
        if (block.lower() == "decay"): continue # in order to comply with weird default ufo param_card format
        ufo_params += yaml_indent + "block {0}\n".format(block)
        ufo_params += "".join([yaml_indent + table_format(nci,param.lha_indices(),
                                                          ncv, param.raw_value(),
                                                          param.name()) for param in param_list])
        ufo_params += yaml_indent + "\n"

    # in order to comply with weird default ufo param_card format
    if "decay" in blocks:
        for param in blocks["decay"]:
            ufo_params += yaml_indent
            ufo_params += "decay "
            ufo_params += table_format(nci, param.lha_indices(), ncv, param.raw_value(), param.name())

    # generate a helpful template for a user specification of coupling orders 
    order_statement = 'Order: {' + ','.join([order.name + ': Any' for order in model.all_orders]) + '}'

    # collect all particles of the model for an example process section
    all_particles = [s_particle(p) for p in model.all_particles]
    all_particles = ",".join([str(p.kf_code()) for p in all_particles if not (p.is_goldstone() or p.is_ghost())])
        
    with open(run_card_path, "w") as outfile:
        outfile.write(run_card_template.substitute(model=model, 
                                                   model_name=model_name, 
                                                   ufo_params=ufo_params,
                                                   order_statement=order_statement,
                                                   all_particles=all_particles))
