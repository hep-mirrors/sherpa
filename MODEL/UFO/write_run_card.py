from ufo_interface import s_parameter
from ufo_interface.templates import run_card_template
from ufo_interface.write_model import sort_by_hierarchy

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

    ufo_params = ""

    for block,param_list in blocks.iteritems():
        if (block.lower() == "decay"): continue # in order to comply with weird default ufo param_card format
        ufo_params += "block {0}\n".format(block)
        ufo_params += "".join([table_format(nci,param.lha_indices(),
                                            ncv, param.raw_value(),
                                            param.name()) for param in param_list])
        ufo_params += "\n"

    # in order to comply with weird default ufo param_card format
    if "decay" in blocks:
        for param in blocks["decay"]:
            ufo_params += "decay "+table_format(nci, param.lha_indices(), ncv, param.raw_value(), param.name())

    sorted_orders = sort_by_hierarchy(model.all_orders)
    # generate a helpful template for a user specification of coupling orders 
    order_statement = 'Order (' + ','.join(['*' for order in sorted_orders]) + ')'
    order_comment   = '# Syntax: "Order (' + ','.join([order.name for order in sorted_orders]) + ');"'
        
    with open(run_card_path, "w") as outfile:
        outfile.write(run_card_template.substitute(model=model, 
                                                   model_name=model_name, 
                                                   ufo_params=ufo_params,
                                                   order_statement=order_statement,
                                                   order_comment=order_comment))
