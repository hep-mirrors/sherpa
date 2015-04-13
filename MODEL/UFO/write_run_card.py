from ufo_interface import s_parameter
from ufo_interface.templates import run_card_template

def write_run_card(model, model_name):

    ext_params = [s_parameter(param) for param in  model.all_parameters if (s_parameter(param).is_external())]
    blocks     = dict()
    for param in ext_params:
        cur_block = param.lha_block().lower()
        if not cur_block in blocks:
            blocks[cur_block]=[par for par in ext_params if par.lha_block().lower()==cur_block]

    ufo_params = ""

    for block,param_list in blocks.iteritems():
        if (block.lower() == "decay"): continue # in order to comply with weird default ufo param_card format
        ufo_params += "block {0}\n\t".format(block)
        ufo_params += "\n\t".join([" ".join([str(ind) for ind in param.lha_indices()])+" "+str(param.raw_value()) for param in param_list])
        ufo_params += "\n"

    # in order to comply with weird default ufo param_card format
    if "decay" in blocks:
        for param in blocks["decay"]:
            ufo_params += "\ndecay {0} {1}".format(param.lha_indices()[0], str(param.raw_value()))
        
    with open("Run.dat", "w") as outfile:
        outfile.write(run_card_template.substitute(model=model, model_name=model_name, ufo_params=ufo_params))
