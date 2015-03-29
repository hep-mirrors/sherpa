from os import path, chdir, getcwd
from subprocess import check_call

mg_exec='/mt/home/kuttimalai/work/trunk/validation/MG5_aMC_v2_2_1/bin/mg5_aMC'

# Guess particle name in MG convention.
# First search the hard-wired dict,
# use ufo particle name if not in dict
def mg_name(pdgid, model):
    pdg_to_mg = sm_dict if model._name!='mssm' else mssm_dict
    if pdgid in pdg_to_mg:
        return pdg_to_mg[pdgid]
    for part in model._ufo_model.all_particles:
        if part.pdg_code == pdgid:
            return part.name
    print model._name
    raise RuntimeError("Cannot guess MG name for particle with pdg id {0}".format(pdgid))

def write_mg_card(test_proc, dir_path):
    with open(path.join(dir_path,'proc_card.dat'), 'w') as f:
        f.write('import model {0}\n'.format(test_proc._model._name))
        proc_string  = ' '.join([mg_name(id, test_proc._model) for id in test_proc._in_ids])
        proc_string += ' > '+' '.join([mg_name(id,test_proc._model) for id in test_proc._out_ids])
        f.write('add process '+proc_string+' QED=99 QCD=99 NP=99 HIW=99 HIG=99\n')
        f.write('output standalone_cpp')

def call_mg(dir_path):
    cur_path = getcwd()
    chdir(dir_path)
    check_call([mg_exec, 'proc_card.dat'])
    chdir(cur_path)
    

sm_dict = {
    # quarks
    1: "d",
    2: "u",
    3: "s",
    4: "c",
    5: "b",
    6: "t",
    -2: "u~",
    -1: "d~",
    -3: "s~",
    -4: "c~",
    -5: "b~",
    -6: "t~",
    # leptons
    11: "e-",
    -11: "e+",
    13: "mu-",
    -13: "mu+",
    15: "ta-",
    -15: "ta+",
    12: "ve",
    -12: "ve~",
    14: "vm", 
    -14: "vm~",
    16: "vt",
    -16: "vt~",
    # gauge bosons
    22: "a",
    21: "g",
    23: "z",
    -24: "w-",
    24: "w+",
    # higgses
    25: "h"
}

mssm_dict = {
    # quarks
    1: "d",
    2: "u",
    3: "s",
    4: "c",
    5: "b",
    6: "t",
    -2: "u~",
    -1: "d~",
    -3: "s~",
    -4: "c~",
    -5: "b~",
    -6: "t~",
    # leptons
    11: "e-",
    -11: "e+",
    13: "mu-",
    -13: "mu+",
    15: "ta-",
    -15: "ta+",
    12: "ve",
    -12: "ve~",
    14: "vm", 
    -14: "vm~",
    16: "vt",
    -16: "vt~",
    # gauge bosons
    22: "a",
    21: "g",
    23: "z",
    -24: "w-",
    24: "w+",
    # higgses
    25: "h1",
    35: "h2",
    36: "h3",
    37: "h+",
    -37: "h-",
    # squarks
    1000001: "dl",
    1000002: "ul",
    1000003: "sl",
    1000004: "cl",
    1000005: "b1",
    1000006: "t1",
    -1000001: "dl~",
    -1000002: "ul~",
    -1000003: "sl~",
    -1000004: "cl~",
    -1000005: "b1~",
    -1000006: "t1~",
    2000001: "dr",
    2000002: "ur",
    2000003: "sr",
    2000004: "cr",
    2000005: "b2",
    2000006: "t2",
    -2000001: "dr~",
    -2000002: "ur~",
    -2000003: "sr~",
    -2000004: "cr~",
    -2000005: "b2~",
    -2000006: "t2~",
    # sleptons
    1000011: "el-", 
    1000012: "sve",
    1000013: "mul-",
    1000014: "svm",
    1000015: "ta1-",
    1000016: "svt",
    2000011: "er-",
    2000013: "mur-",
    2000015: "ta2-",
    -1000011: "el+",
    -1000012: "sve~",
    -1000013: "mul+",
    -1000014: "svm~",
    -1000015: "ta1+",
    -1000016: "svt~",
    -2000011: "er+",
    -2000013: "mur+",
    -2000015: "ta2+",
    # MSSM ferms
    1000021: "go",
    1000022: "n1",
    1000023: "n2",
    1000025: "n3",
    1000035: "n4",
    1000024: "x1+",
    1000037: "x2+",
    -1000024: "x1-",
    -1000037: "x2-"
}
