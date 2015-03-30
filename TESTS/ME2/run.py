#!/usr/bin/python2

from os import path,makedirs,getcwd,chdir
from shutil import copyfile,copy
from subprocess import check_call
from glob import glob

from test_processes import all_test_procs as procs
from mg_helpers import write_mg_card, call_mg
from templates import mg_template,sh_template,run_template

def write_run_card(proc):
    dir_path = path.join(proc.test_dir_name(),'Run.dat')
    proc_string  = "Process "
    proc_string += " ".join([str(ide) for ide in proc._in_ids])
    proc_string += " -> "
    proc_string += " ".join([str(ide) for ide in proc._out_ids])
    with open(dir_path,'w') as f:
        f.write(run_template.substitute(en=(proc._cms/2.0),
                                        in_flav0=proc._in_ids[0],
                                        in_flav1=proc._in_ids[1],
                                        model=proc._model._name,
                                        process_string=proc_string))

def write_param_card(proc):
    copyfile(proc._model._param_card,
             path.join(proc.test_dir_name(),'param_card.dat'))

def write_sh_script(proc):
    f_path = path.join(proc.test_dir_name(),'run_sh.py')
    with open(f_path,'w') as f:
        f.write(sh_template.substitute(in_flavs=proc._in_ids,
                                       out_flavs=proc._out_ids,
                                       threshold=proc._threshold,
                                       procstring=str(proc)))

def call_sh_script(proc):
    cur_dir = getcwd()
    chdir(proc.test_dir_name())
    check_call(['chmod','+x','./run_sh.py'])
    check_call(['./run_sh.py','mg_results.dat'])
    chdir(cur_dir)
    
def mg_prepare(proc):
    # make a directory for the process
    dir_path = path.join('./',proc.test_dir_name())
    if not path.exists(dir_path):
        makedirs(dir_path)
    # write a process card
    write_mg_card(proc,dir_path)
    # generate/copy makefile, param_card, executable src
    copyfile('Makefile',path.join(dir_path,'Makefile'))
    with open(path.join(dir_path,'run_mg.cc'),'w') as f:
        f.write(mg_template.substitute(cms=proc._cms,
                                       in_flav0=proc._in_ids[0],
                                       in_flav1=proc._in_ids[1],
                                       n_flavs=len(proc._out_ids)+len(proc._in_ids),
                                       ncalls=100))

def mg_generate(proc):
    dir_path = path.join('./',proc.test_dir_name())
    # run mg on generated process card
    call_mg(dir_path)
    # copy generated code into dir_path
    src  = glob(path.join(dir_path,'PROC_*/SubProcesses/P0*/*.cc'))
    src += glob(path.join(dir_path,'PROC_*/SubProcesses/P0*/*.h'))
    src += glob(path.join(dir_path,'PROC_*/src/*.cc'))
    src += glob(path.join(dir_path,'PROC_*/src/*.h'))
    for sr in src:
        copy(sr, dir_path)

# compile the generated mg code and the executable
def mg_compile(proc):
    cur_dir=getcwd()
    chdir(path.join('./',proc.test_dir_name()))
    check_call(['make'])
    chdir(cur_dir)

# run the mg executable, writes ps points+MEs
def mg_run(proc):
    cur_dir=getcwd()
    chdir(path.join('./',proc.test_dir_name()))
    check_call(['./run_mg'])
    chdir(cur_dir)


for proc in procs:
    # mg_prepare(proc)
    write_param_card(proc)
    # mg_generate(proc)
    # mg_compile(proc)
    # mg_run(proc)
    write_run_card(proc)
    write_sh_script(proc)
    call_sh_script(proc)

print "\n All ME tests for {0} processes passed\n".format(len(procs))
    
exit(0)
    
