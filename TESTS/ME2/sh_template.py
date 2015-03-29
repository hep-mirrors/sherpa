#!/usr/bin/env python2
## from mpi4py import MPI
import sys
from math import sqrt
sys.path.append('/mt/home/kuttimalai/work/SHERPA/lib/python2.7/site-packages')
sys.argv.append('INIT_ONLY=2')
import Sherpa

n_flav = len(${in_flavs})+len(${out_flavs})

if len(sys.argv)<2:
    raise RuntimeError('Too few arguments given')

me_fname = sys.argv.pop(1)
devs = []

Generator=Sherpa.Sherpa()
Generator.InitializeTheRun(len(sys.argv),sys.argv)
Process=Sherpa.MEProcess(Generator)

# Incoming flavors must be added first!
for fl in ${in_flavs}:
    Process.AddInFlav(fl);
for fl in ${out_flavs}:
    Process.AddOutFlav(fl);

Process.Initialize();
if Process.HasColorIntegrator(): Process.GenerateColorPoint()

with open(me_fname,'r') as f:
    for line in f:
        numbers = line.strip().split()
        if not len(numbers)-1==4*n_flav:
            raise RuntimeError("Number of momenta in in file {0} doesn't match process.".format(me_fname))

        for i in range(n_flav):
            Process.SetMomentum(i, float(numbers[0+i*4]), float(numbers[1+i*4]), float(numbers[2+i*4]), float(numbers[3+i*4]))
            
        f_me = float(numbers[-1])
        s_me = Process.CSMatrixElement()
        devs.append((f_me-s_me)/(f_me))

max_dev = max(devs)
min_dev = min(devs)
spread  = sqrt(sum([de*de for de in devs])/len(devs))
passed  = max([abs(max_dev),abs(min_dev)])<${threshold}

msg  =  "\n--------------------------------------------------------\n"
msg += "Test for ${procstring} {0}passed: \n".format('' if passed else 'not ' )
msg += "{0:45} {1: 1.3e}\n".format('Expected deviation', ${threshold})
msg += "{0:45} {1: 1.3e}\n".format('Maximum deviation', max_dev)
msg += "{0:45} {1: 1.3e}\n".format('Minimum deviation', min_dev)
msg += "{0:45} {1: 1.3e}\n".format('RMS deviation', spread)
msg += "--------------------------------------------------------\n"

if not passed:
    raise Sherpa.Exception(msg)

else:
    print msg
    exit(0)
    
