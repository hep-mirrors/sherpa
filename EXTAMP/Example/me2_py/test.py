#!/usr/bin/env python
from mpi4py import MPI
import sys
import Sherpa

megen = None
if "Comix" in sys.argv[-1]:
    megen = "Comix"
if "External" in sys.argv[-1]:
    megen = "External"

assert(megen is not None)

Generator=Sherpa.Sherpa()

print(sys.argv[-1])

Generator.InitializeTheRun(len(sys.argv),sys.argv)
Process=Sherpa.MEProcess(Generator)

Process.AddInFlav(21);
Process.AddInFlav(21);
Process.AddOutFlav(6);
Process.AddOutFlav(-6);
Process.AddOutFlav(21);
Process.Initialize();

Process.SetMomentum(0,250,0,0,+250)
Process.SetMomentum(1,250,0,0,-250)
Process.SetMomentum(2,203.496,106.782,-0.198869,-2.55226)
Process.SetMomentum(3,182.178,-5.85637,-12.1869,54.8116)
Process.SetMomentum(4,114.326,-100.926,12.3858,-52.2593)


from json import dump
with open(megen, 'w') as outfile:
    dump(Process.CSMatrixElement(), outfile)
        
print('Squared ME: ', Process.CSMatrixElement(), '\n')

exit(0)

