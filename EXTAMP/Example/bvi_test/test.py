#!/usr/bin/env python2
from mpi4py import MPI
import sys
import Sherpa

# Add this to the execution arguments to prevent Sherpa from starting the cross section integration
sys.argv.append('INIT_ONLY=2')

Generator=Sherpa.Sherpa()
try:
    Generator.InitializeTheRun(len(sys.argv),sys.argv)
    Process=Sherpa.MEProcess(Generator)

    Process.Initialize();
    
    # Random momenta
    E_cms = 5000.0
    if 'Amegic' not in Process.GeneratorName():
        wgt = Process.TestPoint(E_cms)

    for i in range(5):
        wgt = Process.TestPoint(E_cms)
        print 'Squared ME: ', Process.CSMatrixElement(), '\n'
        print Process.GetMomenta()[2]
    print Process.GeneratorName()

except Sherpa.SherpaException as exc:
    print exc
    exit(1)
