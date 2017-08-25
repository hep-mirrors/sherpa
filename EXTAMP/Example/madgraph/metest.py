#!/usr/bin/env python2
from mpi4py import MPI
import sys
sys.path.append('/u/th/silvan/work/unlops/lib/python2.7/site-packages')
import Sherpa
sys.argv.append('INIT_ONLY=2')

Generator=Sherpa.Sherpa()
Generator.InitializeTheRun(len(sys.argv),sys.argv)
Process=Sherpa.MEProcess(Generator)

Process.Initialize();

Process.SetMomenta([[250,0,0,250],
                    [250,-0,-0,-250],
                    [250,245.669,-44.8756,11.5211],
                    [250,-245.669,44.8756,-11.5211]])

print 'Squared ME: ', Process.CSMatrixElement(), '\n'

