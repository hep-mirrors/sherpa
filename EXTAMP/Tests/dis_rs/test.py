#!/usr/bin/env python2
from mpi4py import MPI
from json import load, dump

def writepoint(moms, me, fname):
    moms = [[mom[0],mom[1],mom[2],mom[3]] for mom in moms]
    with open(fname, 'w') as outfile:
        dump([moms, me], outfile)

if __name__ == "__main__":

    from argparse import ArgumentParser
    arg_parser = ArgumentParser()
    arg_parser.add_argument('--sherpa',   help='Path to Sherpa installation',   default='../../../')
    arg_parser.add_argument('--args',     help='Arguments passed on to Sherpa', default='')
    arg_parser.add_argument('--generate', help='Only generate test points',     action='store_true')
    arg_parser.add_argument('--openloops',help='Use OpenLoops',                 action='store_true')
    args        = arg_parser.parse_args()
    sherpa_args = ["Sherpa"]+(args.args.split())

    if args.openloops:
        sherpa_args.append("-fRun_ol.dat")
    else :
        sherpa_args.append("-fRun_amegic.dat")
        
    from subprocess import check_output
    pylibdir = check_output([args.sherpa+'/bin/Sherpa-config',
                             '--python-libs']).strip()


    print pylibdir
    from sys import path
    path.append(pylibdir)
    import Sherpa
    Generator=Sherpa.Sherpa()
    Generator.InitializeTheRun(len(sherpa_args), sherpa_args)
    Process=Sherpa.MEProcess(Generator)
    Process.Initialize();

    if(args.generate):
        Process.TestPoint(10.0)
        writepoint(Process.GetMomenta(), Process.CSMatrixElement(),
                   'momenta_0')
        Process.TestPoint(50.0)
        writepoint(Process.GetMomenta(), Process.CSMatrixElement(),
                   'momenta_1')
        Process.TestPoint(100.0)
        writepoint(Process.GetMomenta(), Process.CSMatrixElement(),
                   'momenta_2')
        Process.TestPoint(500.0)
        writepoint(Process.GetMomenta(), Process.CSMatrixElement(),
                   'momenta_3')

    else:
        for fname in ['momenta_0', 'momenta_1', 'momenta_2', 'momenta_3']:
            with open(fname, 'r') as infile:
                [moms, exp_me] = load(infile)
            Process.SetMomenta(moms)
            me = Process.CSMatrixElement()
            me = Process.CSMatrixElement()
            me = Process.CSMatrixElement()
            rdev = (me - exp_me)/(exp_me if exp_me!=0.0 else 1.0)
            print '---------------------'
            print 'Squared ME:          ', me
            print 'Expected:            ', exp_me
            print 'Rel. dev. from exp.: ', rdev
            print '---------------------'

            assert(abs(rdev) < 1.0e-10)
            
    
    exit(0)
    
