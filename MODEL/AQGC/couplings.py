# This file was automatically created by FeynRules 2.0.8
# Mathematica version: 8.0 for Linux x86 (64-bit) (February 23, 2011)
# Date: Tue 11 Nov 2014 15:33:22


from object_library import all_couplings, Coupling

from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot



GC_1 = Coupling(name = 'GC_1',
                value = '-G',
                order = {'QCD':1})

GC_2 = Coupling(name = 'GC_2',
                value = 'complex(0,1)*G**2',
                order = {'QCD':2})

GC_3 = Coupling(name = 'GC_3',
                value = 'cw*complex(0,1)*gw',
                order = {'QED':1})

GC_4 = Coupling(name = 'GC_4',
                value = '-(complex(0,1)*gw**2)',
                order = {'QED':2})

GC_5 = Coupling(name = 'GC_5',
                value = 'cw**2*complex(0,1)*gw**2',
                order = {'QED':2})

GC_6 = Coupling(name = 'GC_6',
                value = 'complex(0,1)*gw*sw',
                order = {'QED':1})

GC_7 = Coupling(name = 'GC_7',
                value = '-2*cw*complex(0,1)*gw**2*sw',
                order = {'QED':2})

GC_8 = Coupling(name = 'GC_8',
                value = 'complex(0,1)*gw**2*sw**2',
                order = {'QED':2})

