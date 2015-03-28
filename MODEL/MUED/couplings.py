# This file was automatically created by FeynRules $Revision: 915 $
# Mathematica version: 7.0 for Mac OS X x86 (64-bit) (November 11, 2008)
# Date: Mon 23 Jan 2012 11:06:17


from object_library import all_couplings, Coupling

from function_library import complexconjugate, re, im, csc, sec, acsc, asec



GC_1 = Coupling(name = 'GC_1',
                value = '-G',
                order = {'QCD':1})

GC_2 = Coupling(name = 'GC_2',
                value = 'complex(0,1)*G',
                order = {'QCD':1})

GC_3 = Coupling(name = 'GC_3',
                value = 'G',
                order = {'QCD':1})

GC_4 = Coupling(name = 'GC_4',
                value = 'complex(0,1)*G**2',
                order = {'QCD':2})

GC_5 = Coupling(name = 'GC_5',
                value = '(3.0*complex(0,1)*G**2)/2.',
                order = {'QCD':2})

GC_6 = Coupling(name = 'GC_6',
                value = '-(complex(0,1)*gw)/2.',
                order = {'QED':1})

GC_7 = Coupling(name = 'GC_7',
                value = '(complex(0,1)*gw)/2.',
                order = {'QED':1})

GC_8 = Coupling(name = 'GC_8',
                value = 'complex(0,1)*gw',
                order = {'QED':1})

GC_9 = Coupling(name = 'GC_9',
                value = '(complex(0,1)*gw)/cmath.sqrt(2)',
                order = {'QED':1})

GC_10 = Coupling(name = 'GC_10',
                 value = '(CKM11*complex(0,1)*gw)/cmath.sqrt(2)',
                 order = {'QED':1})

GC_11 = Coupling(name = 'GC_11',
                 value = '(CKM12*complex(0,1)*gw)/cmath.sqrt(2)',
                 order = {'QED':1})

GC_12 = Coupling(name = 'GC_12',
                 value = '(CKM13*complex(0,1)*gw)/cmath.sqrt(2)',
                 order = {'QED':1})

GC_13 = Coupling(name = 'GC_13',
                 value = '(CKM21*complex(0,1)*gw)/cmath.sqrt(2)',
                 order = {'QED':1})

GC_14 = Coupling(name = 'GC_14',
                 value = '(CKM22*complex(0,1)*gw)/cmath.sqrt(2)',
                 order = {'QED':1})

GC_15 = Coupling(name = 'GC_15',
                 value = '(CKM23*complex(0,1)*gw)/cmath.sqrt(2)',
                 order = {'QED':1})

GC_16 = Coupling(name = 'GC_16',
                 value = '(CKM31*complex(0,1)*gw)/cmath.sqrt(2)',
                 order = {'QED':1})

GC_17 = Coupling(name = 'GC_17',
                 value = '(CKM32*complex(0,1)*gw)/cmath.sqrt(2)',
                 order = {'QED':1})

GC_18 = Coupling(name = 'GC_18',
                 value = '(CKM33*complex(0,1)*gw)/cmath.sqrt(2)',
                 order = {'QED':1})

GC_19 = Coupling(name = 'GC_19',
                 value = 'cw*complex(0,1)*gw',
                 order = {'QED':1})

GC_20 = Coupling(name = 'GC_20',
                 value = '-(complex(0,1)*gw**2)',
                 order = {'QED':2})

GC_21 = Coupling(name = 'GC_21',
                 value = 'complex(0,1)*gw**2',
                 order = {'QED':2})

GC_22 = Coupling(name = 'GC_22',
                 value = '(-3.0*complex(0,1)*gw**2)/2.',
                 order = {'QED':2})

GC_23 = Coupling(name = 'GC_23',
                 value = '(3.0*complex(0,1)*gw**2)/2.',
                 order = {'QED':2})

GC_24 = Coupling(name = 'GC_24',
                 value = 'cw*complex(0,1)*gw**2',
                 order = {'QED':2})

GC_25 = Coupling(name = 'GC_25',
                 value = 'cw**2*complex(0,1)*gw**2',
                 order = {'QED':2})

GC_26 = Coupling(name = 'GC_26',
                 value = '-6.0*complex(0,1)*lam',
                 order = {'QED':2})

GC_27 = Coupling(name = 'GC_27',
                 value = '(ee**2*complex(0,1))/(2.*sw**2)',
                 order = {'QED':2})

GC_28 = Coupling(name = 'GC_28',
                 value = 'complex(0,1)*gw*sw',
                 order = {'QED':1})

GC_29 = Coupling(name = 'GC_29',
                 value = '-2.0*complex(0,1)*gw**2*sw',
                 order = {'QED':2})

GC_30 = Coupling(name = 'GC_30',
                 value = '-2.0*cw*complex(0,1)*gw**2*sw',
                 order = {'QED':2})

GC_31 = Coupling(name = 'GC_31',
                 value = 'complex(0,1)*gw**2*sw**2',
                 order = {'QED':2})

GC_32 = Coupling(name = 'GC_32',
                 value = 'ee**2*complex(0,1) + (cw**2*ee**2*complex(0,1))/(2.*sw**2) + (ee**2*complex(0,1)*sw**2)/(2.*cw**2)',
                 order = {'QED':2})

GC_33 = Coupling(name = 'GC_33',
                 value = '-6.0*complex(0,1)*lam*v',
                 order = {'QED':1})

GC_34 = Coupling(name = 'GC_34',
                 value = '(ee**2*complex(0,1)*v)/(2.*sw**2)',
                 order = {'QED':1})

GC_35 = Coupling(name = 'GC_35',
                 value = 'ee**2*complex(0,1)*v + (cw**2*ee**2*complex(0,1)*v)/(2.*sw**2) + (ee**2*complex(0,1)*sw**2*v)/(2.*cw**2)',
                 order = {'QED':1})

GC_36 = Coupling(name = 'GC_36',
                 value = '(complex(0,1)*g1*y1)/2.',
                 order = {'QED':1})

GC_37 = Coupling(name = 'GC_37',
                 value = '-(complex(0,1)*gw*sw)/2. + (cw*complex(0,1)*g1*y1)/2.',
                 order = {'QED':1})

GC_38 = Coupling(name = 'GC_38',
                 value = '(complex(0,1)*gw*sw)/2. + (cw*complex(0,1)*g1*y1)/2.',
                 order = {'QED':1})

GC_39 = Coupling(name = 'GC_39',
                 value = '-(cw*complex(0,1)*gw)/2. - (complex(0,1)*g1*sw*y1)/2.',
                 order = {'QED':1})

GC_40 = Coupling(name = 'GC_40',
                 value = '(cw*complex(0,1)*gw)/2. - (complex(0,1)*g1*sw*y1)/2.',
                 order = {'QED':1})

GC_41 = Coupling(name = 'GC_41',
                 value = '(complex(0,1)*g1*y2)/2.',
                 order = {'QED':1})

GC_42 = Coupling(name = 'GC_42',
                 value = '(cw*complex(0,1)*g1*y2)/2.',
                 order = {'QED':1})

GC_43 = Coupling(name = 'GC_43',
                 value = '-(complex(0,1)*g1*sw*y2)/2.',
                 order = {'QED':1})

GC_44 = Coupling(name = 'GC_44',
                 value = '(complex(0,1)*g1*y3)/2.',
                 order = {'QED':1})

GC_45 = Coupling(name = 'GC_45',
                 value = '-(complex(0,1)*gw*sw)/2. + (cw*complex(0,1)*g1*y3)/2.',
                 order = {'QED':1})

GC_46 = Coupling(name = 'GC_46',
                 value = '(complex(0,1)*gw*sw)/2. + (cw*complex(0,1)*g1*y3)/2.',
                 order = {'QED':1})

GC_47 = Coupling(name = 'GC_47',
                 value = '-(cw*complex(0,1)*gw)/2. - (complex(0,1)*g1*sw*y3)/2.',
                 order = {'QED':1})

GC_48 = Coupling(name = 'GC_48',
                 value = '(cw*complex(0,1)*gw)/2. - (complex(0,1)*g1*sw*y3)/2.',
                 order = {'QED':1})

GC_49 = Coupling(name = 'GC_49',
                 value = '(complex(0,1)*g1*y4)/2.',
                 order = {'QED':1})

GC_50 = Coupling(name = 'GC_50',
                 value = '(cw*complex(0,1)*g1*y4)/2.',
                 order = {'QED':1})

GC_51 = Coupling(name = 'GC_51',
                 value = '-(complex(0,1)*g1*sw*y4)/2.',
                 order = {'QED':1})

GC_52 = Coupling(name = 'GC_52',
                 value = '(complex(0,1)*g1*y5)/2.',
                 order = {'QED':1})

GC_53 = Coupling(name = 'GC_53',
                 value = '(cw*complex(0,1)*g1*y5)/2.',
                 order = {'QED':1})

GC_54 = Coupling(name = 'GC_54',
                 value = '-(complex(0,1)*g1*sw*y5)/2.',
                 order = {'QED':1})

GC_55 = Coupling(name = 'GC_55',
                 value = '-((complex(0,1)*yb)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_56 = Coupling(name = 'GC_56',
                 value = '-((complex(0,1)*yc)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_57 = Coupling(name = 'GC_57',
                 value = '-((complex(0,1)*yt)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_58 = Coupling(name = 'GC_58',
                 value = '-((complex(0,1)*ytau)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_59 = Coupling(name = 'GC_59',
                 value = '(complex(0,1)*gw*complexconjugate(CKM11))/cmath.sqrt(2)',
                 order = {'QED':1})

GC_60 = Coupling(name = 'GC_60',
                 value = '(complex(0,1)*gw*complexconjugate(CKM12))/cmath.sqrt(2)',
                 order = {'QED':1})

GC_61 = Coupling(name = 'GC_61',
                 value = '(complex(0,1)*gw*complexconjugate(CKM13))/cmath.sqrt(2)',
                 order = {'QED':1})

GC_62 = Coupling(name = 'GC_62',
                 value = '(complex(0,1)*gw*complexconjugate(CKM21))/cmath.sqrt(2)',
                 order = {'QED':1})

GC_63 = Coupling(name = 'GC_63',
                 value = '(complex(0,1)*gw*complexconjugate(CKM22))/cmath.sqrt(2)',
                 order = {'QED':1})

GC_64 = Coupling(name = 'GC_64',
                 value = '(complex(0,1)*gw*complexconjugate(CKM23))/cmath.sqrt(2)',
                 order = {'QED':1})

GC_65 = Coupling(name = 'GC_65',
                 value = '(complex(0,1)*gw*complexconjugate(CKM31))/cmath.sqrt(2)',
                 order = {'QED':1})

GC_66 = Coupling(name = 'GC_66',
                 value = '(complex(0,1)*gw*complexconjugate(CKM32))/cmath.sqrt(2)',
                 order = {'QED':1})

GC_67 = Coupling(name = 'GC_67',
                 value = '(complex(0,1)*gw*complexconjugate(CKM33))/cmath.sqrt(2)',
                 order = {'QED':1})

