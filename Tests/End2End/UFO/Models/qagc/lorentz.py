# This file was automatically created by FeynRules $Revision: 595 $
# Mathematica version: 9.0 for Mac OS X x86 (64-bit) (November 20, 2012)
# Date: Wed 5 Jun 2013 11:59:40


from object_library import all_lorentz, Lorentz

from function_library import complexconjugate, re, im, csc, sec, acsc, asec



SSS1 = Lorentz(name = 'SSS1',
               spins = [ 1, 1, 1 ],
               structure = '1')

FFS1 = Lorentz(name = 'FFS1',
               spins = [ 2, 2, 1 ],
               structure = 'Identity(2,1)')

FFV1 = Lorentz(name = 'FFV1',
               spins = [ 2, 2, 3 ],
               structure = 'Gamma(3,2,1)')

FFV2 = Lorentz(name = 'FFV2',
               spins = [ 2, 2, 3 ],
               structure = 'Gamma(3,2,-1)*ProjM(-1,1)')

FFV3 = Lorentz(name = 'FFV3',
               spins = [ 2, 2, 3 ],
               structure = 'Gamma(3,2,-1)*ProjM(-1,1) - 2*Gamma(3,2,-1)*ProjP(-1,1)')

FFV4 = Lorentz(name = 'FFV4',
               spins = [ 2, 2, 3 ],
               structure = 'Gamma(3,2,-1)*ProjM(-1,1) + 2*Gamma(3,2,-1)*ProjP(-1,1)')

FFV5 = Lorentz(name = 'FFV5',
               spins = [ 2, 2, 3 ],
               structure = 'Gamma(3,2,-1)*ProjM(-1,1) + 4*Gamma(3,2,-1)*ProjP(-1,1)')

VVS1 = Lorentz(name = 'VVS1',
               spins = [ 3, 3, 1 ],
               structure = 'Metric(1,2)')

VVV1 = Lorentz(name = 'VVV1',
               spins = [ 3, 3, 3 ],
               structure = 'P(3,1)*Metric(1,2) - P(3,2)*Metric(1,2) - P(2,1)*Metric(1,3) + P(2,3)*Metric(1,3) + P(1,2)*Metric(2,3) - P(1,3)*Metric(2,3)')

SSSS1 = Lorentz(name = 'SSSS1',
                spins = [ 1, 1, 1, 1 ],
                structure = '1')

VVSS1 = Lorentz(name = 'VVSS1',
                spins = [ 3, 3, 1, 1 ],
                structure = 'Metric(1,2)')

VVSS2 = Lorentz(name = 'VVSS2',
                spins = [ 3, 3, 1, 1 ],
                structure = '-(P(-1,2)*P(-1,4)*P(1,3)*P(2,1)) - P(-1,2)*P(-1,3)*P(1,4)*P(2,1) - P(-1,1)*P(-1,4)*P(1,2)*P(2,3) + P(-1,1)*P(-1,2)*P(1,4)*P(2,3) - P(-1,1)*P(-1,3)*P(1,2)*P(2,4) + P(-1,1)*P(-1,2)*P(1,3)*P(2,4) + P(-2,2)*P(-2,4)*P(-1,1)*P(-1,3)*Metric(1,2) + P(-2,2)*P(-2,3)*P(-1,1)*P(-1,4)*Metric(1,2)')

VVSS3 = Lorentz(name = 'VVSS3',
                spins = [ 3, 3, 1, 1 ],
                structure = 'P(-1,3)*P(-1,4)*P(1,2)*P(2,1) - P(-2,3)*P(-2,4)*P(-1,1)*P(-1,2)*Metric(1,2)')

VVVV1 = Lorentz(name = 'VVVV1',
                spins = [ 3, 3, 3, 3 ],
                structure = 'P(1,4)*P(4,1)*Metric(2,3) - P(-1,1)*P(-1,4)*Metric(1,4)*Metric(2,3)')

VVVV2 = Lorentz(name = 'VVVV2',
                spins = [ 3, 3, 3, 3 ],
                structure = 'Metric(1,4)*Metric(2,3) - Metric(1,3)*Metric(2,4)')

VVVV3 = Lorentz(name = 'VVVV3',
                spins = [ 3, 3, 3, 3 ],
                structure = 'P(3,2)*P(4,1)*Metric(1,2) + P(3,1)*P(4,2)*Metric(1,2) - P(2,1)*P(4,2)*Metric(1,3) - P(2,1)*P(3,2)*Metric(1,4) - P(1,2)*P(4,1)*Metric(2,3) + P(-1,1)*P(-1,2)*Metric(1,4)*Metric(2,3) - P(1,2)*P(3,1)*Metric(2,4) + P(-1,1)*P(-1,2)*Metric(1,3)*Metric(2,4)')

VVVV4 = Lorentz(name = 'VVVV4',
                spins = [ 3, 3, 3, 3 ],
                structure = 'P(2,4)*P(4,2)*Metric(1,3) + P(2,3)*P(3,2)*Metric(1,4) + P(1,4)*P(4,1)*Metric(2,3) - P(-1,2)*P(-1,3)*Metric(1,4)*Metric(2,3) - P(-1,1)*P(-1,4)*Metric(1,4)*Metric(2,3) + P(1,3)*P(3,1)*Metric(2,4) - P(-1,1)*P(-1,3)*Metric(1,3)*Metric(2,4) - P(-1,2)*P(-1,4)*Metric(1,3)*Metric(2,4)')

VVVV5 = Lorentz(name = 'VVVV5',
                spins = [ 3, 3, 3, 3 ],
                structure = 'P(2,4)*P(4,3)*Metric(1,3) + P(2,3)*P(3,4)*Metric(1,4) + P(1,4)*P(4,3)*Metric(2,3) - P(-1,3)*P(-1,4)*Metric(1,4)*Metric(2,3) + P(1,3)*P(3,4)*Metric(2,4) - P(-1,3)*P(-1,4)*Metric(1,3)*Metric(2,4) - P(1,4)*P(2,3)*Metric(3,4) - P(1,3)*P(2,4)*Metric(3,4)')

VVVV6 = Lorentz(name = 'VVVV6',
                spins = [ 3, 3, 3, 3 ],
                structure = 'P(3,2)*P(4,1)*Metric(1,2) + P(3,1)*P(4,2)*Metric(1,2) - P(2,1)*P(4,2)*Metric(1,3) - P(2,4)*P(4,3)*Metric(1,3) - P(2,1)*P(3,2)*Metric(1,4) - P(2,3)*P(3,4)*Metric(1,4) - P(1,2)*P(4,1)*Metric(2,3) - P(1,4)*P(4,3)*Metric(2,3) + P(-1,1)*P(-1,2)*Metric(1,4)*Metric(2,3) + P(-1,3)*P(-1,4)*Metric(1,4)*Metric(2,3) - P(1,2)*P(3,1)*Metric(2,4) - P(1,3)*P(3,4)*Metric(2,4) + P(-1,1)*P(-1,2)*Metric(1,3)*Metric(2,4) + P(-1,3)*P(-1,4)*Metric(1,3)*Metric(2,4) + P(1,4)*P(2,3)*Metric(3,4) + P(1,3)*P(2,4)*Metric(3,4)')

VVVV7 = Lorentz(name = 'VVVV7',
                spins = [ 3, 3, 3, 3 ],
                structure = 'Metric(1,4)*Metric(2,3) + Metric(1,3)*Metric(2,4) - 2*Metric(1,2)*Metric(3,4)')

VVVV8 = Lorentz(name = 'VVVV8',
                spins = [ 3, 3, 3, 3 ],
                structure = 'Metric(1,4)*Metric(2,3) - Metric(1,2)*Metric(3,4)')

VVVV9 = Lorentz(name = 'VVVV9',
                spins = [ 3, 3, 3, 3 ],
                structure = 'Metric(1,3)*Metric(2,4) - Metric(1,2)*Metric(3,4)')

VVVV10 = Lorentz(name = 'VVVV10',
                 spins = [ 3, 3, 3, 3 ],
                 structure = 'Metric(1,4)*Metric(2,3) - (Metric(1,3)*Metric(2,4))/2. - (Metric(1,2)*Metric(3,4))/2.')

VVVV11 = Lorentz(name = 'VVVV11',
                 spins = [ 3, 3, 3, 3 ],
                 structure = 'P(1,2)*P(2,1)*Metric(3,4) - P(-1,1)*P(-1,2)*Metric(1,2)*Metric(3,4)')

VVVV12 = Lorentz(name = 'VVVV12',
                 spins = [ 3, 3, 3, 3 ],
                 structure = 'P(1,4)*P(4,1)*Metric(2,3) - P(-1,1)*P(-1,4)*Metric(1,4)*Metric(2,3) + P(1,3)*P(3,1)*Metric(2,4) - P(-1,1)*P(-1,3)*Metric(1,3)*Metric(2,4) + P(1,2)*P(2,1)*Metric(3,4) - P(-1,1)*P(-1,2)*Metric(1,2)*Metric(3,4)')

VVVV13 = Lorentz(name = 'VVVV13',
                 spins = [ 3, 3, 3, 3 ],
                 structure = 'P(3,4)*P(4,1)*Metric(1,2) + P(2,4)*P(4,1)*Metric(1,3) - P(2,4)*P(3,1)*Metric(1,4) - P(2,1)*P(3,4)*Metric(1,4) + P(1,4)*P(3,1)*Metric(2,4) - P(-1,1)*P(-1,4)*Metric(1,3)*Metric(2,4) + P(1,4)*P(2,1)*Metric(3,4) - P(-1,1)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVV14 = Lorentz(name = 'VVVV14',
                 spins = [ 3, 3, 3, 3 ],
                 structure = 'P(3,2)*P(4,1)*Metric(1,2) - P(3,4)*P(4,1)*Metric(1,2) + P(3,1)*P(4,2)*Metric(1,2) - P(3,1)*P(4,3)*Metric(1,2) + P(2,3)*P(4,1)*Metric(1,3) - P(2,4)*P(4,1)*Metric(1,3) - P(2,1)*P(4,2)*Metric(1,3) + P(2,1)*P(4,3)*Metric(1,3) - P(2,3)*P(3,1)*Metric(1,4) + P(2,4)*P(3,1)*Metric(1,4) - P(2,1)*P(3,2)*Metric(1,4) + P(2,1)*P(3,4)*Metric(1,4) - P(1,2)*P(4,1)*Metric(2,3) - P(1,3)*P(4,1)*Metric(2,3) + P(-1,1)*P(-1,2)*Metric(1,4)*Metric(2,3) + P(-1,1)*P(-1,3)*Metric(1,4)*Metric(2,3) - P(1,2)*P(3,1)*Metric(2,4) - P(1,4)*P(3,1)*Metric(2,4) + P(-1,1)*P(-1,2)*Metric(1,3)*Metric(2,4) + P(-1,1)*P(-1,4)*Metric(1,3)*Metric(2,4) - P(1,3)*P(2,1)*Metric(3,4) - P(1,4)*P(2,1)*Metric(3,4) + P(-1,1)*P(-1,3)*Metric(1,2)*Metric(3,4) + P(-1,1)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVV15 = Lorentz(name = 'VVVV15',
                 spins = [ 3, 3, 3, 3 ],
                 structure = 'P(3,4)*P(4,1)*Metric(1,2) + P(3,4)*P(4,2)*Metric(1,2) + P(3,1)*P(4,3)*Metric(1,2) + P(3,2)*P(4,3)*Metric(1,2) - P(2,3)*P(4,1)*Metric(1,3) + P(2,4)*P(4,1)*Metric(1,3) + P(2,3)*P(4,2)*Metric(1,3) - P(2,1)*P(4,3)*Metric(1,3) + P(2,3)*P(3,1)*Metric(1,4) - P(2,4)*P(3,1)*Metric(1,4) + P(2,4)*P(3,2)*Metric(1,4) - P(2,1)*P(3,4)*Metric(1,4) + P(1,3)*P(4,1)*Metric(2,3) - P(1,3)*P(4,2)*Metric(2,3) + P(1,4)*P(4,2)*Metric(2,3) - P(1,2)*P(4,3)*Metric(2,3) - P(-1,1)*P(-1,3)*Metric(1,4)*Metric(2,3) - P(-1,2)*P(-1,4)*Metric(1,4)*Metric(2,3) + P(1,4)*P(3,1)*Metric(2,4) + P(1,3)*P(3,2)*Metric(2,4) - P(1,4)*P(3,2)*Metric(2,4) - P(1,2)*P(3,4)*Metric(2,4) - P(-1,2)*P(-1,3)*Metric(1,3)*Metric(2,4) - P(-1,1)*P(-1,4)*Metric(1,3)*Metric(2,4) + P(1,3)*P(2,1)*Metric(3,4) + P(1,4)*P(2,1)*Metric(3,4) + P(1,2)*P(2,3)*Metric(3,4) + P(1,2)*P(2,4)*Metric(3,4) - P(-1,1)*P(-1,3)*Metric(1,2)*Metric(3,4) - P(-1,2)*P(-1,3)*Metric(1,2)*Metric(3,4) - P(-1,1)*P(-1,4)*Metric(1,2)*Metric(3,4) - P(-1,2)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVV16 = Lorentz(name = 'VVVV16',
                 spins = [ 3, 3, 3, 3 ],
                 structure = 'P(3,2)*P(4,1)*Metric(1,2) - P(3,4)*P(4,1)*Metric(1,2) + P(3,1)*P(4,2)*Metric(1,2) - P(3,4)*P(4,2)*Metric(1,2) - P(3,1)*P(4,3)*Metric(1,2) - P(3,2)*P(4,3)*Metric(1,2) + P(2,3)*P(4,1)*Metric(1,3) - P(2,4)*P(4,1)*Metric(1,3) - P(2,1)*P(4,2)*Metric(1,3) - P(2,3)*P(4,2)*Metric(1,3) + P(2,1)*P(4,3)*Metric(1,3) - P(2,4)*P(4,3)*Metric(1,3) - P(2,3)*P(3,1)*Metric(1,4) + P(2,4)*P(3,1)*Metric(1,4) - P(2,1)*P(3,2)*Metric(1,4) - P(2,4)*P(3,2)*Metric(1,4) + P(2,1)*P(3,4)*Metric(1,4) - P(2,3)*P(3,4)*Metric(1,4) - P(1,2)*P(4,1)*Metric(2,3) - P(1,3)*P(4,1)*Metric(2,3) + P(1,3)*P(4,2)*Metric(2,3) - P(1,4)*P(4,2)*Metric(2,3) + P(1,2)*P(4,3)*Metric(2,3) - P(1,4)*P(4,3)*Metric(2,3) + P(-1,1)*P(-1,2)*Metric(1,4)*Metric(2,3) + P(-1,1)*P(-1,3)*Metric(1,4)*Metric(2,3) + P(-1,2)*P(-1,4)*Metric(1,4)*Metric(2,3) + P(-1,3)*P(-1,4)*Metric(1,4)*Metric(2,3) - P(1,2)*P(3,1)*Metric(2,4) - P(1,4)*P(3,1)*Metric(2,4) - P(1,3)*P(3,2)*Metric(2,4) + P(1,4)*P(3,2)*Metric(2,4) + P(1,2)*P(3,4)*Metric(2,4) - P(1,3)*P(3,4)*Metric(2,4) + P(-1,1)*P(-1,2)*Metric(1,3)*Metric(2,4) + P(-1,2)*P(-1,3)*Metric(1,3)*Metric(2,4) + P(-1,1)*P(-1,4)*Metric(1,3)*Metric(2,4) + P(-1,3)*P(-1,4)*Metric(1,3)*Metric(2,4) - P(1,3)*P(2,1)*Metric(3,4) - P(1,4)*P(2,1)*Metric(3,4) - P(1,2)*P(2,3)*Metric(3,4) + P(1,4)*P(2,3)*Metric(3,4) - P(1,2)*P(2,4)*Metric(3,4) + P(1,3)*P(2,4)*Metric(3,4) + P(-1,1)*P(-1,3)*Metric(1,2)*Metric(3,4) + P(-1,2)*P(-1,3)*Metric(1,2)*Metric(3,4) + P(-1,1)*P(-1,4)*Metric(1,2)*Metric(3,4) + P(-1,2)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVV17 = Lorentz(name = 'VVVV17',
                 spins = [ 3, 3, 3, 3 ],
                 structure = 'P(3,4)*P(4,3)*Metric(1,2) - P(-1,3)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVV18 = Lorentz(name = 'VVVV18',
                 spins = [ 3, 3, 3, 3 ],
                 structure = 'P(3,4)*P(4,3)*Metric(1,2) + P(1,2)*P(2,1)*Metric(3,4) - P(-1,1)*P(-1,2)*Metric(1,2)*Metric(3,4) - P(-1,3)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVV19 = Lorentz(name = 'VVVV19',
                 spins = [ 3, 3, 3, 3 ],
                 structure = 'P(3,4)*P(4,3)*Metric(1,2) + P(2,4)*P(4,2)*Metric(1,3) + P(2,3)*P(3,2)*Metric(1,4) + P(1,4)*P(4,1)*Metric(2,3) - P(-1,2)*P(-1,3)*Metric(1,4)*Metric(2,3) - P(-1,1)*P(-1,4)*Metric(1,4)*Metric(2,3) + P(1,3)*P(3,1)*Metric(2,4) - P(-1,1)*P(-1,3)*Metric(1,3)*Metric(2,4) - P(-1,2)*P(-1,4)*Metric(1,3)*Metric(2,4) + P(1,2)*P(2,1)*Metric(3,4) - P(-1,1)*P(-1,2)*Metric(1,2)*Metric(3,4) - P(-1,3)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVSS1 = Lorentz(name = 'VVVSS1',
                 spins = [ 3, 3, 3, 1, 1 ],
                 structure = 'P(1,5)*P(2,4)*P(3,1) + P(1,4)*P(2,5)*P(3,1) - P(1,5)*P(2,4)*P(3,2) - P(1,4)*P(2,5)*P(3,2) - P(1,5)*P(2,1)*P(3,4) + P(1,5)*P(2,3)*P(3,4) + P(1,2)*P(2,5)*P(3,4) - P(1,3)*P(2,5)*P(3,4) - P(1,4)*P(2,1)*P(3,5) + P(1,4)*P(2,3)*P(3,5) + P(1,2)*P(2,4)*P(3,5) - P(1,3)*P(2,4)*P(3,5) + P(-1,1)*P(-1,5)*P(3,4)*Metric(1,2) - P(-1,2)*P(-1,5)*P(3,4)*Metric(1,2) + P(-1,1)*P(-1,4)*P(3,5)*Metric(1,2) - P(-1,2)*P(-1,4)*P(3,5)*Metric(1,2) - P(-1,1)*P(-1,5)*P(2,4)*Metric(1,3) + P(-1,3)*P(-1,5)*P(2,4)*Metric(1,3) - P(-1,1)*P(-1,4)*P(2,5)*Metric(1,3) + P(-1,3)*P(-1,4)*P(2,5)*Metric(1,3) + P(-1,2)*P(-1,5)*P(1,4)*Metric(2,3) - P(-1,3)*P(-1,5)*P(1,4)*Metric(2,3) + P(-1,2)*P(-1,4)*P(1,5)*Metric(2,3) - P(-1,3)*P(-1,4)*P(1,5)*Metric(2,3)')

VVVSS2 = Lorentz(name = 'VVVSS2',
                 spins = [ 3, 3, 3, 1, 1 ],
                 structure = 'P(-1,4)*P(-1,5)*P(3,1)*Metric(1,2) - P(-1,4)*P(-1,5)*P(3,2)*Metric(1,2) - P(-1,4)*P(-1,5)*P(2,1)*Metric(1,3) + P(-1,4)*P(-1,5)*P(2,3)*Metric(1,3) + P(-1,4)*P(-1,5)*P(1,2)*Metric(2,3) - P(-1,4)*P(-1,5)*P(1,3)*Metric(2,3)')

VVVVS1 = Lorentz(name = 'VVVVS1',
                 spins = [ 3, 3, 3, 3, 1 ],
                 structure = 'P(1,4)*P(4,1)*Metric(2,3) - P(-1,1)*P(-1,4)*Metric(1,4)*Metric(2,3)')

VVVVS2 = Lorentz(name = 'VVVVS2',
                 spins = [ 3, 3, 3, 3, 1 ],
                 structure = 'P(3,2)*P(4,1)*Metric(1,2) + P(3,1)*P(4,2)*Metric(1,2) - P(2,1)*P(4,2)*Metric(1,3) - P(2,1)*P(3,2)*Metric(1,4) - P(1,2)*P(4,1)*Metric(2,3) + P(-1,1)*P(-1,2)*Metric(1,4)*Metric(2,3) - P(1,2)*P(3,1)*Metric(2,4) + P(-1,1)*P(-1,2)*Metric(1,3)*Metric(2,4)')

VVVVS3 = Lorentz(name = 'VVVVS3',
                 spins = [ 3, 3, 3, 3, 1 ],
                 structure = 'P(2,4)*P(4,2)*Metric(1,3) + P(2,3)*P(3,2)*Metric(1,4) + P(1,4)*P(4,1)*Metric(2,3) - P(-1,2)*P(-1,3)*Metric(1,4)*Metric(2,3) - P(-1,1)*P(-1,4)*Metric(1,4)*Metric(2,3) + P(1,3)*P(3,1)*Metric(2,4) - P(-1,1)*P(-1,3)*Metric(1,3)*Metric(2,4) - P(-1,2)*P(-1,4)*Metric(1,3)*Metric(2,4)')

VVVVS4 = Lorentz(name = 'VVVVS4',
                 spins = [ 3, 3, 3, 3, 1 ],
                 structure = 'P(2,4)*P(4,3)*Metric(1,3) + P(2,3)*P(3,4)*Metric(1,4) + P(1,4)*P(4,3)*Metric(2,3) - P(-1,3)*P(-1,4)*Metric(1,4)*Metric(2,3) + P(1,3)*P(3,4)*Metric(2,4) - P(-1,3)*P(-1,4)*Metric(1,3)*Metric(2,4) - P(1,4)*P(2,3)*Metric(3,4) - P(1,3)*P(2,4)*Metric(3,4)')

VVVVS5 = Lorentz(name = 'VVVVS5',
                 spins = [ 3, 3, 3, 3, 1 ],
                 structure = 'P(3,2)*P(4,1)*Metric(1,2) + P(3,1)*P(4,2)*Metric(1,2) - P(2,1)*P(4,2)*Metric(1,3) - P(2,4)*P(4,3)*Metric(1,3) - P(2,1)*P(3,2)*Metric(1,4) - P(2,3)*P(3,4)*Metric(1,4) - P(1,2)*P(4,1)*Metric(2,3) - P(1,4)*P(4,3)*Metric(2,3) + P(-1,1)*P(-1,2)*Metric(1,4)*Metric(2,3) + P(-1,3)*P(-1,4)*Metric(1,4)*Metric(2,3) - P(1,2)*P(3,1)*Metric(2,4) - P(1,3)*P(3,4)*Metric(2,4) + P(-1,1)*P(-1,2)*Metric(1,3)*Metric(2,4) + P(-1,3)*P(-1,4)*Metric(1,3)*Metric(2,4) + P(1,4)*P(2,3)*Metric(3,4) + P(1,3)*P(2,4)*Metric(3,4)')

VVVVS6 = Lorentz(name = 'VVVVS6',
                 spins = [ 3, 3, 3, 3, 1 ],
                 structure = 'P(1,2)*P(2,1)*Metric(3,4) - P(-1,1)*P(-1,2)*Metric(1,2)*Metric(3,4)')

VVVVS7 = Lorentz(name = 'VVVVS7',
                 spins = [ 3, 3, 3, 3, 1 ],
                 structure = 'P(1,4)*P(4,1)*Metric(2,3) - P(-1,1)*P(-1,4)*Metric(1,4)*Metric(2,3) + P(1,3)*P(3,1)*Metric(2,4) - P(-1,1)*P(-1,3)*Metric(1,3)*Metric(2,4) + P(1,2)*P(2,1)*Metric(3,4) - P(-1,1)*P(-1,2)*Metric(1,2)*Metric(3,4)')

VVVVS8 = Lorentz(name = 'VVVVS8',
                 spins = [ 3, 3, 3, 3, 1 ],
                 structure = 'P(3,4)*P(4,1)*Metric(1,2) + P(2,4)*P(4,1)*Metric(1,3) - P(2,4)*P(3,1)*Metric(1,4) - P(2,1)*P(3,4)*Metric(1,4) + P(1,4)*P(3,1)*Metric(2,4) - P(-1,1)*P(-1,4)*Metric(1,3)*Metric(2,4) + P(1,4)*P(2,1)*Metric(3,4) - P(-1,1)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVVS9 = Lorentz(name = 'VVVVS9',
                 spins = [ 3, 3, 3, 3, 1 ],
                 structure = 'P(3,2)*P(4,1)*Metric(1,2) - P(3,4)*P(4,1)*Metric(1,2) + P(3,1)*P(4,2)*Metric(1,2) - P(3,1)*P(4,3)*Metric(1,2) + P(2,3)*P(4,1)*Metric(1,3) - P(2,4)*P(4,1)*Metric(1,3) - P(2,1)*P(4,2)*Metric(1,3) + P(2,1)*P(4,3)*Metric(1,3) - P(2,3)*P(3,1)*Metric(1,4) + P(2,4)*P(3,1)*Metric(1,4) - P(2,1)*P(3,2)*Metric(1,4) + P(2,1)*P(3,4)*Metric(1,4) - P(1,2)*P(4,1)*Metric(2,3) - P(1,3)*P(4,1)*Metric(2,3) + P(-1,1)*P(-1,2)*Metric(1,4)*Metric(2,3) + P(-1,1)*P(-1,3)*Metric(1,4)*Metric(2,3) - P(1,2)*P(3,1)*Metric(2,4) - P(1,4)*P(3,1)*Metric(2,4) + P(-1,1)*P(-1,2)*Metric(1,3)*Metric(2,4) + P(-1,1)*P(-1,4)*Metric(1,3)*Metric(2,4) - P(1,3)*P(2,1)*Metric(3,4) - P(1,4)*P(2,1)*Metric(3,4) + P(-1,1)*P(-1,3)*Metric(1,2)*Metric(3,4) + P(-1,1)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVVS10 = Lorentz(name = 'VVVVS10',
                  spins = [ 3, 3, 3, 3, 1 ],
                  structure = 'P(3,4)*P(4,1)*Metric(1,2) + P(3,4)*P(4,2)*Metric(1,2) + P(3,1)*P(4,3)*Metric(1,2) + P(3,2)*P(4,3)*Metric(1,2) - P(2,3)*P(4,1)*Metric(1,3) + P(2,4)*P(4,1)*Metric(1,3) + P(2,3)*P(4,2)*Metric(1,3) - P(2,1)*P(4,3)*Metric(1,3) + P(2,3)*P(3,1)*Metric(1,4) - P(2,4)*P(3,1)*Metric(1,4) + P(2,4)*P(3,2)*Metric(1,4) - P(2,1)*P(3,4)*Metric(1,4) + P(1,3)*P(4,1)*Metric(2,3) - P(1,3)*P(4,2)*Metric(2,3) + P(1,4)*P(4,2)*Metric(2,3) - P(1,2)*P(4,3)*Metric(2,3) - P(-1,1)*P(-1,3)*Metric(1,4)*Metric(2,3) - P(-1,2)*P(-1,4)*Metric(1,4)*Metric(2,3) + P(1,4)*P(3,1)*Metric(2,4) + P(1,3)*P(3,2)*Metric(2,4) - P(1,4)*P(3,2)*Metric(2,4) - P(1,2)*P(3,4)*Metric(2,4) - P(-1,2)*P(-1,3)*Metric(1,3)*Metric(2,4) - P(-1,1)*P(-1,4)*Metric(1,3)*Metric(2,4) + P(1,3)*P(2,1)*Metric(3,4) + P(1,4)*P(2,1)*Metric(3,4) + P(1,2)*P(2,3)*Metric(3,4) + P(1,2)*P(2,4)*Metric(3,4) - P(-1,1)*P(-1,3)*Metric(1,2)*Metric(3,4) - P(-1,2)*P(-1,3)*Metric(1,2)*Metric(3,4) - P(-1,1)*P(-1,4)*Metric(1,2)*Metric(3,4) - P(-1,2)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVVS11 = Lorentz(name = 'VVVVS11',
                  spins = [ 3, 3, 3, 3, 1 ],
                  structure = 'P(3,2)*P(4,1)*Metric(1,2) - P(3,4)*P(4,1)*Metric(1,2) + P(3,1)*P(4,2)*Metric(1,2) - P(3,4)*P(4,2)*Metric(1,2) - P(3,1)*P(4,3)*Metric(1,2) - P(3,2)*P(4,3)*Metric(1,2) + P(2,3)*P(4,1)*Metric(1,3) - P(2,4)*P(4,1)*Metric(1,3) - P(2,1)*P(4,2)*Metric(1,3) - P(2,3)*P(4,2)*Metric(1,3) + P(2,1)*P(4,3)*Metric(1,3) - P(2,4)*P(4,3)*Metric(1,3) - P(2,3)*P(3,1)*Metric(1,4) + P(2,4)*P(3,1)*Metric(1,4) - P(2,1)*P(3,2)*Metric(1,4) - P(2,4)*P(3,2)*Metric(1,4) + P(2,1)*P(3,4)*Metric(1,4) - P(2,3)*P(3,4)*Metric(1,4) - P(1,2)*P(4,1)*Metric(2,3) - P(1,3)*P(4,1)*Metric(2,3) + P(1,3)*P(4,2)*Metric(2,3) - P(1,4)*P(4,2)*Metric(2,3) + P(1,2)*P(4,3)*Metric(2,3) - P(1,4)*P(4,3)*Metric(2,3) + P(-1,1)*P(-1,2)*Metric(1,4)*Metric(2,3) + P(-1,1)*P(-1,3)*Metric(1,4)*Metric(2,3) + P(-1,2)*P(-1,4)*Metric(1,4)*Metric(2,3) + P(-1,3)*P(-1,4)*Metric(1,4)*Metric(2,3) - P(1,2)*P(3,1)*Metric(2,4) - P(1,4)*P(3,1)*Metric(2,4) - P(1,3)*P(3,2)*Metric(2,4) + P(1,4)*P(3,2)*Metric(2,4) + P(1,2)*P(3,4)*Metric(2,4) - P(1,3)*P(3,4)*Metric(2,4) + P(-1,1)*P(-1,2)*Metric(1,3)*Metric(2,4) + P(-1,2)*P(-1,3)*Metric(1,3)*Metric(2,4) + P(-1,1)*P(-1,4)*Metric(1,3)*Metric(2,4) + P(-1,3)*P(-1,4)*Metric(1,3)*Metric(2,4) - P(1,3)*P(2,1)*Metric(3,4) - P(1,4)*P(2,1)*Metric(3,4) - P(1,2)*P(2,3)*Metric(3,4) + P(1,4)*P(2,3)*Metric(3,4) - P(1,2)*P(2,4)*Metric(3,4) + P(1,3)*P(2,4)*Metric(3,4) + P(-1,1)*P(-1,3)*Metric(1,2)*Metric(3,4) + P(-1,2)*P(-1,3)*Metric(1,2)*Metric(3,4) + P(-1,1)*P(-1,4)*Metric(1,2)*Metric(3,4) + P(-1,2)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVVS12 = Lorentz(name = 'VVVVS12',
                  spins = [ 3, 3, 3, 3, 1 ],
                  structure = 'P(3,4)*P(4,3)*Metric(1,2) - P(-1,3)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVVS13 = Lorentz(name = 'VVVVS13',
                  spins = [ 3, 3, 3, 3, 1 ],
                  structure = 'P(3,4)*P(4,3)*Metric(1,2) + P(1,2)*P(2,1)*Metric(3,4) - P(-1,1)*P(-1,2)*Metric(1,2)*Metric(3,4) - P(-1,3)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVVS14 = Lorentz(name = 'VVVVS14',
                  spins = [ 3, 3, 3, 3, 1 ],
                  structure = 'P(3,4)*P(4,3)*Metric(1,2) + P(2,4)*P(4,2)*Metric(1,3) + P(2,3)*P(3,2)*Metric(1,4) + P(1,4)*P(4,1)*Metric(2,3) - P(-1,2)*P(-1,3)*Metric(1,4)*Metric(2,3) - P(-1,1)*P(-1,4)*Metric(1,4)*Metric(2,3) + P(1,3)*P(3,1)*Metric(2,4) - P(-1,1)*P(-1,3)*Metric(1,3)*Metric(2,4) - P(-1,2)*P(-1,4)*Metric(1,3)*Metric(2,4) + P(1,2)*P(2,1)*Metric(3,4) - P(-1,1)*P(-1,2)*Metric(1,2)*Metric(3,4) - P(-1,3)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVVV1 = Lorentz(name = 'VVVVV1',
                 spins = [ 3, 3, 3, 3, 3 ],
                 structure = 'P(5,2)*Metric(1,4)*Metric(2,3) - P(5,3)*Metric(1,4)*Metric(2,3) + P(4,2)*Metric(1,5)*Metric(2,3) - P(4,3)*Metric(1,5)*Metric(2,3) - P(5,1)*Metric(1,3)*Metric(2,4) + P(5,3)*Metric(1,3)*Metric(2,4) + P(3,1)*Metric(1,5)*Metric(2,4) - P(3,2)*Metric(1,5)*Metric(2,4) - P(4,1)*Metric(1,3)*Metric(2,5) + P(4,3)*Metric(1,3)*Metric(2,5) + P(3,1)*Metric(1,4)*Metric(2,5) - P(3,2)*Metric(1,4)*Metric(2,5) + P(5,1)*Metric(1,2)*Metric(3,4) - P(5,2)*Metric(1,2)*Metric(3,4) - P(2,1)*Metric(1,5)*Metric(3,4) + P(2,3)*Metric(1,5)*Metric(3,4) + P(1,2)*Metric(2,5)*Metric(3,4) - P(1,3)*Metric(2,5)*Metric(3,4) + P(4,1)*Metric(1,2)*Metric(3,5) - P(4,2)*Metric(1,2)*Metric(3,5) - P(2,1)*Metric(1,4)*Metric(3,5) + P(2,3)*Metric(1,4)*Metric(3,5) + P(1,2)*Metric(2,4)*Metric(3,5) - P(1,3)*Metric(2,4)*Metric(3,5)')

VVVVV2 = Lorentz(name = 'VVVVV2',
                 spins = [ 3, 3, 3, 3, 3 ],
                 structure = 'P(5,1)*Metric(1,3)*Metric(2,4) - P(5,3)*Metric(1,3)*Metric(2,4) - P(3,1)*Metric(1,5)*Metric(2,4) + P(3,5)*Metric(1,5)*Metric(2,4) + P(4,1)*Metric(1,3)*Metric(2,5) - P(4,3)*Metric(1,3)*Metric(2,5) - P(3,1)*Metric(1,4)*Metric(2,5) + P(3,4)*Metric(1,4)*Metric(2,5) + P(5,1)*Metric(1,2)*Metric(3,4) - P(5,2)*Metric(1,2)*Metric(3,4) - P(2,1)*Metric(1,5)*Metric(3,4) + P(2,5)*Metric(1,5)*Metric(3,4) + P(1,2)*Metric(2,5)*Metric(3,4) + P(1,3)*Metric(2,5)*Metric(3,4) - P(1,4)*Metric(2,5)*Metric(3,4) - P(1,5)*Metric(2,5)*Metric(3,4) + P(4,1)*Metric(1,2)*Metric(3,5) - P(4,2)*Metric(1,2)*Metric(3,5) - P(2,1)*Metric(1,4)*Metric(3,5) + P(2,4)*Metric(1,4)*Metric(3,5) + P(1,2)*Metric(2,4)*Metric(3,5) + P(1,3)*Metric(2,4)*Metric(3,5) - P(1,4)*Metric(2,4)*Metric(3,5) - P(1,5)*Metric(2,4)*Metric(3,5)')

VVVVV3 = Lorentz(name = 'VVVVV3',
                 spins = [ 3, 3, 3, 3, 3 ],
                 structure = 'P(3,1)*Metric(1,2)*Metric(4,5) - P(3,2)*Metric(1,2)*Metric(4,5) - P(2,1)*Metric(1,3)*Metric(4,5) + P(2,3)*Metric(1,3)*Metric(4,5) + P(1,2)*Metric(2,3)*Metric(4,5) - P(1,3)*Metric(2,3)*Metric(4,5)')

VVVVV4 = Lorentz(name = 'VVVVV4',
                 spins = [ 3, 3, 3, 3, 3 ],
                 structure = 'P(5,1)*Metric(1,2)*Metric(3,4) - P(5,2)*Metric(1,2)*Metric(3,4) - P(2,1)*Metric(1,5)*Metric(3,4) + P(2,5)*Metric(1,5)*Metric(3,4) + P(1,2)*Metric(2,5)*Metric(3,4) - P(1,5)*Metric(2,5)*Metric(3,4) + P(4,1)*Metric(1,2)*Metric(3,5) - P(4,2)*Metric(1,2)*Metric(3,5) - P(2,1)*Metric(1,4)*Metric(3,5) + P(2,4)*Metric(1,4)*Metric(3,5) + P(1,2)*Metric(2,4)*Metric(3,5) - P(1,4)*Metric(2,4)*Metric(3,5) + P(3,1)*Metric(1,2)*Metric(4,5) - P(3,2)*Metric(1,2)*Metric(4,5) - P(2,1)*Metric(1,3)*Metric(4,5) + P(2,3)*Metric(1,3)*Metric(4,5) + P(1,2)*Metric(2,3)*Metric(4,5) - P(1,3)*Metric(2,3)*Metric(4,5)')

VVVVV5 = Lorentz(name = 'VVVVV5',
                 spins = [ 3, 3, 3, 3, 3 ],
                 structure = 'P(5,1)*Metric(1,4)*Metric(2,3) + P(5,2)*Metric(1,4)*Metric(2,3) - P(5,3)*Metric(1,4)*Metric(2,3) - P(5,4)*Metric(1,4)*Metric(2,3) - P(4,1)*Metric(1,5)*Metric(2,3) - P(4,3)*Metric(1,5)*Metric(2,3) + 2*P(4,5)*Metric(1,5)*Metric(2,3) + P(5,1)*Metric(1,3)*Metric(2,4) + P(5,2)*Metric(1,3)*Metric(2,4) - P(5,3)*Metric(1,3)*Metric(2,4) - P(5,4)*Metric(1,3)*Metric(2,4) - P(3,1)*Metric(1,5)*Metric(2,4) - P(3,4)*Metric(1,5)*Metric(2,4) + 2*P(3,5)*Metric(1,5)*Metric(2,4) - P(4,2)*Metric(1,3)*Metric(2,5) - P(4,3)*Metric(1,3)*Metric(2,5) + 2*P(4,5)*Metric(1,3)*Metric(2,5) - P(3,2)*Metric(1,4)*Metric(2,5) - P(3,4)*Metric(1,4)*Metric(2,5) + 2*P(3,5)*Metric(1,4)*Metric(2,5) + 2*P(5,1)*Metric(1,2)*Metric(3,4) + 2*P(5,2)*Metric(1,2)*Metric(3,4) - 2*P(5,3)*Metric(1,2)*Metric(3,4) - 2*P(5,4)*Metric(1,2)*Metric(3,4) - 2*P(2,1)*Metric(1,5)*Metric(3,4) + P(2,3)*Metric(1,5)*Metric(3,4) + P(2,4)*Metric(1,5)*Metric(3,4) - 2*P(1,2)*Metric(2,5)*Metric(3,4) + P(1,3)*Metric(2,5)*Metric(3,4) + P(1,4)*Metric(2,5)*Metric(3,4) - P(4,1)*Metric(1,2)*Metric(3,5) - P(4,2)*Metric(1,2)*Metric(3,5) + 2*P(4,3)*Metric(1,2)*Metric(3,5) + P(2,1)*Metric(1,4)*Metric(3,5) + P(2,3)*Metric(1,4)*Metric(3,5) - 2*P(2,5)*Metric(1,4)*Metric(3,5) + P(1,2)*Metric(2,4)*Metric(3,5) + P(1,3)*Metric(2,4)*Metric(3,5) - 2*P(1,5)*Metric(2,4)*Metric(3,5) - P(3,1)*Metric(1,2)*Metric(4,5) - P(3,2)*Metric(1,2)*Metric(4,5) + 2*P(3,4)*Metric(1,2)*Metric(4,5) + P(2,1)*Metric(1,3)*Metric(4,5) + P(2,4)*Metric(1,3)*Metric(4,5) - 2*P(2,5)*Metric(1,3)*Metric(4,5) + P(1,2)*Metric(2,3)*Metric(4,5) + P(1,4)*Metric(2,3)*Metric(4,5) - 2*P(1,5)*Metric(2,3)*Metric(4,5)')

VVVVV6 = Lorentz(name = 'VVVVV6',
                 spins = [ 3, 3, 3, 3, 3 ],
                 structure = 'P(5,1)*Metric(1,4)*Metric(2,3) + P(5,2)*Metric(1,4)*Metric(2,3) - P(5,3)*Metric(1,4)*Metric(2,3) - P(5,4)*Metric(1,4)*Metric(2,3) - P(4,1)*Metric(1,5)*Metric(2,3) + P(4,5)*Metric(1,5)*Metric(2,3) + P(5,1)*Metric(1,3)*Metric(2,4) + P(5,2)*Metric(1,3)*Metric(2,4) - P(5,3)*Metric(1,3)*Metric(2,4) - P(5,4)*Metric(1,3)*Metric(2,4) - P(3,1)*Metric(1,5)*Metric(2,4) + P(3,5)*Metric(1,5)*Metric(2,4) - P(4,2)*Metric(1,3)*Metric(2,5) + P(4,5)*Metric(1,3)*Metric(2,5) - P(3,2)*Metric(1,4)*Metric(2,5) + P(3,5)*Metric(1,4)*Metric(2,5) + P(2,3)*Metric(1,4)*Metric(3,5) - P(2,5)*Metric(1,4)*Metric(3,5) + P(1,3)*Metric(2,4)*Metric(3,5) - P(1,5)*Metric(2,4)*Metric(3,5) + P(2,4)*Metric(1,3)*Metric(4,5) - P(2,5)*Metric(1,3)*Metric(4,5) + P(1,4)*Metric(2,3)*Metric(4,5) - P(1,5)*Metric(2,3)*Metric(4,5)')

VVVVV7 = Lorentz(name = 'VVVVV7',
                 spins = [ 3, 3, 3, 3, 3 ],
                 structure = 'P(5,3)*Metric(1,4)*Metric(2,3) - P(5,4)*Metric(1,4)*Metric(2,3) + P(4,3)*Metric(1,5)*Metric(2,3) - P(4,5)*Metric(1,5)*Metric(2,3) - P(5,3)*Metric(1,3)*Metric(2,4) + P(5,4)*Metric(1,3)*Metric(2,4) + P(3,4)*Metric(1,5)*Metric(2,4) - P(3,5)*Metric(1,5)*Metric(2,4) - P(4,3)*Metric(1,3)*Metric(2,5) + P(4,5)*Metric(1,3)*Metric(2,5) - P(3,4)*Metric(1,4)*Metric(2,5) + P(3,5)*Metric(1,4)*Metric(2,5) - 2*P(5,1)*Metric(1,2)*Metric(3,4) + 2*P(5,2)*Metric(1,2)*Metric(3,4) + 2*P(2,1)*Metric(1,5)*Metric(3,4) - P(2,3)*Metric(1,5)*Metric(3,4) - P(2,4)*Metric(1,5)*Metric(3,4) - 2*P(1,2)*Metric(2,5)*Metric(3,4) + P(1,3)*Metric(2,5)*Metric(3,4) + P(1,4)*Metric(2,5)*Metric(3,4) - 2*P(4,1)*Metric(1,2)*Metric(3,5) + 2*P(4,2)*Metric(1,2)*Metric(3,5) + 2*P(2,1)*Metric(1,4)*Metric(3,5) - P(2,3)*Metric(1,4)*Metric(3,5) - P(2,5)*Metric(1,4)*Metric(3,5) - 2*P(1,2)*Metric(2,4)*Metric(3,5) + P(1,3)*Metric(2,4)*Metric(3,5) + P(1,5)*Metric(2,4)*Metric(3,5) - 2*P(3,1)*Metric(1,2)*Metric(4,5) + 2*P(3,2)*Metric(1,2)*Metric(4,5) + 2*P(2,1)*Metric(1,3)*Metric(4,5) - P(2,4)*Metric(1,3)*Metric(4,5) - P(2,5)*Metric(1,3)*Metric(4,5) - 2*P(1,2)*Metric(2,3)*Metric(4,5) + P(1,4)*Metric(2,3)*Metric(4,5) + P(1,5)*Metric(2,3)*Metric(4,5)')

VVVVV8 = Lorentz(name = 'VVVVV8',
                 spins = [ 3, 3, 3, 3, 3 ],
                 structure = 'P(5,2)*Metric(1,4)*Metric(2,3) + P(5,3)*Metric(1,4)*Metric(2,3) - 2*P(5,4)*Metric(1,4)*Metric(2,3) + P(4,2)*Metric(1,5)*Metric(2,3) + P(4,3)*Metric(1,5)*Metric(2,3) - 2*P(4,5)*Metric(1,5)*Metric(2,3) - 2*P(5,1)*Metric(1,3)*Metric(2,4) + P(5,3)*Metric(1,3)*Metric(2,4) + P(5,4)*Metric(1,3)*Metric(2,4) + 2*P(3,1)*Metric(1,5)*Metric(2,4) - P(3,2)*Metric(1,5)*Metric(2,4) - P(3,5)*Metric(1,5)*Metric(2,4) - 2*P(4,1)*Metric(1,3)*Metric(2,5) + P(4,3)*Metric(1,3)*Metric(2,5) + P(4,5)*Metric(1,3)*Metric(2,5) + 2*P(3,1)*Metric(1,4)*Metric(2,5) - P(3,2)*Metric(1,4)*Metric(2,5) - P(3,4)*Metric(1,4)*Metric(2,5) - 2*P(5,1)*Metric(1,2)*Metric(3,4) + P(5,2)*Metric(1,2)*Metric(3,4) + P(5,4)*Metric(1,2)*Metric(3,4) + 2*P(2,1)*Metric(1,5)*Metric(3,4) - P(2,3)*Metric(1,5)*Metric(3,4) - P(2,5)*Metric(1,5)*Metric(3,4) - P(1,2)*Metric(2,5)*Metric(3,4) - P(1,3)*Metric(2,5)*Metric(3,4) + P(1,4)*Metric(2,5)*Metric(3,4) + P(1,5)*Metric(2,5)*Metric(3,4) - 2*P(4,1)*Metric(1,2)*Metric(3,5) + P(4,2)*Metric(1,2)*Metric(3,5) + P(4,5)*Metric(1,2)*Metric(3,5) + 2*P(2,1)*Metric(1,4)*Metric(3,5) - P(2,3)*Metric(1,4)*Metric(3,5) - P(2,4)*Metric(1,4)*Metric(3,5) - P(1,2)*Metric(2,4)*Metric(3,5) - P(1,3)*Metric(2,4)*Metric(3,5) + P(1,4)*Metric(2,4)*Metric(3,5) + P(1,5)*Metric(2,4)*Metric(3,5) + 2*P(3,2)*Metric(1,2)*Metric(4,5) - P(3,4)*Metric(1,2)*Metric(4,5) - P(3,5)*Metric(1,2)*Metric(4,5) + 2*P(2,3)*Metric(1,3)*Metric(4,5) - P(2,4)*Metric(1,3)*Metric(4,5) - P(2,5)*Metric(1,3)*Metric(4,5) - 2*P(1,2)*Metric(2,3)*Metric(4,5) - 2*P(1,3)*Metric(2,3)*Metric(4,5) + 2*P(1,4)*Metric(2,3)*Metric(4,5) + 2*P(1,5)*Metric(2,3)*Metric(4,5)')

VVVVSS1 = Lorentz(name = 'VVVVSS1',
                  spins = [ 3, 3, 3, 3, 1, 1 ],
                  structure = 'P(1,4)*P(4,1)*Metric(2,3) - P(-1,1)*P(-1,4)*Metric(1,4)*Metric(2,3)')

VVVVSS2 = Lorentz(name = 'VVVVSS2',
                  spins = [ 3, 3, 3, 3, 1, 1 ],
                  structure = 'P(3,2)*P(4,1)*Metric(1,2) + P(3,1)*P(4,2)*Metric(1,2) - P(2,1)*P(4,2)*Metric(1,3) - P(2,1)*P(3,2)*Metric(1,4) - P(1,2)*P(4,1)*Metric(2,3) + P(-1,1)*P(-1,2)*Metric(1,4)*Metric(2,3) - P(1,2)*P(3,1)*Metric(2,4) + P(-1,1)*P(-1,2)*Metric(1,3)*Metric(2,4)')

VVVVSS3 = Lorentz(name = 'VVVVSS3',
                  spins = [ 3, 3, 3, 3, 1, 1 ],
                  structure = 'P(2,4)*P(4,2)*Metric(1,3) + P(2,3)*P(3,2)*Metric(1,4) + P(1,4)*P(4,1)*Metric(2,3) - P(-1,2)*P(-1,3)*Metric(1,4)*Metric(2,3) - P(-1,1)*P(-1,4)*Metric(1,4)*Metric(2,3) + P(1,3)*P(3,1)*Metric(2,4) - P(-1,1)*P(-1,3)*Metric(1,3)*Metric(2,4) - P(-1,2)*P(-1,4)*Metric(1,3)*Metric(2,4)')

VVVVSS4 = Lorentz(name = 'VVVVSS4',
                  spins = [ 3, 3, 3, 3, 1, 1 ],
                  structure = 'P(2,4)*P(4,3)*Metric(1,3) + P(2,3)*P(3,4)*Metric(1,4) + P(1,4)*P(4,3)*Metric(2,3) - P(-1,3)*P(-1,4)*Metric(1,4)*Metric(2,3) + P(1,3)*P(3,4)*Metric(2,4) - P(-1,3)*P(-1,4)*Metric(1,3)*Metric(2,4) - P(1,4)*P(2,3)*Metric(3,4) - P(1,3)*P(2,4)*Metric(3,4)')

VVVVSS5 = Lorentz(name = 'VVVVSS5',
                  spins = [ 3, 3, 3, 3, 1, 1 ],
                  structure = 'P(3,2)*P(4,1)*Metric(1,2) + P(3,1)*P(4,2)*Metric(1,2) - P(2,1)*P(4,2)*Metric(1,3) - P(2,4)*P(4,3)*Metric(1,3) - P(2,1)*P(3,2)*Metric(1,4) - P(2,3)*P(3,4)*Metric(1,4) - P(1,2)*P(4,1)*Metric(2,3) - P(1,4)*P(4,3)*Metric(2,3) + P(-1,1)*P(-1,2)*Metric(1,4)*Metric(2,3) + P(-1,3)*P(-1,4)*Metric(1,4)*Metric(2,3) - P(1,2)*P(3,1)*Metric(2,4) - P(1,3)*P(3,4)*Metric(2,4) + P(-1,1)*P(-1,2)*Metric(1,3)*Metric(2,4) + P(-1,3)*P(-1,4)*Metric(1,3)*Metric(2,4) + P(1,4)*P(2,3)*Metric(3,4) + P(1,3)*P(2,4)*Metric(3,4)')

VVVVSS6 = Lorentz(name = 'VVVVSS6',
                  spins = [ 3, 3, 3, 3, 1, 1 ],
                  structure = 'P(3,6)*P(4,5)*Metric(1,2) + P(3,5)*P(4,6)*Metric(1,2) - (P(2,6)*P(4,5)*Metric(1,3))/2. - (P(2,5)*P(4,6)*Metric(1,3))/2. - (P(2,6)*P(3,5)*Metric(1,4))/2. - (P(2,5)*P(3,6)*Metric(1,4))/2. - (P(1,6)*P(4,5)*Metric(2,3))/2. - (P(1,5)*P(4,6)*Metric(2,3))/2. - (P(1,6)*P(3,5)*Metric(2,4))/2. - (P(1,5)*P(3,6)*Metric(2,4))/2. + P(1,6)*P(2,5)*Metric(3,4) + P(1,5)*P(2,6)*Metric(3,4)')

VVVVSS7 = Lorentz(name = 'VVVVSS7',
                  spins = [ 3, 3, 3, 3, 1, 1 ],
                  structure = 'P(3,6)*P(4,5)*Metric(1,2) + P(3,5)*P(4,6)*Metric(1,2) + P(2,6)*P(4,5)*Metric(1,3) + P(2,5)*P(4,6)*Metric(1,3) - 2*P(2,6)*P(3,5)*Metric(1,4) - 2*P(2,5)*P(3,6)*Metric(1,4) - 2*P(1,6)*P(4,5)*Metric(2,3) - 2*P(1,5)*P(4,6)*Metric(2,3) + P(1,6)*P(3,5)*Metric(2,4) + P(1,5)*P(3,6)*Metric(2,4) + P(1,6)*P(2,5)*Metric(3,4) + P(1,5)*P(2,6)*Metric(3,4)')

VVVVSS8 = Lorentz(name = 'VVVVSS8',
                  spins = [ 3, 3, 3, 3, 1, 1 ],
                  structure = 'P(1,2)*P(2,1)*Metric(3,4) - P(-1,1)*P(-1,2)*Metric(1,2)*Metric(3,4)')

VVVVSS9 = Lorentz(name = 'VVVVSS9',
                  spins = [ 3, 3, 3, 3, 1, 1 ],
                  structure = 'P(1,4)*P(4,1)*Metric(2,3) - P(-1,1)*P(-1,4)*Metric(1,4)*Metric(2,3) + P(1,3)*P(3,1)*Metric(2,4) - P(-1,1)*P(-1,3)*Metric(1,3)*Metric(2,4) + P(1,2)*P(2,1)*Metric(3,4) - P(-1,1)*P(-1,2)*Metric(1,2)*Metric(3,4)')

VVVVSS10 = Lorentz(name = 'VVVVSS10',
                   spins = [ 3, 3, 3, 3, 1, 1 ],
                   structure = 'P(3,4)*P(4,1)*Metric(1,2) + P(2,4)*P(4,1)*Metric(1,3) - P(2,4)*P(3,1)*Metric(1,4) - P(2,1)*P(3,4)*Metric(1,4) + P(1,4)*P(3,1)*Metric(2,4) - P(-1,1)*P(-1,4)*Metric(1,3)*Metric(2,4) + P(1,4)*P(2,1)*Metric(3,4) - P(-1,1)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVVSS11 = Lorentz(name = 'VVVVSS11',
                   spins = [ 3, 3, 3, 3, 1, 1 ],
                   structure = 'P(3,2)*P(4,1)*Metric(1,2) - P(3,4)*P(4,1)*Metric(1,2) + P(3,1)*P(4,2)*Metric(1,2) - P(3,1)*P(4,3)*Metric(1,2) + P(2,3)*P(4,1)*Metric(1,3) - P(2,4)*P(4,1)*Metric(1,3) - P(2,1)*P(4,2)*Metric(1,3) + P(2,1)*P(4,3)*Metric(1,3) - P(2,3)*P(3,1)*Metric(1,4) + P(2,4)*P(3,1)*Metric(1,4) - P(2,1)*P(3,2)*Metric(1,4) + P(2,1)*P(3,4)*Metric(1,4) - P(1,2)*P(4,1)*Metric(2,3) - P(1,3)*P(4,1)*Metric(2,3) + P(-1,1)*P(-1,2)*Metric(1,4)*Metric(2,3) + P(-1,1)*P(-1,3)*Metric(1,4)*Metric(2,3) - P(1,2)*P(3,1)*Metric(2,4) - P(1,4)*P(3,1)*Metric(2,4) + P(-1,1)*P(-1,2)*Metric(1,3)*Metric(2,4) + P(-1,1)*P(-1,4)*Metric(1,3)*Metric(2,4) - P(1,3)*P(2,1)*Metric(3,4) - P(1,4)*P(2,1)*Metric(3,4) + P(-1,1)*P(-1,3)*Metric(1,2)*Metric(3,4) + P(-1,1)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVVSS12 = Lorentz(name = 'VVVVSS12',
                   spins = [ 3, 3, 3, 3, 1, 1 ],
                   structure = 'P(3,4)*P(4,1)*Metric(1,2) + P(3,4)*P(4,2)*Metric(1,2) + P(3,1)*P(4,3)*Metric(1,2) + P(3,2)*P(4,3)*Metric(1,2) - P(2,3)*P(4,1)*Metric(1,3) + P(2,4)*P(4,1)*Metric(1,3) + P(2,3)*P(4,2)*Metric(1,3) - P(2,1)*P(4,3)*Metric(1,3) + P(2,3)*P(3,1)*Metric(1,4) - P(2,4)*P(3,1)*Metric(1,4) + P(2,4)*P(3,2)*Metric(1,4) - P(2,1)*P(3,4)*Metric(1,4) + P(1,3)*P(4,1)*Metric(2,3) - P(1,3)*P(4,2)*Metric(2,3) + P(1,4)*P(4,2)*Metric(2,3) - P(1,2)*P(4,3)*Metric(2,3) - P(-1,1)*P(-1,3)*Metric(1,4)*Metric(2,3) - P(-1,2)*P(-1,4)*Metric(1,4)*Metric(2,3) + P(1,4)*P(3,1)*Metric(2,4) + P(1,3)*P(3,2)*Metric(2,4) - P(1,4)*P(3,2)*Metric(2,4) - P(1,2)*P(3,4)*Metric(2,4) - P(-1,2)*P(-1,3)*Metric(1,3)*Metric(2,4) - P(-1,1)*P(-1,4)*Metric(1,3)*Metric(2,4) + P(1,3)*P(2,1)*Metric(3,4) + P(1,4)*P(2,1)*Metric(3,4) + P(1,2)*P(2,3)*Metric(3,4) + P(1,2)*P(2,4)*Metric(3,4) - P(-1,1)*P(-1,3)*Metric(1,2)*Metric(3,4) - P(-1,2)*P(-1,3)*Metric(1,2)*Metric(3,4) - P(-1,1)*P(-1,4)*Metric(1,2)*Metric(3,4) - P(-1,2)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVVSS13 = Lorentz(name = 'VVVVSS13',
                   spins = [ 3, 3, 3, 3, 1, 1 ],
                   structure = 'P(3,2)*P(4,1)*Metric(1,2) - P(3,4)*P(4,1)*Metric(1,2) + P(3,1)*P(4,2)*Metric(1,2) - P(3,4)*P(4,2)*Metric(1,2) - P(3,1)*P(4,3)*Metric(1,2) - P(3,2)*P(4,3)*Metric(1,2) + P(2,3)*P(4,1)*Metric(1,3) - P(2,4)*P(4,1)*Metric(1,3) - P(2,1)*P(4,2)*Metric(1,3) - P(2,3)*P(4,2)*Metric(1,3) + P(2,1)*P(4,3)*Metric(1,3) - P(2,4)*P(4,3)*Metric(1,3) - P(2,3)*P(3,1)*Metric(1,4) + P(2,4)*P(3,1)*Metric(1,4) - P(2,1)*P(3,2)*Metric(1,4) - P(2,4)*P(3,2)*Metric(1,4) + P(2,1)*P(3,4)*Metric(1,4) - P(2,3)*P(3,4)*Metric(1,4) - P(1,2)*P(4,1)*Metric(2,3) - P(1,3)*P(4,1)*Metric(2,3) + P(1,3)*P(4,2)*Metric(2,3) - P(1,4)*P(4,2)*Metric(2,3) + P(1,2)*P(4,3)*Metric(2,3) - P(1,4)*P(4,3)*Metric(2,3) + P(-1,1)*P(-1,2)*Metric(1,4)*Metric(2,3) + P(-1,1)*P(-1,3)*Metric(1,4)*Metric(2,3) + P(-1,2)*P(-1,4)*Metric(1,4)*Metric(2,3) + P(-1,3)*P(-1,4)*Metric(1,4)*Metric(2,3) - P(1,2)*P(3,1)*Metric(2,4) - P(1,4)*P(3,1)*Metric(2,4) - P(1,3)*P(3,2)*Metric(2,4) + P(1,4)*P(3,2)*Metric(2,4) + P(1,2)*P(3,4)*Metric(2,4) - P(1,3)*P(3,4)*Metric(2,4) + P(-1,1)*P(-1,2)*Metric(1,3)*Metric(2,4) + P(-1,2)*P(-1,3)*Metric(1,3)*Metric(2,4) + P(-1,1)*P(-1,4)*Metric(1,3)*Metric(2,4) + P(-1,3)*P(-1,4)*Metric(1,3)*Metric(2,4) - P(1,3)*P(2,1)*Metric(3,4) - P(1,4)*P(2,1)*Metric(3,4) - P(1,2)*P(2,3)*Metric(3,4) + P(1,4)*P(2,3)*Metric(3,4) - P(1,2)*P(2,4)*Metric(3,4) + P(1,3)*P(2,4)*Metric(3,4) + P(-1,1)*P(-1,3)*Metric(1,2)*Metric(3,4) + P(-1,2)*P(-1,3)*Metric(1,2)*Metric(3,4) + P(-1,1)*P(-1,4)*Metric(1,2)*Metric(3,4) + P(-1,2)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVVSS14 = Lorentz(name = 'VVVVSS14',
                   spins = [ 3, 3, 3, 3, 1, 1 ],
                   structure = 'P(3,4)*P(4,3)*Metric(1,2) - P(-1,3)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVVSS15 = Lorentz(name = 'VVVVSS15',
                   spins = [ 3, 3, 3, 3, 1, 1 ],
                   structure = 'P(3,4)*P(4,3)*Metric(1,2) + P(1,2)*P(2,1)*Metric(3,4) - P(-1,1)*P(-1,2)*Metric(1,2)*Metric(3,4) - P(-1,3)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVVSS16 = Lorentz(name = 'VVVVSS16',
                   spins = [ 3, 3, 3, 3, 1, 1 ],
                   structure = 'P(3,4)*P(4,3)*Metric(1,2) + P(2,4)*P(4,2)*Metric(1,3) + P(2,3)*P(3,2)*Metric(1,4) + P(1,4)*P(4,1)*Metric(2,3) - P(-1,2)*P(-1,3)*Metric(1,4)*Metric(2,3) - P(-1,1)*P(-1,4)*Metric(1,4)*Metric(2,3) + P(1,3)*P(3,1)*Metric(2,4) - P(-1,1)*P(-1,3)*Metric(1,3)*Metric(2,4) - P(-1,2)*P(-1,4)*Metric(1,3)*Metric(2,4) + P(1,2)*P(2,1)*Metric(3,4) - P(-1,1)*P(-1,2)*Metric(1,2)*Metric(3,4) - P(-1,3)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVVSS17 = Lorentz(name = 'VVVVSS17',
                   spins = [ 3, 3, 3, 3, 1, 1 ],
                   structure = 'P(-1,5)*P(-1,6)*Metric(1,4)*Metric(2,3) + P(-1,5)*P(-1,6)*Metric(1,3)*Metric(2,4) - 2*P(-1,5)*P(-1,6)*Metric(1,2)*Metric(3,4)')

VVVVSS18 = Lorentz(name = 'VVVVSS18',
                   spins = [ 3, 3, 3, 3, 1, 1 ],
                   structure = 'P(-1,5)*P(-1,6)*Metric(1,4)*Metric(2,3) - (P(-1,5)*P(-1,6)*Metric(1,3)*Metric(2,4))/2. - (P(-1,5)*P(-1,6)*Metric(1,2)*Metric(3,4))/2.')

VVVVVS1 = Lorentz(name = 'VVVVVS1',
                  spins = [ 3, 3, 3, 3, 3, 1 ],
                  structure = 'P(5,2)*Metric(1,4)*Metric(2,3) - P(5,3)*Metric(1,4)*Metric(2,3) + P(4,2)*Metric(1,5)*Metric(2,3) - P(4,3)*Metric(1,5)*Metric(2,3) - P(5,1)*Metric(1,3)*Metric(2,4) + P(5,3)*Metric(1,3)*Metric(2,4) + P(3,1)*Metric(1,5)*Metric(2,4) - P(3,2)*Metric(1,5)*Metric(2,4) - P(4,1)*Metric(1,3)*Metric(2,5) + P(4,3)*Metric(1,3)*Metric(2,5) + P(3,1)*Metric(1,4)*Metric(2,5) - P(3,2)*Metric(1,4)*Metric(2,5) + P(5,1)*Metric(1,2)*Metric(3,4) - P(5,2)*Metric(1,2)*Metric(3,4) - P(2,1)*Metric(1,5)*Metric(3,4) + P(2,3)*Metric(1,5)*Metric(3,4) + P(1,2)*Metric(2,5)*Metric(3,4) - P(1,3)*Metric(2,5)*Metric(3,4) + P(4,1)*Metric(1,2)*Metric(3,5) - P(4,2)*Metric(1,2)*Metric(3,5) - P(2,1)*Metric(1,4)*Metric(3,5) + P(2,3)*Metric(1,4)*Metric(3,5) + P(1,2)*Metric(2,4)*Metric(3,5) - P(1,3)*Metric(2,4)*Metric(3,5)')

VVVVVS2 = Lorentz(name = 'VVVVVS2',
                  spins = [ 3, 3, 3, 3, 3, 1 ],
                  structure = 'P(5,1)*Metric(1,3)*Metric(2,4) - P(5,3)*Metric(1,3)*Metric(2,4) - P(3,1)*Metric(1,5)*Metric(2,4) + P(3,5)*Metric(1,5)*Metric(2,4) + P(4,1)*Metric(1,3)*Metric(2,5) - P(4,3)*Metric(1,3)*Metric(2,5) - P(3,1)*Metric(1,4)*Metric(2,5) + P(3,4)*Metric(1,4)*Metric(2,5) + P(5,1)*Metric(1,2)*Metric(3,4) - P(5,2)*Metric(1,2)*Metric(3,4) - P(2,1)*Metric(1,5)*Metric(3,4) + P(2,5)*Metric(1,5)*Metric(3,4) + P(1,2)*Metric(2,5)*Metric(3,4) + P(1,3)*Metric(2,5)*Metric(3,4) - P(1,4)*Metric(2,5)*Metric(3,4) - P(1,5)*Metric(2,5)*Metric(3,4) + P(4,1)*Metric(1,2)*Metric(3,5) - P(4,2)*Metric(1,2)*Metric(3,5) - P(2,1)*Metric(1,4)*Metric(3,5) + P(2,4)*Metric(1,4)*Metric(3,5) + P(1,2)*Metric(2,4)*Metric(3,5) + P(1,3)*Metric(2,4)*Metric(3,5) - P(1,4)*Metric(2,4)*Metric(3,5) - P(1,5)*Metric(2,4)*Metric(3,5)')

VVVVVS3 = Lorentz(name = 'VVVVVS3',
                  spins = [ 3, 3, 3, 3, 3, 1 ],
                  structure = 'P(3,1)*Metric(1,2)*Metric(4,5) - P(3,2)*Metric(1,2)*Metric(4,5) - P(2,1)*Metric(1,3)*Metric(4,5) + P(2,3)*Metric(1,3)*Metric(4,5) + P(1,2)*Metric(2,3)*Metric(4,5) - P(1,3)*Metric(2,3)*Metric(4,5)')

VVVVVS4 = Lorentz(name = 'VVVVVS4',
                  spins = [ 3, 3, 3, 3, 3, 1 ],
                  structure = 'P(5,1)*Metric(1,2)*Metric(3,4) - P(5,2)*Metric(1,2)*Metric(3,4) - P(2,1)*Metric(1,5)*Metric(3,4) + P(2,5)*Metric(1,5)*Metric(3,4) + P(1,2)*Metric(2,5)*Metric(3,4) - P(1,5)*Metric(2,5)*Metric(3,4) + P(4,1)*Metric(1,2)*Metric(3,5) - P(4,2)*Metric(1,2)*Metric(3,5) - P(2,1)*Metric(1,4)*Metric(3,5) + P(2,4)*Metric(1,4)*Metric(3,5) + P(1,2)*Metric(2,4)*Metric(3,5) - P(1,4)*Metric(2,4)*Metric(3,5) + P(3,1)*Metric(1,2)*Metric(4,5) - P(3,2)*Metric(1,2)*Metric(4,5) - P(2,1)*Metric(1,3)*Metric(4,5) + P(2,3)*Metric(1,3)*Metric(4,5) + P(1,2)*Metric(2,3)*Metric(4,5) - P(1,3)*Metric(2,3)*Metric(4,5)')

VVVVVS5 = Lorentz(name = 'VVVVVS5',
                  spins = [ 3, 3, 3, 3, 3, 1 ],
                  structure = 'P(5,1)*Metric(1,4)*Metric(2,3) + P(5,2)*Metric(1,4)*Metric(2,3) - P(5,3)*Metric(1,4)*Metric(2,3) - P(5,4)*Metric(1,4)*Metric(2,3) - P(4,1)*Metric(1,5)*Metric(2,3) - P(4,3)*Metric(1,5)*Metric(2,3) + 2*P(4,5)*Metric(1,5)*Metric(2,3) + P(5,1)*Metric(1,3)*Metric(2,4) + P(5,2)*Metric(1,3)*Metric(2,4) - P(5,3)*Metric(1,3)*Metric(2,4) - P(5,4)*Metric(1,3)*Metric(2,4) - P(3,1)*Metric(1,5)*Metric(2,4) - P(3,4)*Metric(1,5)*Metric(2,4) + 2*P(3,5)*Metric(1,5)*Metric(2,4) - P(4,2)*Metric(1,3)*Metric(2,5) - P(4,3)*Metric(1,3)*Metric(2,5) + 2*P(4,5)*Metric(1,3)*Metric(2,5) - P(3,2)*Metric(1,4)*Metric(2,5) - P(3,4)*Metric(1,4)*Metric(2,5) + 2*P(3,5)*Metric(1,4)*Metric(2,5) + 2*P(5,1)*Metric(1,2)*Metric(3,4) + 2*P(5,2)*Metric(1,2)*Metric(3,4) - 2*P(5,3)*Metric(1,2)*Metric(3,4) - 2*P(5,4)*Metric(1,2)*Metric(3,4) - 2*P(2,1)*Metric(1,5)*Metric(3,4) + P(2,3)*Metric(1,5)*Metric(3,4) + P(2,4)*Metric(1,5)*Metric(3,4) - 2*P(1,2)*Metric(2,5)*Metric(3,4) + P(1,3)*Metric(2,5)*Metric(3,4) + P(1,4)*Metric(2,5)*Metric(3,4) - P(4,1)*Metric(1,2)*Metric(3,5) - P(4,2)*Metric(1,2)*Metric(3,5) + 2*P(4,3)*Metric(1,2)*Metric(3,5) + P(2,1)*Metric(1,4)*Metric(3,5) + P(2,3)*Metric(1,4)*Metric(3,5) - 2*P(2,5)*Metric(1,4)*Metric(3,5) + P(1,2)*Metric(2,4)*Metric(3,5) + P(1,3)*Metric(2,4)*Metric(3,5) - 2*P(1,5)*Metric(2,4)*Metric(3,5) - P(3,1)*Metric(1,2)*Metric(4,5) - P(3,2)*Metric(1,2)*Metric(4,5) + 2*P(3,4)*Metric(1,2)*Metric(4,5) + P(2,1)*Metric(1,3)*Metric(4,5) + P(2,4)*Metric(1,3)*Metric(4,5) - 2*P(2,5)*Metric(1,3)*Metric(4,5) + P(1,2)*Metric(2,3)*Metric(4,5) + P(1,4)*Metric(2,3)*Metric(4,5) - 2*P(1,5)*Metric(2,3)*Metric(4,5)')

VVVVVS6 = Lorentz(name = 'VVVVVS6',
                  spins = [ 3, 3, 3, 3, 3, 1 ],
                  structure = 'P(5,1)*Metric(1,4)*Metric(2,3) + P(5,2)*Metric(1,4)*Metric(2,3) - P(5,3)*Metric(1,4)*Metric(2,3) - P(5,4)*Metric(1,4)*Metric(2,3) - P(4,1)*Metric(1,5)*Metric(2,3) + P(4,5)*Metric(1,5)*Metric(2,3) + P(5,1)*Metric(1,3)*Metric(2,4) + P(5,2)*Metric(1,3)*Metric(2,4) - P(5,3)*Metric(1,3)*Metric(2,4) - P(5,4)*Metric(1,3)*Metric(2,4) - P(3,1)*Metric(1,5)*Metric(2,4) + P(3,5)*Metric(1,5)*Metric(2,4) - P(4,2)*Metric(1,3)*Metric(2,5) + P(4,5)*Metric(1,3)*Metric(2,5) - P(3,2)*Metric(1,4)*Metric(2,5) + P(3,5)*Metric(1,4)*Metric(2,5) + P(2,3)*Metric(1,4)*Metric(3,5) - P(2,5)*Metric(1,4)*Metric(3,5) + P(1,3)*Metric(2,4)*Metric(3,5) - P(1,5)*Metric(2,4)*Metric(3,5) + P(2,4)*Metric(1,3)*Metric(4,5) - P(2,5)*Metric(1,3)*Metric(4,5) + P(1,4)*Metric(2,3)*Metric(4,5) - P(1,5)*Metric(2,3)*Metric(4,5)')

VVVVVS7 = Lorentz(name = 'VVVVVS7',
                  spins = [ 3, 3, 3, 3, 3, 1 ],
                  structure = 'P(5,3)*Metric(1,4)*Metric(2,3) - P(5,4)*Metric(1,4)*Metric(2,3) + P(4,3)*Metric(1,5)*Metric(2,3) - P(4,5)*Metric(1,5)*Metric(2,3) - P(5,3)*Metric(1,3)*Metric(2,4) + P(5,4)*Metric(1,3)*Metric(2,4) + P(3,4)*Metric(1,5)*Metric(2,4) - P(3,5)*Metric(1,5)*Metric(2,4) - P(4,3)*Metric(1,3)*Metric(2,5) + P(4,5)*Metric(1,3)*Metric(2,5) - P(3,4)*Metric(1,4)*Metric(2,5) + P(3,5)*Metric(1,4)*Metric(2,5) - 2*P(5,1)*Metric(1,2)*Metric(3,4) + 2*P(5,2)*Metric(1,2)*Metric(3,4) + 2*P(2,1)*Metric(1,5)*Metric(3,4) - P(2,3)*Metric(1,5)*Metric(3,4) - P(2,4)*Metric(1,5)*Metric(3,4) - 2*P(1,2)*Metric(2,5)*Metric(3,4) + P(1,3)*Metric(2,5)*Metric(3,4) + P(1,4)*Metric(2,5)*Metric(3,4) - 2*P(4,1)*Metric(1,2)*Metric(3,5) + 2*P(4,2)*Metric(1,2)*Metric(3,5) + 2*P(2,1)*Metric(1,4)*Metric(3,5) - P(2,3)*Metric(1,4)*Metric(3,5) - P(2,5)*Metric(1,4)*Metric(3,5) - 2*P(1,2)*Metric(2,4)*Metric(3,5) + P(1,3)*Metric(2,4)*Metric(3,5) + P(1,5)*Metric(2,4)*Metric(3,5) - 2*P(3,1)*Metric(1,2)*Metric(4,5) + 2*P(3,2)*Metric(1,2)*Metric(4,5) + 2*P(2,1)*Metric(1,3)*Metric(4,5) - P(2,4)*Metric(1,3)*Metric(4,5) - P(2,5)*Metric(1,3)*Metric(4,5) - 2*P(1,2)*Metric(2,3)*Metric(4,5) + P(1,4)*Metric(2,3)*Metric(4,5) + P(1,5)*Metric(2,3)*Metric(4,5)')

VVVVVS8 = Lorentz(name = 'VVVVVS8',
                  spins = [ 3, 3, 3, 3, 3, 1 ],
                  structure = 'P(5,2)*Metric(1,4)*Metric(2,3) + P(5,3)*Metric(1,4)*Metric(2,3) - 2*P(5,4)*Metric(1,4)*Metric(2,3) + P(4,2)*Metric(1,5)*Metric(2,3) + P(4,3)*Metric(1,5)*Metric(2,3) - 2*P(4,5)*Metric(1,5)*Metric(2,3) - 2*P(5,1)*Metric(1,3)*Metric(2,4) + P(5,3)*Metric(1,3)*Metric(2,4) + P(5,4)*Metric(1,3)*Metric(2,4) + 2*P(3,1)*Metric(1,5)*Metric(2,4) - P(3,2)*Metric(1,5)*Metric(2,4) - P(3,5)*Metric(1,5)*Metric(2,4) - 2*P(4,1)*Metric(1,3)*Metric(2,5) + P(4,3)*Metric(1,3)*Metric(2,5) + P(4,5)*Metric(1,3)*Metric(2,5) + 2*P(3,1)*Metric(1,4)*Metric(2,5) - P(3,2)*Metric(1,4)*Metric(2,5) - P(3,4)*Metric(1,4)*Metric(2,5) - 2*P(5,1)*Metric(1,2)*Metric(3,4) + P(5,2)*Metric(1,2)*Metric(3,4) + P(5,4)*Metric(1,2)*Metric(3,4) + 2*P(2,1)*Metric(1,5)*Metric(3,4) - P(2,3)*Metric(1,5)*Metric(3,4) - P(2,5)*Metric(1,5)*Metric(3,4) - P(1,2)*Metric(2,5)*Metric(3,4) - P(1,3)*Metric(2,5)*Metric(3,4) + P(1,4)*Metric(2,5)*Metric(3,4) + P(1,5)*Metric(2,5)*Metric(3,4) - 2*P(4,1)*Metric(1,2)*Metric(3,5) + P(4,2)*Metric(1,2)*Metric(3,5) + P(4,5)*Metric(1,2)*Metric(3,5) + 2*P(2,1)*Metric(1,4)*Metric(3,5) - P(2,3)*Metric(1,4)*Metric(3,5) - P(2,4)*Metric(1,4)*Metric(3,5) - P(1,2)*Metric(2,4)*Metric(3,5) - P(1,3)*Metric(2,4)*Metric(3,5) + P(1,4)*Metric(2,4)*Metric(3,5) + P(1,5)*Metric(2,4)*Metric(3,5) + 2*P(3,2)*Metric(1,2)*Metric(4,5) - P(3,4)*Metric(1,2)*Metric(4,5) - P(3,5)*Metric(1,2)*Metric(4,5) + 2*P(2,3)*Metric(1,3)*Metric(4,5) - P(2,4)*Metric(1,3)*Metric(4,5) - P(2,5)*Metric(1,3)*Metric(4,5) - 2*P(1,2)*Metric(2,3)*Metric(4,5) - 2*P(1,3)*Metric(2,3)*Metric(4,5) + 2*P(1,4)*Metric(2,3)*Metric(4,5) + 2*P(1,5)*Metric(2,3)*Metric(4,5)')

VVVVVV1 = Lorentz(name = 'VVVVVV1',
                  spins = [ 3, 3, 3, 3, 3, 3 ],
                  structure = 'Metric(1,6)*Metric(2,4)*Metric(3,5) + Metric(1,4)*Metric(2,6)*Metric(3,5) + Metric(1,5)*Metric(2,4)*Metric(3,6) + Metric(1,4)*Metric(2,5)*Metric(3,6) + Metric(1,6)*Metric(2,3)*Metric(4,5) + Metric(1,3)*Metric(2,6)*Metric(4,5) - 4*Metric(1,2)*Metric(3,6)*Metric(4,5) + Metric(1,5)*Metric(2,3)*Metric(4,6) + Metric(1,3)*Metric(2,5)*Metric(4,6) - 4*Metric(1,2)*Metric(3,5)*Metric(4,6)')

VVVVVV2 = Lorentz(name = 'VVVVVV2',
                  spins = [ 3, 3, 3, 3, 3, 3 ],
                  structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) + Metric(1,5)*Metric(2,6)*Metric(3,4) - (Metric(1,6)*Metric(2,4)*Metric(3,5))/2. - (Metric(1,4)*Metric(2,6)*Metric(3,5))/2. - (Metric(1,5)*Metric(2,4)*Metric(3,6))/2. - (Metric(1,4)*Metric(2,5)*Metric(3,6))/2. - (Metric(1,6)*Metric(2,3)*Metric(4,5))/2. - (Metric(1,3)*Metric(2,6)*Metric(4,5))/2. + Metric(1,2)*Metric(3,6)*Metric(4,5) - (Metric(1,5)*Metric(2,3)*Metric(4,6))/2. - (Metric(1,3)*Metric(2,5)*Metric(4,6))/2. + Metric(1,2)*Metric(3,5)*Metric(4,6)')

VVVVVV3 = Lorentz(name = 'VVVVVV3',
                  spins = [ 3, 3, 3, 3, 3, 3 ],
                  structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) + Metric(1,6)*Metric(2,4)*Metric(3,5) + 2*Metric(1,6)*Metric(2,3)*Metric(4,5) - Metric(1,3)*Metric(2,6)*Metric(4,5) - Metric(1,2)*Metric(3,6)*Metric(4,5) - Metric(1,5)*Metric(2,3)*Metric(4,6) - Metric(1,4)*Metric(2,3)*Metric(5,6)')

VVVVVV4 = Lorentz(name = 'VVVVVV4',
                  spins = [ 3, 3, 3, 3, 3, 3 ],
                  structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) + Metric(1,5)*Metric(2,6)*Metric(3,4) + Metric(1,6)*Metric(2,4)*Metric(3,5) + Metric(1,4)*Metric(2,6)*Metric(3,5) + Metric(1,5)*Metric(2,4)*Metric(3,6) + Metric(1,4)*Metric(2,5)*Metric(3,6) + Metric(1,6)*Metric(2,3)*Metric(4,5) + Metric(1,3)*Metric(2,6)*Metric(4,5) - 4*Metric(1,2)*Metric(3,6)*Metric(4,5) + Metric(1,5)*Metric(2,3)*Metric(4,6) + Metric(1,3)*Metric(2,5)*Metric(4,6) - 4*Metric(1,2)*Metric(3,5)*Metric(4,6) + Metric(1,4)*Metric(2,3)*Metric(5,6) + Metric(1,3)*Metric(2,4)*Metric(5,6) - 4*Metric(1,2)*Metric(3,4)*Metric(5,6)')

VVVVVV5 = Lorentz(name = 'VVVVVV5',
                  spins = [ 3, 3, 3, 3, 3, 3 ],
                  structure = 'Metric(1,4)*Metric(2,3)*Metric(5,6) + Metric(1,3)*Metric(2,4)*Metric(5,6) - 2*Metric(1,2)*Metric(3,4)*Metric(5,6)')

VVVVVV6 = Lorentz(name = 'VVVVVV6',
                  spins = [ 3, 3, 3, 3, 3, 3 ],
                  structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) + Metric(1,5)*Metric(2,6)*Metric(3,4) - Metric(1,2)*Metric(3,6)*Metric(4,5) - Metric(1,2)*Metric(3,5)*Metric(4,6) + Metric(1,4)*Metric(2,3)*Metric(5,6) + Metric(1,3)*Metric(2,4)*Metric(5,6) - 2*Metric(1,2)*Metric(3,4)*Metric(5,6)')

VVVVVV7 = Lorentz(name = 'VVVVVV7',
                  spins = [ 3, 3, 3, 3, 3, 3 ],
                  structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) + Metric(1,5)*Metric(2,6)*Metric(3,4) - (Metric(1,6)*Metric(2,4)*Metric(3,5))/6. - (Metric(1,4)*Metric(2,6)*Metric(3,5))/6. - (Metric(1,5)*Metric(2,4)*Metric(3,6))/6. - (Metric(1,4)*Metric(2,5)*Metric(3,6))/6. - (Metric(1,6)*Metric(2,3)*Metric(4,5))/6. - (Metric(1,3)*Metric(2,6)*Metric(4,5))/6. + Metric(1,2)*Metric(3,6)*Metric(4,5) - (Metric(1,5)*Metric(2,3)*Metric(4,6))/6. - (Metric(1,3)*Metric(2,5)*Metric(4,6))/6. + Metric(1,2)*Metric(3,5)*Metric(4,6) - (2*Metric(1,4)*Metric(2,3)*Metric(5,6))/3. - (2*Metric(1,3)*Metric(2,4)*Metric(5,6))/3. - (4*Metric(1,2)*Metric(3,4)*Metric(5,6))/3.')

VVVVVV8 = Lorentz(name = 'VVVVVV8',
                  spins = [ 3, 3, 3, 3, 3, 3 ],
                  structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) + Metric(1,5)*Metric(2,6)*Metric(3,4) + Metric(1,6)*Metric(2,4)*Metric(3,5) + Metric(1,4)*Metric(2,6)*Metric(3,5) + Metric(1,5)*Metric(2,4)*Metric(3,6) + Metric(1,4)*Metric(2,5)*Metric(3,6) - (2*Metric(1,6)*Metric(2,3)*Metric(4,5))/3. - (2*Metric(1,3)*Metric(2,6)*Metric(4,5))/3. - (2*Metric(1,2)*Metric(3,6)*Metric(4,5))/3. - (2*Metric(1,5)*Metric(2,3)*Metric(4,6))/3. - (2*Metric(1,3)*Metric(2,5)*Metric(4,6))/3. - (2*Metric(1,2)*Metric(3,5)*Metric(4,6))/3. - (2*Metric(1,4)*Metric(2,3)*Metric(5,6))/3. - (2*Metric(1,3)*Metric(2,4)*Metric(5,6))/3. - (2*Metric(1,2)*Metric(3,4)*Metric(5,6))/3.')

VVVVVV9 = Lorentz(name = 'VVVVVV9',
                  spins = [ 3, 3, 3, 3, 3, 3 ],
                  structure = 'Metric(1,6)*Metric(2,3)*Metric(4,5) - (Metric(1,3)*Metric(2,6)*Metric(4,5))/2. - (Metric(1,2)*Metric(3,6)*Metric(4,5))/2. + Metric(1,5)*Metric(2,3)*Metric(4,6) - (Metric(1,3)*Metric(2,5)*Metric(4,6))/2. - (Metric(1,2)*Metric(3,5)*Metric(4,6))/2. + Metric(1,4)*Metric(2,3)*Metric(5,6) - (Metric(1,3)*Metric(2,4)*Metric(5,6))/2. - (Metric(1,2)*Metric(3,4)*Metric(5,6))/2.')

VVVVVV10 = Lorentz(name = 'VVVVVV10',
                   spins = [ 3, 3, 3, 3, 3, 3 ],
                   structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) - (Metric(1,5)*Metric(2,6)*Metric(3,4))/4. + Metric(1,6)*Metric(2,4)*Metric(3,5) - (Metric(1,4)*Metric(2,6)*Metric(3,5))/4. - (Metric(1,5)*Metric(2,4)*Metric(3,6))/4. - (Metric(1,4)*Metric(2,5)*Metric(3,6))/4. - (Metric(1,3)*Metric(2,5)*Metric(4,6))/4. - (Metric(1,2)*Metric(3,5)*Metric(4,6))/4. - (Metric(1,3)*Metric(2,4)*Metric(5,6))/4. - (Metric(1,2)*Metric(3,4)*Metric(5,6))/4.')

VVVVVV11 = Lorentz(name = 'VVVVVV11',
                   spins = [ 3, 3, 3, 3, 3, 3 ],
                   structure = 'Metric(1,6)*Metric(2,4)*Metric(3,5) + Metric(1,4)*Metric(2,6)*Metric(3,5) + Metric(1,5)*Metric(2,4)*Metric(3,6) + Metric(1,4)*Metric(2,5)*Metric(3,6) + Metric(1,6)*Metric(2,3)*Metric(4,5) + Metric(1,3)*Metric(2,6)*Metric(4,5) + Metric(1,5)*Metric(2,3)*Metric(4,6) + Metric(1,3)*Metric(2,5)*Metric(4,6) - 5*Metric(1,4)*Metric(2,3)*Metric(5,6) - 5*Metric(1,3)*Metric(2,4)*Metric(5,6) + 2*Metric(1,2)*Metric(3,4)*Metric(5,6)')

VVVVVSS1 = Lorentz(name = 'VVVVVSS1',
                   spins = [ 3, 3, 3, 3, 3, 1, 1 ],
                   structure = 'P(5,2)*Metric(1,4)*Metric(2,3) - P(5,3)*Metric(1,4)*Metric(2,3) + P(4,2)*Metric(1,5)*Metric(2,3) - P(4,3)*Metric(1,5)*Metric(2,3) - P(5,1)*Metric(1,3)*Metric(2,4) + P(5,3)*Metric(1,3)*Metric(2,4) + P(3,1)*Metric(1,5)*Metric(2,4) - P(3,2)*Metric(1,5)*Metric(2,4) - P(4,1)*Metric(1,3)*Metric(2,5) + P(4,3)*Metric(1,3)*Metric(2,5) + P(3,1)*Metric(1,4)*Metric(2,5) - P(3,2)*Metric(1,4)*Metric(2,5) + P(5,1)*Metric(1,2)*Metric(3,4) - P(5,2)*Metric(1,2)*Metric(3,4) - P(2,1)*Metric(1,5)*Metric(3,4) + P(2,3)*Metric(1,5)*Metric(3,4) + P(1,2)*Metric(2,5)*Metric(3,4) - P(1,3)*Metric(2,5)*Metric(3,4) + P(4,1)*Metric(1,2)*Metric(3,5) - P(4,2)*Metric(1,2)*Metric(3,5) - P(2,1)*Metric(1,4)*Metric(3,5) + P(2,3)*Metric(1,4)*Metric(3,5) + P(1,2)*Metric(2,4)*Metric(3,5) - P(1,3)*Metric(2,4)*Metric(3,5)')

VVVVVSS2 = Lorentz(name = 'VVVVVSS2',
                   spins = [ 3, 3, 3, 3, 3, 1, 1 ],
                   structure = 'P(5,1)*Metric(1,3)*Metric(2,4) - P(5,3)*Metric(1,3)*Metric(2,4) - P(3,1)*Metric(1,5)*Metric(2,4) + P(3,5)*Metric(1,5)*Metric(2,4) + P(4,1)*Metric(1,3)*Metric(2,5) - P(4,3)*Metric(1,3)*Metric(2,5) - P(3,1)*Metric(1,4)*Metric(2,5) + P(3,4)*Metric(1,4)*Metric(2,5) + P(5,1)*Metric(1,2)*Metric(3,4) - P(5,2)*Metric(1,2)*Metric(3,4) - P(2,1)*Metric(1,5)*Metric(3,4) + P(2,5)*Metric(1,5)*Metric(3,4) + P(1,2)*Metric(2,5)*Metric(3,4) + P(1,3)*Metric(2,5)*Metric(3,4) - P(1,4)*Metric(2,5)*Metric(3,4) - P(1,5)*Metric(2,5)*Metric(3,4) + P(4,1)*Metric(1,2)*Metric(3,5) - P(4,2)*Metric(1,2)*Metric(3,5) - P(2,1)*Metric(1,4)*Metric(3,5) + P(2,4)*Metric(1,4)*Metric(3,5) + P(1,2)*Metric(2,4)*Metric(3,5) + P(1,3)*Metric(2,4)*Metric(3,5) - P(1,4)*Metric(2,4)*Metric(3,5) - P(1,5)*Metric(2,4)*Metric(3,5)')

VVVVVSS3 = Lorentz(name = 'VVVVVSS3',
                   spins = [ 3, 3, 3, 3, 3, 1, 1 ],
                   structure = 'P(3,1)*Metric(1,2)*Metric(4,5) - P(3,2)*Metric(1,2)*Metric(4,5) - P(2,1)*Metric(1,3)*Metric(4,5) + P(2,3)*Metric(1,3)*Metric(4,5) + P(1,2)*Metric(2,3)*Metric(4,5) - P(1,3)*Metric(2,3)*Metric(4,5)')

VVVVVSS4 = Lorentz(name = 'VVVVVSS4',
                   spins = [ 3, 3, 3, 3, 3, 1, 1 ],
                   structure = 'P(5,1)*Metric(1,2)*Metric(3,4) - P(5,2)*Metric(1,2)*Metric(3,4) - P(2,1)*Metric(1,5)*Metric(3,4) + P(2,5)*Metric(1,5)*Metric(3,4) + P(1,2)*Metric(2,5)*Metric(3,4) - P(1,5)*Metric(2,5)*Metric(3,4) + P(4,1)*Metric(1,2)*Metric(3,5) - P(4,2)*Metric(1,2)*Metric(3,5) - P(2,1)*Metric(1,4)*Metric(3,5) + P(2,4)*Metric(1,4)*Metric(3,5) + P(1,2)*Metric(2,4)*Metric(3,5) - P(1,4)*Metric(2,4)*Metric(3,5) + P(3,1)*Metric(1,2)*Metric(4,5) - P(3,2)*Metric(1,2)*Metric(4,5) - P(2,1)*Metric(1,3)*Metric(4,5) + P(2,3)*Metric(1,3)*Metric(4,5) + P(1,2)*Metric(2,3)*Metric(4,5) - P(1,3)*Metric(2,3)*Metric(4,5)')

VVVVVSS5 = Lorentz(name = 'VVVVVSS5',
                   spins = [ 3, 3, 3, 3, 3, 1, 1 ],
                   structure = 'P(5,1)*Metric(1,4)*Metric(2,3) + P(5,2)*Metric(1,4)*Metric(2,3) - P(5,3)*Metric(1,4)*Metric(2,3) - P(5,4)*Metric(1,4)*Metric(2,3) - P(4,1)*Metric(1,5)*Metric(2,3) - P(4,3)*Metric(1,5)*Metric(2,3) + 2*P(4,5)*Metric(1,5)*Metric(2,3) + P(5,1)*Metric(1,3)*Metric(2,4) + P(5,2)*Metric(1,3)*Metric(2,4) - P(5,3)*Metric(1,3)*Metric(2,4) - P(5,4)*Metric(1,3)*Metric(2,4) - P(3,1)*Metric(1,5)*Metric(2,4) - P(3,4)*Metric(1,5)*Metric(2,4) + 2*P(3,5)*Metric(1,5)*Metric(2,4) - P(4,2)*Metric(1,3)*Metric(2,5) - P(4,3)*Metric(1,3)*Metric(2,5) + 2*P(4,5)*Metric(1,3)*Metric(2,5) - P(3,2)*Metric(1,4)*Metric(2,5) - P(3,4)*Metric(1,4)*Metric(2,5) + 2*P(3,5)*Metric(1,4)*Metric(2,5) + 2*P(5,1)*Metric(1,2)*Metric(3,4) + 2*P(5,2)*Metric(1,2)*Metric(3,4) - 2*P(5,3)*Metric(1,2)*Metric(3,4) - 2*P(5,4)*Metric(1,2)*Metric(3,4) - 2*P(2,1)*Metric(1,5)*Metric(3,4) + P(2,3)*Metric(1,5)*Metric(3,4) + P(2,4)*Metric(1,5)*Metric(3,4) - 2*P(1,2)*Metric(2,5)*Metric(3,4) + P(1,3)*Metric(2,5)*Metric(3,4) + P(1,4)*Metric(2,5)*Metric(3,4) - P(4,1)*Metric(1,2)*Metric(3,5) - P(4,2)*Metric(1,2)*Metric(3,5) + 2*P(4,3)*Metric(1,2)*Metric(3,5) + P(2,1)*Metric(1,4)*Metric(3,5) + P(2,3)*Metric(1,4)*Metric(3,5) - 2*P(2,5)*Metric(1,4)*Metric(3,5) + P(1,2)*Metric(2,4)*Metric(3,5) + P(1,3)*Metric(2,4)*Metric(3,5) - 2*P(1,5)*Metric(2,4)*Metric(3,5) - P(3,1)*Metric(1,2)*Metric(4,5) - P(3,2)*Metric(1,2)*Metric(4,5) + 2*P(3,4)*Metric(1,2)*Metric(4,5) + P(2,1)*Metric(1,3)*Metric(4,5) + P(2,4)*Metric(1,3)*Metric(4,5) - 2*P(2,5)*Metric(1,3)*Metric(4,5) + P(1,2)*Metric(2,3)*Metric(4,5) + P(1,4)*Metric(2,3)*Metric(4,5) - 2*P(1,5)*Metric(2,3)*Metric(4,5)')

VVVVVSS6 = Lorentz(name = 'VVVVVSS6',
                   spins = [ 3, 3, 3, 3, 3, 1, 1 ],
                   structure = 'P(5,1)*Metric(1,4)*Metric(2,3) + P(5,2)*Metric(1,4)*Metric(2,3) - P(5,3)*Metric(1,4)*Metric(2,3) - P(5,4)*Metric(1,4)*Metric(2,3) - P(4,1)*Metric(1,5)*Metric(2,3) + P(4,5)*Metric(1,5)*Metric(2,3) + P(5,1)*Metric(1,3)*Metric(2,4) + P(5,2)*Metric(1,3)*Metric(2,4) - P(5,3)*Metric(1,3)*Metric(2,4) - P(5,4)*Metric(1,3)*Metric(2,4) - P(3,1)*Metric(1,5)*Metric(2,4) + P(3,5)*Metric(1,5)*Metric(2,4) - P(4,2)*Metric(1,3)*Metric(2,5) + P(4,5)*Metric(1,3)*Metric(2,5) - P(3,2)*Metric(1,4)*Metric(2,5) + P(3,5)*Metric(1,4)*Metric(2,5) + P(2,3)*Metric(1,4)*Metric(3,5) - P(2,5)*Metric(1,4)*Metric(3,5) + P(1,3)*Metric(2,4)*Metric(3,5) - P(1,5)*Metric(2,4)*Metric(3,5) + P(2,4)*Metric(1,3)*Metric(4,5) - P(2,5)*Metric(1,3)*Metric(4,5) + P(1,4)*Metric(2,3)*Metric(4,5) - P(1,5)*Metric(2,3)*Metric(4,5)')

VVVVVSS7 = Lorentz(name = 'VVVVVSS7',
                   spins = [ 3, 3, 3, 3, 3, 1, 1 ],
                   structure = 'P(5,3)*Metric(1,4)*Metric(2,3) - P(5,4)*Metric(1,4)*Metric(2,3) + P(4,3)*Metric(1,5)*Metric(2,3) - P(4,5)*Metric(1,5)*Metric(2,3) - P(5,3)*Metric(1,3)*Metric(2,4) + P(5,4)*Metric(1,3)*Metric(2,4) + P(3,4)*Metric(1,5)*Metric(2,4) - P(3,5)*Metric(1,5)*Metric(2,4) - P(4,3)*Metric(1,3)*Metric(2,5) + P(4,5)*Metric(1,3)*Metric(2,5) - P(3,4)*Metric(1,4)*Metric(2,5) + P(3,5)*Metric(1,4)*Metric(2,5) - 2*P(5,1)*Metric(1,2)*Metric(3,4) + 2*P(5,2)*Metric(1,2)*Metric(3,4) + 2*P(2,1)*Metric(1,5)*Metric(3,4) - P(2,3)*Metric(1,5)*Metric(3,4) - P(2,4)*Metric(1,5)*Metric(3,4) - 2*P(1,2)*Metric(2,5)*Metric(3,4) + P(1,3)*Metric(2,5)*Metric(3,4) + P(1,4)*Metric(2,5)*Metric(3,4) - 2*P(4,1)*Metric(1,2)*Metric(3,5) + 2*P(4,2)*Metric(1,2)*Metric(3,5) + 2*P(2,1)*Metric(1,4)*Metric(3,5) - P(2,3)*Metric(1,4)*Metric(3,5) - P(2,5)*Metric(1,4)*Metric(3,5) - 2*P(1,2)*Metric(2,4)*Metric(3,5) + P(1,3)*Metric(2,4)*Metric(3,5) + P(1,5)*Metric(2,4)*Metric(3,5) - 2*P(3,1)*Metric(1,2)*Metric(4,5) + 2*P(3,2)*Metric(1,2)*Metric(4,5) + 2*P(2,1)*Metric(1,3)*Metric(4,5) - P(2,4)*Metric(1,3)*Metric(4,5) - P(2,5)*Metric(1,3)*Metric(4,5) - 2*P(1,2)*Metric(2,3)*Metric(4,5) + P(1,4)*Metric(2,3)*Metric(4,5) + P(1,5)*Metric(2,3)*Metric(4,5)')

VVVVVSS8 = Lorentz(name = 'VVVVVSS8',
                   spins = [ 3, 3, 3, 3, 3, 1, 1 ],
                   structure = 'P(5,2)*Metric(1,4)*Metric(2,3) + P(5,3)*Metric(1,4)*Metric(2,3) - 2*P(5,4)*Metric(1,4)*Metric(2,3) + P(4,2)*Metric(1,5)*Metric(2,3) + P(4,3)*Metric(1,5)*Metric(2,3) - 2*P(4,5)*Metric(1,5)*Metric(2,3) - 2*P(5,1)*Metric(1,3)*Metric(2,4) + P(5,3)*Metric(1,3)*Metric(2,4) + P(5,4)*Metric(1,3)*Metric(2,4) + 2*P(3,1)*Metric(1,5)*Metric(2,4) - P(3,2)*Metric(1,5)*Metric(2,4) - P(3,5)*Metric(1,5)*Metric(2,4) - 2*P(4,1)*Metric(1,3)*Metric(2,5) + P(4,3)*Metric(1,3)*Metric(2,5) + P(4,5)*Metric(1,3)*Metric(2,5) + 2*P(3,1)*Metric(1,4)*Metric(2,5) - P(3,2)*Metric(1,4)*Metric(2,5) - P(3,4)*Metric(1,4)*Metric(2,5) - 2*P(5,1)*Metric(1,2)*Metric(3,4) + P(5,2)*Metric(1,2)*Metric(3,4) + P(5,4)*Metric(1,2)*Metric(3,4) + 2*P(2,1)*Metric(1,5)*Metric(3,4) - P(2,3)*Metric(1,5)*Metric(3,4) - P(2,5)*Metric(1,5)*Metric(3,4) - P(1,2)*Metric(2,5)*Metric(3,4) - P(1,3)*Metric(2,5)*Metric(3,4) + P(1,4)*Metric(2,5)*Metric(3,4) + P(1,5)*Metric(2,5)*Metric(3,4) - 2*P(4,1)*Metric(1,2)*Metric(3,5) + P(4,2)*Metric(1,2)*Metric(3,5) + P(4,5)*Metric(1,2)*Metric(3,5) + 2*P(2,1)*Metric(1,4)*Metric(3,5) - P(2,3)*Metric(1,4)*Metric(3,5) - P(2,4)*Metric(1,4)*Metric(3,5) - P(1,2)*Metric(2,4)*Metric(3,5) - P(1,3)*Metric(2,4)*Metric(3,5) + P(1,4)*Metric(2,4)*Metric(3,5) + P(1,5)*Metric(2,4)*Metric(3,5) + 2*P(3,2)*Metric(1,2)*Metric(4,5) - P(3,4)*Metric(1,2)*Metric(4,5) - P(3,5)*Metric(1,2)*Metric(4,5) + 2*P(2,3)*Metric(1,3)*Metric(4,5) - P(2,4)*Metric(1,3)*Metric(4,5) - P(2,5)*Metric(1,3)*Metric(4,5) - 2*P(1,2)*Metric(2,3)*Metric(4,5) - 2*P(1,3)*Metric(2,3)*Metric(4,5) + 2*P(1,4)*Metric(2,3)*Metric(4,5) + 2*P(1,5)*Metric(2,3)*Metric(4,5)')

VVVVVVS1 = Lorentz(name = 'VVVVVVS1',
                   spins = [ 3, 3, 3, 3, 3, 3, 1 ],
                   structure = 'Metric(1,6)*Metric(2,4)*Metric(3,5) + Metric(1,4)*Metric(2,6)*Metric(3,5) + Metric(1,5)*Metric(2,4)*Metric(3,6) + Metric(1,4)*Metric(2,5)*Metric(3,6) + Metric(1,6)*Metric(2,3)*Metric(4,5) + Metric(1,3)*Metric(2,6)*Metric(4,5) - 4*Metric(1,2)*Metric(3,6)*Metric(4,5) + Metric(1,5)*Metric(2,3)*Metric(4,6) + Metric(1,3)*Metric(2,5)*Metric(4,6) - 4*Metric(1,2)*Metric(3,5)*Metric(4,6)')

VVVVVVS2 = Lorentz(name = 'VVVVVVS2',
                   spins = [ 3, 3, 3, 3, 3, 3, 1 ],
                   structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) + Metric(1,5)*Metric(2,6)*Metric(3,4) - (Metric(1,6)*Metric(2,4)*Metric(3,5))/2. - (Metric(1,4)*Metric(2,6)*Metric(3,5))/2. - (Metric(1,5)*Metric(2,4)*Metric(3,6))/2. - (Metric(1,4)*Metric(2,5)*Metric(3,6))/2. - (Metric(1,6)*Metric(2,3)*Metric(4,5))/2. - (Metric(1,3)*Metric(2,6)*Metric(4,5))/2. + Metric(1,2)*Metric(3,6)*Metric(4,5) - (Metric(1,5)*Metric(2,3)*Metric(4,6))/2. - (Metric(1,3)*Metric(2,5)*Metric(4,6))/2. + Metric(1,2)*Metric(3,5)*Metric(4,6)')

VVVVVVS3 = Lorentz(name = 'VVVVVVS3',
                   spins = [ 3, 3, 3, 3, 3, 3, 1 ],
                   structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) + Metric(1,6)*Metric(2,4)*Metric(3,5) + 2*Metric(1,6)*Metric(2,3)*Metric(4,5) - Metric(1,3)*Metric(2,6)*Metric(4,5) - Metric(1,2)*Metric(3,6)*Metric(4,5) - Metric(1,5)*Metric(2,3)*Metric(4,6) - Metric(1,4)*Metric(2,3)*Metric(5,6)')

VVVVVVS4 = Lorentz(name = 'VVVVVVS4',
                   spins = [ 3, 3, 3, 3, 3, 3, 1 ],
                   structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) + Metric(1,5)*Metric(2,6)*Metric(3,4) + Metric(1,6)*Metric(2,4)*Metric(3,5) + Metric(1,4)*Metric(2,6)*Metric(3,5) + Metric(1,5)*Metric(2,4)*Metric(3,6) + Metric(1,4)*Metric(2,5)*Metric(3,6) + Metric(1,6)*Metric(2,3)*Metric(4,5) + Metric(1,3)*Metric(2,6)*Metric(4,5) - 4*Metric(1,2)*Metric(3,6)*Metric(4,5) + Metric(1,5)*Metric(2,3)*Metric(4,6) + Metric(1,3)*Metric(2,5)*Metric(4,6) - 4*Metric(1,2)*Metric(3,5)*Metric(4,6) + Metric(1,4)*Metric(2,3)*Metric(5,6) + Metric(1,3)*Metric(2,4)*Metric(5,6) - 4*Metric(1,2)*Metric(3,4)*Metric(5,6)')

VVVVVVS5 = Lorentz(name = 'VVVVVVS5',
                   spins = [ 3, 3, 3, 3, 3, 3, 1 ],
                   structure = 'Metric(1,4)*Metric(2,3)*Metric(5,6) + Metric(1,3)*Metric(2,4)*Metric(5,6) - 2*Metric(1,2)*Metric(3,4)*Metric(5,6)')

VVVVVVS6 = Lorentz(name = 'VVVVVVS6',
                   spins = [ 3, 3, 3, 3, 3, 3, 1 ],
                   structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) + Metric(1,5)*Metric(2,6)*Metric(3,4) - Metric(1,2)*Metric(3,6)*Metric(4,5) - Metric(1,2)*Metric(3,5)*Metric(4,6) + Metric(1,4)*Metric(2,3)*Metric(5,6) + Metric(1,3)*Metric(2,4)*Metric(5,6) - 2*Metric(1,2)*Metric(3,4)*Metric(5,6)')

VVVVVVS7 = Lorentz(name = 'VVVVVVS7',
                   spins = [ 3, 3, 3, 3, 3, 3, 1 ],
                   structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) + Metric(1,5)*Metric(2,6)*Metric(3,4) - (Metric(1,6)*Metric(2,4)*Metric(3,5))/6. - (Metric(1,4)*Metric(2,6)*Metric(3,5))/6. - (Metric(1,5)*Metric(2,4)*Metric(3,6))/6. - (Metric(1,4)*Metric(2,5)*Metric(3,6))/6. - (Metric(1,6)*Metric(2,3)*Metric(4,5))/6. - (Metric(1,3)*Metric(2,6)*Metric(4,5))/6. + Metric(1,2)*Metric(3,6)*Metric(4,5) - (Metric(1,5)*Metric(2,3)*Metric(4,6))/6. - (Metric(1,3)*Metric(2,5)*Metric(4,6))/6. + Metric(1,2)*Metric(3,5)*Metric(4,6) - (2*Metric(1,4)*Metric(2,3)*Metric(5,6))/3. - (2*Metric(1,3)*Metric(2,4)*Metric(5,6))/3. - (4*Metric(1,2)*Metric(3,4)*Metric(5,6))/3.')

VVVVVVS8 = Lorentz(name = 'VVVVVVS8',
                   spins = [ 3, 3, 3, 3, 3, 3, 1 ],
                   structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) + Metric(1,5)*Metric(2,6)*Metric(3,4) + Metric(1,6)*Metric(2,4)*Metric(3,5) + Metric(1,4)*Metric(2,6)*Metric(3,5) + Metric(1,5)*Metric(2,4)*Metric(3,6) + Metric(1,4)*Metric(2,5)*Metric(3,6) - (2*Metric(1,6)*Metric(2,3)*Metric(4,5))/3. - (2*Metric(1,3)*Metric(2,6)*Metric(4,5))/3. - (2*Metric(1,2)*Metric(3,6)*Metric(4,5))/3. - (2*Metric(1,5)*Metric(2,3)*Metric(4,6))/3. - (2*Metric(1,3)*Metric(2,5)*Metric(4,6))/3. - (2*Metric(1,2)*Metric(3,5)*Metric(4,6))/3. - (2*Metric(1,4)*Metric(2,3)*Metric(5,6))/3. - (2*Metric(1,3)*Metric(2,4)*Metric(5,6))/3. - (2*Metric(1,2)*Metric(3,4)*Metric(5,6))/3.')

VVVVVVS9 = Lorentz(name = 'VVVVVVS9',
                   spins = [ 3, 3, 3, 3, 3, 3, 1 ],
                   structure = 'Metric(1,6)*Metric(2,3)*Metric(4,5) - (Metric(1,3)*Metric(2,6)*Metric(4,5))/2. - (Metric(1,2)*Metric(3,6)*Metric(4,5))/2. + Metric(1,5)*Metric(2,3)*Metric(4,6) - (Metric(1,3)*Metric(2,5)*Metric(4,6))/2. - (Metric(1,2)*Metric(3,5)*Metric(4,6))/2. + Metric(1,4)*Metric(2,3)*Metric(5,6) - (Metric(1,3)*Metric(2,4)*Metric(5,6))/2. - (Metric(1,2)*Metric(3,4)*Metric(5,6))/2.')

VVVVVVS10 = Lorentz(name = 'VVVVVVS10',
                    spins = [ 3, 3, 3, 3, 3, 3, 1 ],
                    structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) - (Metric(1,5)*Metric(2,6)*Metric(3,4))/4. + Metric(1,6)*Metric(2,4)*Metric(3,5) - (Metric(1,4)*Metric(2,6)*Metric(3,5))/4. - (Metric(1,5)*Metric(2,4)*Metric(3,6))/4. - (Metric(1,4)*Metric(2,5)*Metric(3,6))/4. - (Metric(1,3)*Metric(2,5)*Metric(4,6))/4. - (Metric(1,2)*Metric(3,5)*Metric(4,6))/4. - (Metric(1,3)*Metric(2,4)*Metric(5,6))/4. - (Metric(1,2)*Metric(3,4)*Metric(5,6))/4.')

VVVVVVS11 = Lorentz(name = 'VVVVVVS11',
                    spins = [ 3, 3, 3, 3, 3, 3, 1 ],
                    structure = 'Metric(1,6)*Metric(2,4)*Metric(3,5) + Metric(1,4)*Metric(2,6)*Metric(3,5) + Metric(1,5)*Metric(2,4)*Metric(3,6) + Metric(1,4)*Metric(2,5)*Metric(3,6) + Metric(1,6)*Metric(2,3)*Metric(4,5) + Metric(1,3)*Metric(2,6)*Metric(4,5) + Metric(1,5)*Metric(2,3)*Metric(4,6) + Metric(1,3)*Metric(2,5)*Metric(4,6) - 5*Metric(1,4)*Metric(2,3)*Metric(5,6) - 5*Metric(1,3)*Metric(2,4)*Metric(5,6) + 2*Metric(1,2)*Metric(3,4)*Metric(5,6)')

VVVVVVSS1 = Lorentz(name = 'VVVVVVSS1',
                    spins = [ 3, 3, 3, 3, 3, 3, 1, 1 ],
                    structure = 'Metric(1,6)*Metric(2,4)*Metric(3,5) + Metric(1,4)*Metric(2,6)*Metric(3,5) + Metric(1,5)*Metric(2,4)*Metric(3,6) + Metric(1,4)*Metric(2,5)*Metric(3,6) + Metric(1,6)*Metric(2,3)*Metric(4,5) + Metric(1,3)*Metric(2,6)*Metric(4,5) - 4*Metric(1,2)*Metric(3,6)*Metric(4,5) + Metric(1,5)*Metric(2,3)*Metric(4,6) + Metric(1,3)*Metric(2,5)*Metric(4,6) - 4*Metric(1,2)*Metric(3,5)*Metric(4,6)')

VVVVVVSS2 = Lorentz(name = 'VVVVVVSS2',
                    spins = [ 3, 3, 3, 3, 3, 3, 1, 1 ],
                    structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) + Metric(1,5)*Metric(2,6)*Metric(3,4) - (Metric(1,6)*Metric(2,4)*Metric(3,5))/2. - (Metric(1,4)*Metric(2,6)*Metric(3,5))/2. - (Metric(1,5)*Metric(2,4)*Metric(3,6))/2. - (Metric(1,4)*Metric(2,5)*Metric(3,6))/2. - (Metric(1,6)*Metric(2,3)*Metric(4,5))/2. - (Metric(1,3)*Metric(2,6)*Metric(4,5))/2. + Metric(1,2)*Metric(3,6)*Metric(4,5) - (Metric(1,5)*Metric(2,3)*Metric(4,6))/2. - (Metric(1,3)*Metric(2,5)*Metric(4,6))/2. + Metric(1,2)*Metric(3,5)*Metric(4,6)')

VVVVVVSS3 = Lorentz(name = 'VVVVVVSS3',
                    spins = [ 3, 3, 3, 3, 3, 3, 1, 1 ],
                    structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) + Metric(1,6)*Metric(2,4)*Metric(3,5) + 2*Metric(1,6)*Metric(2,3)*Metric(4,5) - Metric(1,3)*Metric(2,6)*Metric(4,5) - Metric(1,2)*Metric(3,6)*Metric(4,5) - Metric(1,5)*Metric(2,3)*Metric(4,6) - Metric(1,4)*Metric(2,3)*Metric(5,6)')

VVVVVVSS4 = Lorentz(name = 'VVVVVVSS4',
                    spins = [ 3, 3, 3, 3, 3, 3, 1, 1 ],
                    structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) + Metric(1,5)*Metric(2,6)*Metric(3,4) + Metric(1,6)*Metric(2,4)*Metric(3,5) + Metric(1,4)*Metric(2,6)*Metric(3,5) + Metric(1,5)*Metric(2,4)*Metric(3,6) + Metric(1,4)*Metric(2,5)*Metric(3,6) + Metric(1,6)*Metric(2,3)*Metric(4,5) + Metric(1,3)*Metric(2,6)*Metric(4,5) - 4*Metric(1,2)*Metric(3,6)*Metric(4,5) + Metric(1,5)*Metric(2,3)*Metric(4,6) + Metric(1,3)*Metric(2,5)*Metric(4,6) - 4*Metric(1,2)*Metric(3,5)*Metric(4,6) + Metric(1,4)*Metric(2,3)*Metric(5,6) + Metric(1,3)*Metric(2,4)*Metric(5,6) - 4*Metric(1,2)*Metric(3,4)*Metric(5,6)')

VVVVVVSS5 = Lorentz(name = 'VVVVVVSS5',
                    spins = [ 3, 3, 3, 3, 3, 3, 1, 1 ],
                    structure = 'Metric(1,4)*Metric(2,3)*Metric(5,6) + Metric(1,3)*Metric(2,4)*Metric(5,6) - 2*Metric(1,2)*Metric(3,4)*Metric(5,6)')

VVVVVVSS6 = Lorentz(name = 'VVVVVVSS6',
                    spins = [ 3, 3, 3, 3, 3, 3, 1, 1 ],
                    structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) + Metric(1,5)*Metric(2,6)*Metric(3,4) - Metric(1,2)*Metric(3,6)*Metric(4,5) - Metric(1,2)*Metric(3,5)*Metric(4,6) + Metric(1,4)*Metric(2,3)*Metric(5,6) + Metric(1,3)*Metric(2,4)*Metric(5,6) - 2*Metric(1,2)*Metric(3,4)*Metric(5,6)')

VVVVVVSS7 = Lorentz(name = 'VVVVVVSS7',
                    spins = [ 3, 3, 3, 3, 3, 3, 1, 1 ],
                    structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) + Metric(1,5)*Metric(2,6)*Metric(3,4) - (Metric(1,6)*Metric(2,4)*Metric(3,5))/6. - (Metric(1,4)*Metric(2,6)*Metric(3,5))/6. - (Metric(1,5)*Metric(2,4)*Metric(3,6))/6. - (Metric(1,4)*Metric(2,5)*Metric(3,6))/6. - (Metric(1,6)*Metric(2,3)*Metric(4,5))/6. - (Metric(1,3)*Metric(2,6)*Metric(4,5))/6. + Metric(1,2)*Metric(3,6)*Metric(4,5) - (Metric(1,5)*Metric(2,3)*Metric(4,6))/6. - (Metric(1,3)*Metric(2,5)*Metric(4,6))/6. + Metric(1,2)*Metric(3,5)*Metric(4,6) - (2*Metric(1,4)*Metric(2,3)*Metric(5,6))/3. - (2*Metric(1,3)*Metric(2,4)*Metric(5,6))/3. - (4*Metric(1,2)*Metric(3,4)*Metric(5,6))/3.')

VVVVVVSS8 = Lorentz(name = 'VVVVVVSS8',
                    spins = [ 3, 3, 3, 3, 3, 3, 1, 1 ],
                    structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) + Metric(1,5)*Metric(2,6)*Metric(3,4) + Metric(1,6)*Metric(2,4)*Metric(3,5) + Metric(1,4)*Metric(2,6)*Metric(3,5) + Metric(1,5)*Metric(2,4)*Metric(3,6) + Metric(1,4)*Metric(2,5)*Metric(3,6) - (2*Metric(1,6)*Metric(2,3)*Metric(4,5))/3. - (2*Metric(1,3)*Metric(2,6)*Metric(4,5))/3. - (2*Metric(1,2)*Metric(3,6)*Metric(4,5))/3. - (2*Metric(1,5)*Metric(2,3)*Metric(4,6))/3. - (2*Metric(1,3)*Metric(2,5)*Metric(4,6))/3. - (2*Metric(1,2)*Metric(3,5)*Metric(4,6))/3. - (2*Metric(1,4)*Metric(2,3)*Metric(5,6))/3. - (2*Metric(1,3)*Metric(2,4)*Metric(5,6))/3. - (2*Metric(1,2)*Metric(3,4)*Metric(5,6))/3.')

VVVVVVSS9 = Lorentz(name = 'VVVVVVSS9',
                    spins = [ 3, 3, 3, 3, 3, 3, 1, 1 ],
                    structure = 'Metric(1,6)*Metric(2,3)*Metric(4,5) - (Metric(1,3)*Metric(2,6)*Metric(4,5))/2. - (Metric(1,2)*Metric(3,6)*Metric(4,5))/2. + Metric(1,5)*Metric(2,3)*Metric(4,6) - (Metric(1,3)*Metric(2,5)*Metric(4,6))/2. - (Metric(1,2)*Metric(3,5)*Metric(4,6))/2. + Metric(1,4)*Metric(2,3)*Metric(5,6) - (Metric(1,3)*Metric(2,4)*Metric(5,6))/2. - (Metric(1,2)*Metric(3,4)*Metric(5,6))/2.')

VVVVVVSS10 = Lorentz(name = 'VVVVVVSS10',
                     spins = [ 3, 3, 3, 3, 3, 3, 1, 1 ],
                     structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) - (Metric(1,5)*Metric(2,6)*Metric(3,4))/4. + Metric(1,6)*Metric(2,4)*Metric(3,5) - (Metric(1,4)*Metric(2,6)*Metric(3,5))/4. - (Metric(1,5)*Metric(2,4)*Metric(3,6))/4. - (Metric(1,4)*Metric(2,5)*Metric(3,6))/4. - (Metric(1,3)*Metric(2,5)*Metric(4,6))/4. - (Metric(1,2)*Metric(3,5)*Metric(4,6))/4. - (Metric(1,3)*Metric(2,4)*Metric(5,6))/4. - (Metric(1,2)*Metric(3,4)*Metric(5,6))/4.')

VVVVVVSS11 = Lorentz(name = 'VVVVVVSS11',
                     spins = [ 3, 3, 3, 3, 3, 3, 1, 1 ],
                     structure = 'Metric(1,6)*Metric(2,4)*Metric(3,5) + Metric(1,4)*Metric(2,6)*Metric(3,5) + Metric(1,5)*Metric(2,4)*Metric(3,6) + Metric(1,4)*Metric(2,5)*Metric(3,6) + Metric(1,6)*Metric(2,3)*Metric(4,5) + Metric(1,3)*Metric(2,6)*Metric(4,5) + Metric(1,5)*Metric(2,3)*Metric(4,6) + Metric(1,3)*Metric(2,5)*Metric(4,6) - 5*Metric(1,4)*Metric(2,3)*Metric(5,6) - 5*Metric(1,3)*Metric(2,4)*Metric(5,6) + 2*Metric(1,2)*Metric(3,4)*Metric(5,6)')

