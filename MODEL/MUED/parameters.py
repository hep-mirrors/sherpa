# This file was automatically created by FeynRules $Revision: 915 $
# Mathematica version: 7.0 for Mac OS X x86 (64-bit) (November 11, 2008)
# Date: Mon 23 Jan 2012 11:06:17



from object_library import all_parameters, Parameter


from function_library import complexconjugate, re, im, csc, sec, acsc, asec

# This is a default parameter object representing 0.
ZERO = Parameter(name = 'ZERO',
                 nature = 'internal',
                 type = 'real',
                 value = '0.0',
                 texname = '0')

# User-defined parameters.
cabi = Parameter(name = 'cabi',
                 nature = 'external',
                 type = 'real',
                 value = 0.488,
                 texname = '\\theta _c',
                 lhablock = 'CKMBLOCK',
                 lhacode = [ 1 ])

s12 = Parameter(name = 's12',
                nature = 'external',
                type = 'real',
                value = 0.221,
                texname = '\\text{s12}',
                lhablock = 'CKMINPUTS',
                lhacode = [ 1 ])

s23 = Parameter(name = 's23',
                nature = 'external',
                type = 'real',
                value = 0.041,
                texname = '\\text{s23}',
                lhablock = 'CKMINPUTS',
                lhacode = [ 2 ])

s13 = Parameter(name = 's13',
                nature = 'external',
                type = 'real',
                value = 0.0035,
                texname = '\\text{s13}',
                lhablock = 'CKMINPUTS',
                lhacode = [ 3 ])

amass = Parameter(name = 'amass',
                  nature = 'external',
                  type = 'real',
                  value = 0.1172,
                  texname = '\\text{amass}',
                  lhablock = 'MASSCORR',
                  lhacode = [ 1 ])

R = Parameter(name = 'R',
              nature = 'external',
              type = 'real',
              value = 0.002,
              texname = 'R',
              lhablock = 'MUEDINPUTS',
              lhacode = [ 1 ])

LR = Parameter(name = 'LR',
               nature = 'external',
               type = 'real',
               value = 20,
               texname = '\\text{LR}',
               lhablock = 'MUEDINPUTS',
               lhacode = [ 2 ])

aEWM1 = Parameter(name = 'aEWM1',
                  nature = 'external',
                  type = 'real',
                  value = 128,
                  texname = '\\text{aEWM1}',
                  lhablock = 'SMINPUTS',
                  lhacode = [ 1 ])

Gf = Parameter(name = 'Gf',
               nature = 'external',
               type = 'real',
               value = 0.000011663900000000002,
               texname = '\\text{Gf}',
               lhablock = 'SMINPUTS',
               lhacode = [ 2 ])

aS = Parameter(name = 'aS',
               nature = 'external',
               type = 'real',
               value = 0.1172,
               texname = '\\text{aS}',
               lhablock = 'SMINPUTS',
               lhacode = [ 3 ])

ZM = Parameter(name = 'ZM',
               nature = 'external',
               type = 'real',
               value = 91.1876,
               texname = '\\text{ZM}',
               lhablock = 'SMINPUTS',
               lhacode = [ 4 ])

sw = Parameter(name = 'sw',
               nature = 'external',
               type = 'real',
               value = 0.48076,
               texname = 's_w',
               lhablock = 'SMINPUTS',
               lhacode = [ 5 ])

ymc = Parameter(name = 'ymc',
                nature = 'external',
                type = 'real',
                value = 1.42,
                texname = '\\text{ymc}',
                lhablock = 'YUKAWA',
                lhacode = [ 4 ])

ymb = Parameter(name = 'ymb',
                nature = 'external',
                type = 'real',
                value = 4.2,
                texname = '\\text{ymb}',
                lhablock = 'YUKAWA',
                lhacode = [ 5 ])

ymt = Parameter(name = 'ymt',
                nature = 'external',
                type = 'real',
                value = 175,
                texname = '\\text{ymt}',
                lhablock = 'YUKAWA',
                lhacode = [ 6 ])

ymtau = Parameter(name = 'ymtau',
                  nature = 'external',
                  type = 'real',
                  value = 1.777,
                  texname = '\\text{ymtau}',
                  lhablock = 'YUKAWA',
                  lhacode = [ 15 ])

y1 = Parameter(name = 'y1',
               nature = 'external',
               type = 'real',
               value = -1,
               texname = '\\text{y1}',
               lhablock = 'YVALUES',
               lhacode = [ 1 ])

y2 = Parameter(name = 'y2',
               nature = 'external',
               type = 'real',
               value = -2,
               texname = '\\text{y2}',
               lhablock = 'YVALUES',
               lhacode = [ 2 ])

y3 = Parameter(name = 'y3',
               nature = 'external',
               type = 'real',
               value = 0.3333333333333333,
               texname = '\\text{y3}',
               lhablock = 'YVALUES',
               lhacode = [ 3 ])

y4 = Parameter(name = 'y4',
               nature = 'external',
               type = 'real',
               value = 1.3333333333333333,
               texname = '\\text{y4}',
               lhablock = 'YVALUES',
               lhacode = [ 4 ])

y5 = Parameter(name = 'y5',
               nature = 'external',
               type = 'real',
               value = -0.6666666666666666,
               texname = '\\text{y5}',
               lhablock = 'YVALUES',
               lhacode = [ 5 ])

MM = Parameter(name = 'MM',
               nature = 'external',
               type = 'real',
               value = 0.1057,
               texname = '\\text{MM}',
               lhablock = 'MASS',
               lhacode = [ 13 ])

MTA = Parameter(name = 'MTA',
                nature = 'external',
                type = 'real',
                value = 1.77,
                texname = '\\text{MTA}',
                lhablock = 'MASS',
                lhacode = [ 15 ])

MU = Parameter(name = 'MU',
               nature = 'external',
               type = 'real',
               value = 0.003,
               texname = 'M',
               lhablock = 'MASS',
               lhacode = [ 2 ])

MC = Parameter(name = 'MC',
               nature = 'external',
               type = 'real',
               value = 1.42,
               texname = '\\text{MC}',
               lhablock = 'MASS',
               lhacode = [ 4 ])

MT = Parameter(name = 'MT',
               nature = 'external',
               type = 'real',
               value = 175,
               texname = '\\text{MT}',
               lhablock = 'MASS',
               lhacode = [ 6 ])

MD = Parameter(name = 'MD',
               nature = 'external',
               type = 'real',
               value = 0.007,
               texname = '\\text{MD}',
               lhablock = 'MASS',
               lhacode = [ 1 ])

MS = Parameter(name = 'MS',
               nature = 'external',
               type = 'real',
               value = 0.2,
               texname = '\\text{MS}',
               lhablock = 'MASS',
               lhacode = [ 3 ])

MB = Parameter(name = 'MB',
               nature = 'external',
               type = 'real',
               value = 4.2,
               texname = '\\text{MB}',
               lhablock = 'MASS',
               lhacode = [ 5 ])

MZ = Parameter(name = 'MZ',
               nature = 'external',
               type = 'real',
               value = 91.1876,
               texname = '\\text{MZ}',
               lhablock = 'MASS',
               lhacode = [ 23 ])

MH = Parameter(name = 'MH',
               nature = 'external',
               type = 'real',
               value = 120,
               texname = '\\text{MH}',
               lhablock = 'MASS',
               lhacode = [ 25 ])

WT = Parameter(name = 'WT',
               nature = 'external',
               type = 'real',
               value = 1.59,
               texname = '\\text{WT}',
               lhablock = 'DECAY',
               lhacode = [ 6 ])

WZ = Parameter(name = 'WZ',
               nature = 'external',
               type = 'real',
               value = 2.49444,
               texname = '\\text{WZ}',
               lhablock = 'DECAY',
               lhacode = [ 23 ])

WW = Parameter(name = 'WW',
               nature = 'external',
               type = 'real',
               value = 2.08895,
               texname = '\\text{WW}',
               lhablock = 'DECAY',
               lhacode = [ 24 ])

WH = Parameter(name = 'WH',
               nature = 'external',
               type = 'real',
               value = 0.004244,
               texname = '\\text{WH}',
               lhablock = 'DECAY',
               lhacode = [ 25 ])

DWe = Parameter(name = 'DWe',
                nature = 'external',
                type = 'real',
                value = 0.001463,
                texname = '\\text{DWe}',
                lhablock = 'DECAY',
                lhacode = [ 5100011 ])

DWmu = Parameter(name = 'DWmu',
                 nature = 'external',
                 type = 'real',
                 value = 0.001463,
                 texname = '\\text{DWmu}',
                 lhablock = 'DECAY',
                 lhacode = [ 5100013 ])

DWtau = Parameter(name = 'DWtau',
                  nature = 'external',
                  type = 'real',
                  value = 0.001463,
                  texname = '\\text{DWtau}',
                  lhablock = 'DECAY',
                  lhacode = [ 5100015 ])

DWn1 = Parameter(name = 'DWn1',
                 nature = 'external',
                 type = 'real',
                 value = 0.001463,
                 texname = '\\text{DWn1}',
                 lhablock = 'DECAY',
                 lhacode = [ 5100012 ])

DWn2 = Parameter(name = 'DWn2',
                 nature = 'external',
                 type = 'real',
                 value = 0.001463,
                 texname = '\\text{DWn2}',
                 lhablock = 'DECAY',
                 lhacode = [ 5100014 ])

DWn3 = Parameter(name = 'DWn3',
                 nature = 'external',
                 type = 'real',
                 value = 0.001463,
                 texname = '\\text{DWn3}',
                 lhablock = 'DECAY',
                 lhacode = [ 5100016 ])

SWe = Parameter(name = 'SWe',
                nature = 'external',
                type = 'real',
                value = 0.00062291,
                texname = '\\text{SWe}',
                lhablock = 'DECAY',
                lhacode = [ 6100011 ])

SWmu = Parameter(name = 'SWmu',
                 nature = 'external',
                 type = 'real',
                 value = 0.000062291,
                 texname = '\\text{SWmu}',
                 lhablock = 'DECAY',
                 lhacode = [ 6100013 ])

SWtau = Parameter(name = 'SWtau',
                  nature = 'external',
                  type = 'real',
                  value = 0.00062291,
                  texname = '\\text{SWtau}',
                  lhablock = 'DECAY',
                  lhacode = [ 6100015 ])

SWu = Parameter(name = 'SWu',
                nature = 'external',
                type = 'real',
                value = 0.082179,
                texname = '\\text{SWu}',
                lhablock = 'DECAY',
                lhacode = [ 6100002 ])

SWc = Parameter(name = 'SWc',
                nature = 'external',
                type = 'real',
                value = 0.082179,
                texname = '\\text{SWc}',
                lhablock = 'DECAY',
                lhacode = [ 6100004 ])

SWd = Parameter(name = 'SWd',
                nature = 'external',
                type = 'real',
                value = 0.019815,
                texname = '\\text{SWd}',
                lhablock = 'DECAY',
                lhacode = [ 6100001 ])

SWs = Parameter(name = 'SWs',
                nature = 'external',
                type = 'real',
                value = 0.019754,
                texname = '\\text{SWs}',
                lhablock = 'DECAY',
                lhacode = [ 6100003 ])

SWb = Parameter(name = 'SWb',
                nature = 'external',
                type = 'real',
                value = 0.019761,
                texname = '\\text{SWb}',
                lhablock = 'DECAY',
                lhacode = [ 6100005 ])

DWu = Parameter(name = 'DWu',
                nature = 'external',
                type = 'real',
                value = 0.2512,
                texname = '\\text{DWu}',
                lhablock = 'DECAY',
                lhacode = [ 5100002 ])

DWc = Parameter(name = 'DWc',
                nature = 'external',
                type = 'real',
                value = 0.2512,
                texname = '\\text{DWc}',
                lhablock = 'DECAY',
                lhacode = [ 5100004 ])

DWtop = Parameter(name = 'DWtop',
                  nature = 'external',
                  type = 'real',
                  value = 0.10096,
                  texname = '\\text{DWtop}',
                  lhablock = 'DECAY',
                  lhacode = [ 5100006 ])

DWd = Parameter(name = 'DWd',
                nature = 'external',
                type = 'real',
                value = 0.2512,
                texname = '\\text{DWd}',
                lhablock = 'DECAY',
                lhacode = [ 5100001 ])

DWs = Parameter(name = 'DWs',
                nature = 'external',
                type = 'real',
                value = 0.25101,
                texname = '\\text{DWs}',
                lhablock = 'DECAY',
                lhacode = [ 5100003 ])

DWb = Parameter(name = 'DWb',
                nature = 'external',
                type = 'real',
                value = 0.054041,
                texname = '\\text{DWb}',
                lhablock = 'DECAY',
                lhacode = [ 5100005 ])

WZ1 = Parameter(name = 'WZ1',
                nature = 'external',
                type = 'real',
                value = 0.0756786,
                texname = '\\text{WZ1}',
                lhablock = 'DECAY',
                lhacode = [ 5100023 ])

WW1 = Parameter(name = 'WW1',
                nature = 'external',
                type = 'real',
                value = 0.07345485,
                texname = '\\text{WW1}',
                lhablock = 'DECAY',
                lhacode = [ 5100024 ])

WG1 = Parameter(name = 'WG1',
                nature = 'external',
                type = 'real',
                value = 4.207153,
                texname = '\\text{WG1}',
                lhablock = 'DECAY',
                lhacode = [ 5100021 ])

WH1 = Parameter(name = 'WH1',
                nature = 'external',
                type = 'real',
                value = 0.004244,
                texname = '\\text{WH1}',
                lhablock = 'DECAY',
                lhacode = [ 5100030 ])

aEW = Parameter(name = 'aEW',
                nature = 'internal',
                type = 'real',
                value = '1./aEWM1',
                texname = '\\text{aEW}')

sw2 = Parameter(name = 'sw2',
                nature = 'internal',
                type = 'real',
                value = 'sw**2',
                texname = '\\text{sw2}')

G = Parameter(name = 'G',
              nature = 'internal',
              type = 'real',
              value = '2*cmath.sqrt(aS)*cmath.sqrt(cmath.pi)',
              texname = 'G')

c12 = Parameter(name = 'c12',
                nature = 'internal',
                type = 'real',
                value = 'cmath.sqrt(1.0 - s12**2)',
                texname = '\\text{c12}')

c23 = Parameter(name = 'c23',
                nature = 'internal',
                type = 'real',
                value = 'cmath.sqrt(1.0 - s23**2)',
                texname = '\\text{c23}')

c13 = Parameter(name = 'c13',
                nature = 'internal',
                type = 'real',
                value = 'cmath.sqrt(1.0 - s13**2)',
                texname = '\\text{c13}')

CKM13 = Parameter(name = 'CKM13',
                  nature = 'internal',
                  type = 'complex',
                  value = 's13',
                  texname = '\\text{CKM13}')

YukD11 = Parameter(name = 'YukD11',
                   nature = 'internal',
                   type = 'complex',
                   value = '0',
                   texname = '\\text{YukD11}')

YukD12 = Parameter(name = 'YukD12',
                   nature = 'internal',
                   type = 'complex',
                   value = '0',
                   texname = '\\text{YukD12}')

YukD13 = Parameter(name = 'YukD13',
                   nature = 'internal',
                   type = 'complex',
                   value = '0',
                   texname = '\\text{YukD13}')

YukD21 = Parameter(name = 'YukD21',
                   nature = 'internal',
                   type = 'complex',
                   value = '0',
                   texname = '\\text{YukD21}')

YukD22 = Parameter(name = 'YukD22',
                   nature = 'internal',
                   type = 'complex',
                   value = '0',
                   texname = '\\text{YukD22}')

YukD23 = Parameter(name = 'YukD23',
                   nature = 'internal',
                   type = 'complex',
                   value = '0',
                   texname = '\\text{YukD23}')

YukD31 = Parameter(name = 'YukD31',
                   nature = 'internal',
                   type = 'complex',
                   value = '0',
                   texname = '\\text{YukD31}')

YukD32 = Parameter(name = 'YukD32',
                   nature = 'internal',
                   type = 'complex',
                   value = '0',
                   texname = '\\text{YukD32}')

gmass = Parameter(name = 'gmass',
                  nature = 'internal',
                  type = 'real',
                  value = '2.0*cmath.sqrt(amass)*cmath.sqrt(cmath.pi)',
                  texname = '\\text{gmass}')

mn = Parameter(name = 'mn',
               nature = 'internal',
               type = 'real',
               value = 'cmath.log(LR**2)/(16.*cmath.pi**2*R)',
               texname = '\\text{mn}')

CKM11 = Parameter(name = 'CKM11',
                  nature = 'internal',
                  type = 'complex',
                  value = 'c12*c13',
                  texname = '\\text{CKM11}')

CKM12 = Parameter(name = 'CKM12',
                  nature = 'internal',
                  type = 'complex',
                  value = 'c13*s12',
                  texname = '\\text{CKM12}')

CKM21 = Parameter(name = 'CKM21',
                  nature = 'internal',
                  type = 'complex',
                  value = '-(c23*s12) - c12*s13*s23',
                  texname = '\\text{CKM21}')

CKM22 = Parameter(name = 'CKM22',
                  nature = 'internal',
                  type = 'complex',
                  value = 'c12*c23 - s12*s13*s23',
                  texname = '\\text{CKM22}')

CKM23 = Parameter(name = 'CKM23',
                  nature = 'internal',
                  type = 'complex',
                  value = 'c13*s23',
                  texname = '\\text{CKM23}')

CKM31 = Parameter(name = 'CKM31',
                  nature = 'internal',
                  type = 'complex',
                  value = '-(c12*c23*s13) + s12*s23',
                  texname = '\\text{CKM31}')

CKM32 = Parameter(name = 'CKM32',
                  nature = 'internal',
                  type = 'complex',
                  value = '-(c23*s12*s13) - c12*s23',
                  texname = '\\text{CKM32}')

CKM33 = Parameter(name = 'CKM33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'c13*c23',
                  texname = '\\text{CKM33}')

cw = Parameter(name = 'cw',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(1 - sw2)',
               texname = 'c_w')

MG1 = Parameter(name = 'MG1',
                nature = 'internal',
                type = 'real',
                value = 'cmath.sqrt(R**(-2) - (0.0011569056214749236*gmass**2)/R**2 + (23*gmass**2*cmath.log(LR**2))/(32.*cmath.pi**2*R**2))',
                texname = '\\text{MG1}')

ee = Parameter(name = 'ee',
               nature = 'internal',
               type = 'real',
               value = '2.0*cmath.sqrt(aEW)*cmath.sqrt(cmath.pi)',
               texname = 'e')

MW = Parameter(name = 'MW',
               nature = 'internal',
               type = 'real',
               value = 'cw*ZM',
               texname = '\\text{MW}')

g1 = Parameter(name = 'g1',
               nature = 'internal',
               type = 'real',
               value = 'ee/cw',
               texname = 'g_1')

gw = Parameter(name = 'gw',
               nature = 'internal',
               type = 'real',
               value = 'ee/sw',
               texname = 'g_w')

md2 = Parameter(name = 'md2',
                nature = 'internal',
                type = 'real',
                value = '(g1**2/4. + 3*gmass**2)*mn',
                texname = '\\text{md2}')

me2 = Parameter(name = 'me2',
                nature = 'internal',
                type = 'real',
                value = '(9*g1**2*mn)/4.',
                texname = '\\text{me2}')

mL2 = Parameter(name = 'mL2',
                nature = 'internal',
                type = 'real',
                value = '((9.0*g1**2)/16. + (27.0*gw**2)/16.)*mn',
                texname = '\\text{mL2}')

mQ2 = Parameter(name = 'mQ2',
                nature = 'internal',
                type = 'real',
                value = '(g1**2/16. + 3.0*gmass**2 + (27.0*gw**2)/16.)*mn',
                texname = '\\text{mQ2}')

mu2 = Parameter(name = 'mu2',
                nature = 'internal',
                type = 'real',
                value = '(g1**2 + 3*gmass**2)*mn',
                texname = '\\text{mu2}')

v = Parameter(name = 'v',
              nature = 'internal',
              type = 'real',
              value = '(2.0*MW*sw)/ee',
              texname = 'v')

cH2 = Parameter(name = 'cH2',
                nature = 'internal',
                type = 'real',
                value = '(mn*((3.0*g1**2)/4. + (3.0*gw**2)/2. - MH**2/v**2))/R',
                texname = '\\text{cH2}')

DMc = Parameter(name = 'DMc',
                nature = 'internal',
                type = 'real',
                value = 'mQ2 + 1.0/R',
                texname = '\\text{DMc}')

DMd = Parameter(name = 'DMd',
                nature = 'internal',
                type = 'real',
                value = 'mQ2 + 1.0/R',
                texname = '\\text{DMd}')

DMe = Parameter(name = 'DMe',
                nature = 'internal',
                type = 'real',
                value = 'mL2 + 1.0/R',
                texname = '\\text{DMe}')

DMmu = Parameter(name = 'DMmu',
                 nature = 'internal',
                 type = 'real',
                 value = 'mL2 + 1.0/R',
                 texname = '\\text{DMmu}')

DMn1 = Parameter(name = 'DMn1',
                 nature = 'internal',
                 type = 'real',
                 value = 'mL2 + 1.0/R',
                 texname = '\\text{DMn1}')

DMn2 = Parameter(name = 'DMn2',
                 nature = 'internal',
                 type = 'real',
                 value = 'mL2 + 1.0/R',
                 texname = '\\text{DMn2}')

DMn3 = Parameter(name = 'DMn3',
                 nature = 'internal',
                 type = 'real',
                 value = 'mL2 + 1.0/R',
                 texname = '\\text{DMn3}')

DMs = Parameter(name = 'DMs',
                nature = 'internal',
                type = 'real',
                value = 'mQ2 + 1.0/R',
                texname = '\\text{DMs}')

DMtau = Parameter(name = 'DMtau',
                  nature = 'internal',
                  type = 'real',
                  value = 'mL2 + 1.0/R',
                  texname = '\\text{DMtau}')

DMu = Parameter(name = 'DMu',
                nature = 'internal',
                type = 'real',
                value = 'mQ2 + 1.0/R',
                texname = '\\text{DMu}')

MM11 = Parameter(name = 'MM11',
                 nature = 'internal',
                 type = 'real',
                 value = 'R**(-2) - (0.015039773079174007*g1**2)/R**2 + (g1**2*v**2)/4. - (g1**2*cmath.log(LR**2))/(96.*cmath.pi**2*R**2)',
                 texname = '\\text{MM11}')

MM12 = Parameter(name = 'MM12',
                 nature = 'internal',
                 type = 'real',
                 value = '(g1*gw*v**2)/4.',
                 texname = '\\text{MM12}')

MM22 = Parameter(name = 'MM22',
                 nature = 'internal',
                 type = 'real',
                 value = 'R**(-2) - (0.0019281760357915392*gw**2)/R**2 + (gw**2*v**2)/4. + (15*gw**2*cmath.log(LR**2))/(32.*cmath.pi**2*R**2)',
                 texname = '\\text{MM22}')

mQt = Parameter(name = 'mQt',
                nature = 'internal',
                type = 'real',
                value = '(-3.0*mn*ymt**2)/(2.*v**2)',
                texname = '\\text{mQt}')

mt0 = Parameter(name = 'mt0',
                nature = 'internal',
                type = 'real',
                value = '(-3.0*mn*ymt**2)/v**2',
                texname = '\\text{mt0}')

SMc = Parameter(name = 'SMc',
                nature = 'internal',
                type = 'real',
                value = 'mu2 + 1.0/R',
                texname = '\\text{SMc}')

SMd = Parameter(name = 'SMd',
                nature = 'internal',
                type = 'real',
                value = 'md2 + 1.0/R',
                texname = '\\text{SMd}')

SMe = Parameter(name = 'SMe',
                nature = 'internal',
                type = 'real',
                value = 'me2 + 1.0/R',
                texname = '\\text{SMe}')

SMmu = Parameter(name = 'SMmu',
                 nature = 'internal',
                 type = 'real',
                 value = 'me2 + 1.0/R',
                 texname = '\\text{SMmu}')

SMs = Parameter(name = 'SMs',
                nature = 'internal',
                type = 'real',
                value = 'md2 + 1.0/R',
                texname = '\\text{SMs}')

SMtau = Parameter(name = 'SMtau',
                  nature = 'internal',
                  type = 'real',
                  value = 'me2 + 1.0/R',
                  texname = '\\text{SMtau}')

SMu = Parameter(name = 'SMu',
                nature = 'internal',
                type = 'real',
                value = 'mu2 + 1.0/R',
                texname = '\\text{SMu}')

YukD33 = Parameter(name = 'YukD33',
                   nature = 'internal',
                   type = 'complex',
                   value = '(ymb*cmath.sqrt(2.0))/v',
                   texname = '\\text{YukD33}')

lam = Parameter(name = 'lam',
                nature = 'internal',
                type = 'real',
                value = 'MH**2/(2.*v**2)',
                texname = '\\text{lam}')

yb = Parameter(name = 'yb',
               nature = 'internal',
               type = 'real',
               value = '(ymb*cmath.sqrt(2.0))/v',
               texname = '\\text{yb}')

yc = Parameter(name = 'yc',
               nature = 'internal',
               type = 'real',
               value = '(ymc*cmath.sqrt(2.0))/v',
               texname = '\\text{yc}')

yt = Parameter(name = 'yt',
               nature = 'internal',
               type = 'real',
               value = '(ymt*cmath.sqrt(2.0))/v',
               texname = '\\text{yt}')

ytau = Parameter(name = 'ytau',
                 nature = 'internal',
                 type = 'real',
                 value = '(ymtau*cmath.sqrt(2.0))/v',
                 texname = '\\text{ytau}')

Delta = Parameter(name = 'Delta',
                  nature = 'internal',
                  type = 'real',
                  value = 'cmath.sqrt(4.0*MM12**2 + (MM11 - MM22)**2)',
                  texname = '\\text{Delta}')

DMb = Parameter(name = 'DMb',
                nature = 'internal',
                type = 'real',
                value = '(-md2 + mQ2 + mQt + cmath.sqrt((md2 + mQ2 + mQt + 2.0/R)**2 + 4.0*ymb**2))/2.',
                texname = '\\text{DMb}')

DMtop = Parameter(name = 'DMtop',
                  nature = 'internal',
                  type = 'real',
                  value = '(-mQ2 - mQt + mt0 + mu2 + cmath.sqrt((mQ2 + mQt + mt0 + mu2 + 2.0/R)**2 + 4.0*ymt**2))/2.',
                  texname = '\\text{DMtop}')

MH1 = Parameter(name = 'MH1',
                nature = 'internal',
                type = 'real',
                value = 'cmath.sqrt(cH2 - MH**2 + R**(-2))',
                texname = '\\text{MH1}')

muH = Parameter(name = 'muH',
                nature = 'internal',
                type = 'real',
                value = 'cmath.sqrt(lam*v**2)',
                texname = '\\mu ')

MW1 = Parameter(name = 'MW1',
                nature = 'internal',
                type = 'real',
                value = 'cmath.sqrt(MM22)',
                texname = '\\text{MW1}')

SMb = Parameter(name = 'SMb',
                nature = 'internal',
                type = 'real',
                value = '(md2 - mQ2 - mQt + cmath.sqrt((md2 + mQ2 + mQt + 2.0/R)**2 + 4.0*ymb**2))/2.',
                texname = '\\text{SMb}')

SMtop = Parameter(name = 'SMtop',
                  nature = 'internal',
                  type = 'real',
                  value = '(mQ2 + mQt - mt0 - mu2 + cmath.sqrt((mQ2 + mQt + mt0 + mu2 + 2.0/R)**2 + 4.0*ymt**2))/2.',
                  texname = '\\text{SMtop}')

MB1 = Parameter(name = 'MB1',
                nature = 'internal',
                type = 'real',
                value = 'cmath.sqrt(-Delta + MM11 + MM22)/cmath.sqrt(2.0)',
                texname = '\\text{MB1}')

MZ1 = Parameter(name = 'MZ1',
                nature = 'internal',
                type = 'real',
                value = 'cmath.sqrt(Delta + MM11 + MM22)/cmath.sqrt(2.0)',
                texname = '\\text{MZ1}')

