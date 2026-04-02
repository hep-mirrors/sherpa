# This file was automatically created by FeynRules $Revision: 915 $
# Mathematica version: 7.0 for Mac OS X x86 (64-bit) (November 11, 2008)
# Date: Mon 23 Jan 2012 11:06:17


from .object_library import all_vertices, Vertex
from . import particles as P
from . import couplings as C
from . import lorentz as L


# V_1 = Vertex(name = 'V_1',
#              particles = [ P.H0, P.H0, P.H0, P.H0 ],
#              color = [ '1' ],
#              lorentz = [ L.SSSS1 ],
#              couplings = {(0,0):C.GC_26})

# V_2 = Vertex(name = 'V_2',
#              particles = [ P.H0, P.H0, P.H0 ],
#              color = [ '1' ],
#              lorentz = [ L.SSS1 ],
#              couplings = {(0,0):C.GC_33})

# V_3 = Vertex(name = 'V_3',
#              particles = [ P.G0, P.G0, P.G0 ],
#              color = [ 'f(1,2,3)' ],
#              lorentz = [ L.VVV1 ],
#              couplings = {(0,0):C.GC_1})

# V_4 = Vertex(name = 'V_4',
#              particles = [ P.G0, P.G0, P.G0, P.G0 ],
#              color = [ 'f(-1,1,2)*f(3,4,-1)', 'f(-1,1,3)*f(2,4,-1)', 'f(-1,1,4)*f(2,3,-1)' ],
#              lorentz = [ L.VVVV1, L.VVVV3, L.VVVV4 ],
#              couplings = {(1,1):C.GC_4,(0,0):C.GC_4,(2,2):C.GC_4})


# V_8 = Vertex(name = 'V_8',
#              particles = [ P.A0, P.W0__minus__, P.W0__plus__ ],
#              color = [ '1' ],
#              lorentz = [ L.VVV1 ],
#              couplings = {(0,0):C.GC_28})


# V_15 = Vertex(name = 'V_15',
#               particles = [ P.W0__minus__, P.W0__plus__, P.H0, P.H0 ],
#               color = [ '1' ],
#               lorentz = [ L.VVSS1 ],
#               couplings = {(0,0):C.GC_27})

# V_16 = Vertex(name = 'V_16',
#               particles = [ P.W0__minus__, P.W0__plus__, P.H0 ],
#               color = [ '1' ],
#               lorentz = [ L.VVS1 ],
#               couplings = {(0,0):C.GC_34})

# V_17 = Vertex(name = 'V_17',
#               particles = [ P.A0, P.A0, P.W0__minus__, P.W0__plus__ ],
#               color = [ '1' ],
#               lorentz = [ L.VVVV2 ],
#               couplings = {(0,0):C.GC_31})

# V_18 = Vertex(name = 'V_18',
#               particles = [ P.W0__minus__, P.W0__plus__, P.Z0 ],
#               color = [ '1' ],
#               lorentz = [ L.VVV1 ],
#               couplings = {(0,0):C.GC_19})


# V_21 = Vertex(name = 'V_21',
#               particles = [ P.W0__minus__, P.W0__minus__, P.W0__plus__, P.W0__plus__ ],
#               color = [ '1' ],
#               lorentz = [ L.VVVV2 ],
#               couplings = {(0,0):C.GC_20})


# V_23 = Vertex(name = 'V_23',
#               particles = [ P.b0__tilde__, P.b0, P.H0 ],
#               color = [ 'Identity(1,2)' ],
#               lorentz = [ L.FFS1 ],
#               couplings = {(0,0):C.GC_55})

# V_24 = Vertex(name = 'V_24',
#               particles = [ P.tt0__plus__, P.tt0__minus__, P.H0 ],
#               color = [ '1' ],
#               lorentz = [ L.FFS1 ],
#               couplings = {(0,0):C.GC_58})

# V_25 = Vertex(name = 'V_25',
#               particles = [ P.c0__tilde__, P.c0, P.H0 ],
#               color = [ 'Identity(1,2)' ],
#               lorentz = [ L.FFS1 ],
#               couplings = {(0,0):C.GC_56})

# V_26 = Vertex(name = 'V_26',
#               particles = [ P.t0__tilde__, P.t0, P.H0 ],
#               color = [ 'Identity(1,2)' ],
#               lorentz = [ L.FFS1 ],
#               couplings = {(0,0):C.GC_57})


# V_28 = Vertex(name = 'V_28',
#               particles = [ P.A0, P.W0__minus__, P.W0__plus__, P.Z0 ],
#               color = [ '1' ],
#               lorentz = [ L.VVVV5 ],
#               couplings = {(0,0):C.GC_30})

# V_29 = Vertex(name = 'V_29',
#               particles = [ P.Z0, P.Z0, P.H0, P.H0 ],
#               color = [ '1' ],
#               lorentz = [ L.VVSS1 ],
#               couplings = {(0,0):C.GC_32})

# V_30 = Vertex(name = 'V_30',
#               particles = [ P.Z0, P.Z0, P.H0 ],
#               color = [ '1' ],
#               lorentz = [ L.VVS1 ],
#               couplings = {(0,0):C.GC_35})


# V_32 = Vertex(name = 'V_32',
#               particles = [ P.W0__minus__, P.W0__plus__, P.Z0, P.Z0 ],
#               color = [ '1' ],
#               lorentz = [ L.VVVV2 ],
#               couplings = {(0,0):C.GC_25})


# V_48 = Vertex(name = 'V_48',
#               particles = [ P.b0__tilde__, P.c0, P.W0__minus__ ],
#               color = [ 'Identity(1,2)' ],
#               lorentz = [ L.FFV2 ],
#               couplings = {(0,0):C.GC_15})

# V_49 = Vertex(name = 'V_49',
#               particles = [ P.d0__tilde__, P.c0, P.W0__minus__ ],
#               color = [ 'Identity(1,2)' ],
#               lorentz = [ L.FFV2 ],
#               couplings = {(0,0):C.GC_13})

# V_50 = Vertex(name = 'V_50',
#               particles = [ P.s0__tilde__, P.c0, P.W0__minus__ ],
#               color = [ 'Identity(1,2)' ],
#               lorentz = [ L.FFV2 ],
#               couplings = {(0,0):C.GC_14})

# V_51 = Vertex(name = 'V_51',
#               particles = [ P.b0__tilde__, P.t0, P.W0__minus__ ],
#               color = [ 'Identity(1,2)' ],
#               lorentz = [ L.FFV2 ],
#               couplings = {(0,0):C.GC_18})

# V_52 = Vertex(name = 'V_52',
#               particles = [ P.d0__tilde__, P.t0, P.W0__minus__ ],
#               color = [ 'Identity(1,2)' ],
#               lorentz = [ L.FFV2 ],
#               couplings = {(0,0):C.GC_16})

# V_53 = Vertex(name = 'V_53',
#               particles = [ P.s0__tilde__, P.t0, P.W0__minus__ ],
#               color = [ 'Identity(1,2)' ],
#               lorentz = [ L.FFV2 ],
#               couplings = {(0,0):C.GC_17})

# V_54 = Vertex(name = 'V_54',
#               particles = [ P.b0__tilde__, P.u0, P.W0__minus__ ],
#               color = [ 'Identity(1,2)' ],
#               lorentz = [ L.FFV2 ],
#               couplings = {(0,0):C.GC_12})

# V_55 = Vertex(name = 'V_55',
#               particles = [ P.d0__tilde__, P.u0, P.W0__minus__ ],
#               color = [ 'Identity(1,2)' ],
#               lorentz = [ L.FFV2 ],
#               couplings = {(0,0):C.GC_10})

# V_56 = Vertex(name = 'V_56',
#               particles = [ P.s0__tilde__, P.u0, P.W0__minus__ ],
#               color = [ 'Identity(1,2)' ],
#               lorentz = [ L.FFV2 ],
#               couplings = {(0,0):C.GC_11})
# V_93 = Vertex(name = 'V_93',
#               particles = [ P.c0__tilde__, P.b0, P.W0__plus__ ],
#               color = [ 'Identity(1,2)' ],
#               lorentz = [ L.FFV2 ],
#               couplings = {(0,0):C.GC_64})

# V_94 = Vertex(name = 'V_94',
#               particles = [ P.t0__tilde__, P.b0, P.W0__plus__ ],
#               color = [ 'Identity(1,2)' ],
#               lorentz = [ L.FFV2 ],
#               couplings = {(0,0):C.GC_67})

# V_95 = Vertex(name = 'V_95',
#               particles = [ P.u0__tilde__, P.b0, P.W0__plus__ ],
#               color = [ 'Identity(1,2)' ],
#               lorentz = [ L.FFV2 ],
#               couplings = {(0,0):C.GC_61})

# V_96 = Vertex(name = 'V_96',
#               particles = [ P.c0__tilde__, P.d0, P.W0__plus__ ],
#               color = [ 'Identity(1,2)' ],
#               lorentz = [ L.FFV2 ],
#               couplings = {(0,0):C.GC_62})

# V_97 = Vertex(name = 'V_97',
#               particles = [ P.t0__tilde__, P.d0, P.W0__plus__ ],
#               color = [ 'Identity(1,2)' ],
#               lorentz = [ L.FFV2 ],
#               couplings = {(0,0):C.GC_65})

# V_98 = Vertex(name = 'V_98',
#               particles = [ P.u0__tilde__, P.d0, P.W0__plus__ ],
#               color = [ 'Identity(1,2)' ],
#               lorentz = [ L.FFV2 ],
#               couplings = {(0,0):C.GC_59})

# V_99 = Vertex(name = 'V_99',
#               particles = [ P.c0__tilde__, P.s0, P.W0__plus__ ],
#               color = [ 'Identity(1,2)' ],
#               lorentz = [ L.FFV2 ],
#               couplings = {(0,0):C.GC_63})

# V_100 = Vertex(name = 'V_100',
#                particles = [ P.t0__tilde__, P.s0, P.W0__plus__ ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_66})

# V_101 = Vertex(name = 'V_101',
#                particles = [ P.u0__tilde__, P.s0, P.W0__plus__ ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_60})

# V_111 = Vertex(name = 'V_111',
#                particles = [ P.b0__tilde__, P.b0, P.Z0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2, L.FFV3 ],
#                couplings = {(0,0):C.GC_47,(0,1):C.GC_54})

# V_112 = Vertex(name = 'V_112',
#                particles = [ P.d0__tilde__, P.d0, P.Z0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2, L.FFV3 ],
#                couplings = {(0,0):C.GC_47,(0,1):C.GC_54})

# V_113 = Vertex(name = 'V_113',
#                particles = [ P.s0__tilde__, P.s0, P.Z0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2, L.FFV3 ],
#                couplings = {(0,0):C.GC_47,(0,1):C.GC_54})




# V_116 = Vertex(name = 'V_116',
#                particles = [ P.Ds1__tilde__, P.Ds1, P.Z0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_47})

# V_117 = Vertex(name = 'V_117',
#                particles = [ P.Dc1__tilde__, P.Dc1, P.Z0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_48})

# V_118 = Vertex(name = 'V_118',
#                particles = [ P.Dt1__tilde__, P.Dt1, P.Z0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_48})

# V_119 = Vertex(name = 'V_119',
#                particles = [ P.Du1__tilde__, P.Du1, P.Z0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_48})

# V_120 = Vertex(name = 'V_120',
#                particles = [ P.c0__tilde__, P.c0, P.Z0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2, L.FFV3 ],
#                couplings = {(0,0):C.GC_48,(0,1):C.GC_51})

# V_121 = Vertex(name = 'V_121',
#                particles = [ P.t0__tilde__, P.t0, P.Z0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2, L.FFV3 ],
#                couplings = {(0,0):C.GC_48,(0,1):C.GC_51})

# V_122 = Vertex(name = 'V_122',
#                particles = [ P.u0__tilde__, P.u0, P.Z0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2, L.FFV3 ],
#                couplings = {(0,0):C.GC_48,(0,1):C.GC_51})


# # V_135 = Vertex(name = 'V_135',
# #                particles = [ P.e0__plus__, P.e0__minus__, P.A0 ],
# #                color = [ '1' ],
# #                lorentz = [ L.FFV2, L.FFV3 ],
# #                couplings = {(0,0):C.GC_37,(0,1):C.GC_42})

# V_136 = Vertex(name = 'V_136',
#                particles = [ P.m0__plus__, P.m0__minus__, P.A0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2, L.FFV3 ],
#                couplings = {(0,0):C.GC_37,(0,1):C.GC_42})

# V_137 = Vertex(name = 'V_137',
#                particles = [ P.tt0__plus__, P.tt0__minus__, P.A0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2, L.FFV3 ],
#                couplings = {(0,0):C.GC_37,(0,1):C.GC_42})

# # V_138 = Vertex(name = 'V_138',
#                particles = [ P.e1L__plus__, P.e0__minus__, P.A1 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_36})

# V_139 = Vertex(name = 'V_139',
#                particles = [ P.m1L__plus__, P.m0__minus__, P.A1 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_36})

# V_140 = Vertex(name = 'V_140',
#                particles = [ P.t1L__plus__, P.tt0__minus__, P.A1 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_36})

# V_141 = Vertex(name = 'V_141',
#                particles = [ P.e0__plus__, P.e1L__minus__, P.A1 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_36})

# V_142 = Vertex(name = 'V_142',
#                particles = [ P.m0__plus__, P.m1L__minus__, P.A1 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_36})

# V_143 = Vertex(name = 'V_143',
#                particles = [ P.tt0__plus__, P.t1L__minus__, P.A1 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_36})

# V_144 = Vertex(name = 'V_144',
#                particles = [ P.e1L__plus__, P.e1L__minus__, P.A0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_37})

# V_145 = Vertex(name = 'V_145',
#                particles = [ P.m1L__plus__, P.m1L__minus__, P.A0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_37})

# V_146 = Vertex(name = 'V_146',
#                particles = [ P.t1L__plus__, P.t1L__minus__, P.A0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_37})

# V_147 = Vertex(name = 'V_147',
#                particles = [ P.ve1__tilde__, P.ve1, P.A0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_38})

# V_148 = Vertex(name = 'V_148',
#                particles = [ P.vm1__tilde__, P.vm1, P.A0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_38})

# V_149 = Vertex(name = 'V_149',
#                particles = [ P.vt1__tilde__, P.vt1, P.A0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_38})

# V_150 = Vertex(name = 'V_150',
#                particles = [ P.ve0__tilde__, P.ve1, P.A1 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_36})

# V_151 = Vertex(name = 'V_151',
#                particles = [ P.vm0__tilde__, P.vm1, P.A1 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_36})

# V_152 = Vertex(name = 'V_152',
#                particles = [ P.vt0__tilde__, P.vt1, P.A1 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_36})

# V_153 = Vertex(name = 'V_153',
#                particles = [ P.ve1__tilde__, P.ve0, P.A1 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_36})

# V_154 = Vertex(name = 'V_154',
#                particles = [ P.vm1__tilde__, P.vm0, P.A1 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_36})

# V_155 = Vertex(name = 'V_155',
#                particles = [ P.vt1__tilde__, P.vt0, P.A1 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_36})

# V_156 = Vertex(name = 'V_156',
#                particles = [ P.ve0__tilde__, P.ve0, P.A0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_38})

# V_157 = Vertex(name = 'V_157',
#                particles = [ P.vm0__tilde__, P.vm0, P.A0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_38})

# V_158 = Vertex(name = 'V_158',
#                particles = [ P.vt0__tilde__, P.vt0, P.A0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_38})

# V_159 = Vertex(name = 'V_159',
#                particles = [ P.e1L__plus__, P.ve1, P.W0__minus__ ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_9})

# V_160 = Vertex(name = 'V_160',
#                particles = [ P.m1L__plus__, P.vm1, P.W0__minus__ ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_9})

# V_161 = Vertex(name = 'V_161',
#                particles = [ P.t1L__plus__, P.vt1, P.W0__minus__ ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_9})

# V_162 = Vertex(name = 'V_162',
#                particles = [ P.e0__plus__, P.ve0, P.W0__minus__ ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_9})

# V_163 = Vertex(name = 'V_163',
#                particles = [ P.m0__plus__, P.vm0, P.W0__minus__ ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_9})

# V_164 = Vertex(name = 'V_164',
#                particles = [ P.tt0__plus__, P.vt0, P.W0__minus__ ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_9})


# V_177 = Vertex(name = 'V_177',
#                particles = [ P.ve0__tilde__, P.e0__minus__, P.W0__plus__ ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_9})

# V_178 = Vertex(name = 'V_178',
#                particles = [ P.vm0__tilde__, P.m0__minus__, P.W0__plus__ ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_9})

# V_179 = Vertex(name = 'V_179',
#                particles = [ P.vt0__tilde__, P.tt0__minus__, P.W0__plus__ ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_9})

# V_180 = Vertex(name = 'V_180',
#                particles = [ P.ve1__tilde__, P.e1L__minus__, P.W0__plus__ ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_9})

# V_181 = Vertex(name = 'V_181',
#                particles = [ P.vm1__tilde__, P.m1L__minus__, P.W0__plus__ ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_9})

# V_182 = Vertex(name = 'V_182',
#                particles = [ P.vt1__tilde__, P.t1L__minus__, P.W0__plus__ ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_9})

# V_183 = Vertex(name = 'V_183',
#                particles = [ P.e0__plus__, P.e0__minus__, P.Z0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2, L.FFV3 ],
#                couplings = {(0,0):C.GC_39,(0,1):C.GC_43})

# V_184 = Vertex(name = 'V_184',
#                particles = [ P.m0__plus__, P.m0__minus__, P.Z0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2, L.FFV3 ],
#                couplings = {(0,0):C.GC_39,(0,1):C.GC_43})

# V_185 = Vertex(name = 'V_185',
#                particles = [ P.tt0__plus__, P.tt0__minus__, P.Z0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2, L.FFV3 ],
#                couplings = {(0,0):C.GC_39,(0,1):C.GC_43})

# V_186 = Vertex(name = 'V_186',
#                particles = [ P.e1L__plus__, P.e1L__minus__, P.Z0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_39})

# V_187 = Vertex(name = 'V_187',
#                particles = [ P.m1L__plus__, P.m1L__minus__, P.Z0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_39})

# V_188 = Vertex(name = 'V_188',
#                particles = [ P.t1L__plus__, P.t1L__minus__, P.Z0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_39})

# V_189 = Vertex(name = 'V_189',
#                particles = [ P.ve1__tilde__, P.ve1, P.Z0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_40})

# V_190 = Vertex(name = 'V_190',
#                particles = [ P.vm1__tilde__, P.vm1, P.Z0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_40})

# V_191 = Vertex(name = 'V_191',
#                particles = [ P.vt1__tilde__, P.vt1, P.Z0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_40})

# V_192 = Vertex(name = 'V_192',
#                particles = [ P.ve0__tilde__, P.ve0, P.Z0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_40})

# V_193 = Vertex(name = 'V_193',
#                particles = [ P.vm0__tilde__, P.vm0, P.Z0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_40})

# V_194 = Vertex(name = 'V_194',
#                particles = [ P.vt0__tilde__, P.vt0, P.Z0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_40})


# V_207 = Vertex(name = 'V_207',
#                particles = [ P.Sc1__tilde__, P.Sc1, P.A0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_50})

# V_208 = Vertex(name = 'V_208',
#                particles = [ P.St1__tilde__, P.St1, P.A0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_50})

# V_209 = Vertex(name = 'V_209',
#                particles = [ P.Su1__tilde__, P.Su1, P.A0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_50})

# V_210 = Vertex(name = 'V_210',
#                particles = [ P.Sc1__tilde__, P.Sc1, P.G0 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_2})

# V_211 = Vertex(name = 'V_211',
#                particles = [ P.St1__tilde__, P.St1, P.G0 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_2})

# V_212 = Vertex(name = 'V_212',
#                particles = [ P.Su1__tilde__, P.Su1, P.G0 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_2})

# V_213 = Vertex(name = 'V_213',
#                particles = [ P.Sc1__tilde__, P.Sc1, P.Z0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_51})

# V_214 = Vertex(name = 'V_214',
#                particles = [ P.St1__tilde__, P.St1, P.Z0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_51})

# V_215 = Vertex(name = 'V_215',
#                particles = [ P.Su1__tilde__, P.Su1, P.Z0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_51})

# V_216 = Vertex(name = 'V_216',
#                particles = [ P.Sb1__tilde__, P.Sb1, P.A0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_53})

# V_217 = Vertex(name = 'V_217',
#                particles = [ P.Sd1__tilde__, P.Sd1, P.A0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_53})

# V_218 = Vertex(name = 'V_218',
#                particles = [ P.Ss1__tilde__, P.Ss1, P.A0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_53})

# V_219 = Vertex(name = 'V_219',
#                particles = [ P.Sb1__tilde__, P.Sb1, P.G0 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_2})

# V_220 = Vertex(name = 'V_220',
#                particles = [ P.Sd1__tilde__, P.Sd1, P.G0 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_2})

# V_221 = Vertex(name = 'V_221',
#                particles = [ P.Ss1__tilde__, P.Ss1, P.G0 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_2})

# V_222 = Vertex(name = 'V_222',
#                particles = [ P.Sb1__tilde__, P.Sb1, P.Z0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_54})

# V_223 = Vertex(name = 'V_223',
#                particles = [ P.Sd1__tilde__, P.Sd1, P.Z0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_54})

# V_224 = Vertex(name = 'V_224',
#                particles = [ P.Ss1__tilde__, P.Ss1, P.Z0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_54})

# V_225 = Vertex(name = 'V_225',
#                particles = [ P.e1R__plus__, P.e1R__minus__, P.A0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_42})

# V_226 = Vertex(name = 'V_226',
#                particles = [ P.m1R__plus__, P.m1R__minus__, P.A0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_42})

# V_227 = Vertex(name = 'V_227',
#                particles = [ P.t1R__plus__, P.t1R__minus__, P.A0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_42})

# V_228 = Vertex(name = 'V_228',
#                particles = [ P.e1R__plus__, P.e1R__minus__, P.Z0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_43})

# V_229 = Vertex(name = 'V_229',
#                particles = [ P.m1R__plus__, P.m1R__minus__, P.Z0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_43})

# V_230 = Vertex(name = 'V_230',
#                particles = [ P.t1R__plus__, P.t1R__minus__, P.Z0 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_43})

# V_231 = Vertex(name = 'V_231',
#                particles = [ P.b0__tilde__, P.b0, P.A0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2, L.FFV3 ],
#                couplings = {(0,0):C.GC_45,(0,1):C.GC_53})

# V_232 = Vertex(name = 'V_232',
#                particles = [ P.d0__tilde__, P.d0, P.A0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2, L.FFV3 ],
#                couplings = {(0,0):C.GC_45,(0,1):C.GC_53})

# V_233 = Vertex(name = 'V_233',
#                particles = [ P.s0__tilde__, P.s0, P.A0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2, L.FFV3 ],
#                couplings = {(0,0):C.GC_45,(0,1):C.GC_53})



# V_236 = Vertex(name = 'V_236',
#                particles = [ P.Ds1__tilde__, P.s0, P.A1 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_44})

# V_237 = Vertex(name = 'V_237',
#                particles = [ P.b0__tilde__, P.Db1, P.A1 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_44})

# V_238 = Vertex(name = 'V_238',
#                particles = [ P.d0__tilde__, P.Dd1, P.A1 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_44})

# V_239 = Vertex(name = 'V_239',
#                particles = [ P.s0__tilde__, P.Ds1, P.A1 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_44})





# V_242 = Vertex(name = 'V_242',
#                particles = [ P.Ds1__tilde__, P.Ds1, P.A0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_45})

# V_243 = Vertex(name = 'V_243',
#                particles = [ P.Dc1__tilde__, P.Dc1, P.A0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_46})

# V_244 = Vertex(name = 'V_244',
#                particles = [ P.Dt1__tilde__, P.Dt1, P.A0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_46})

# V_245 = Vertex(name = 'V_245',
#                particles = [ P.Du1__tilde__, P.Du1, P.A0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_46})

# V_246 = Vertex(name = 'V_246',
#                particles = [ P.c0__tilde__, P.Dc1, P.A1 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_44})

# V_247 = Vertex(name = 'V_247',
#                particles = [ P.t0__tilde__, P.Dt1, P.A1 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_44})

# V_248 = Vertex(name = 'V_248',
#                particles = [ P.u0__tilde__, P.Du1, P.A1 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_44})

# V_249 = Vertex(name = 'V_249',
#                particles = [ P.Dc1__tilde__, P.c0, P.A1 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_44})

# V_250 = Vertex(name = 'V_250',
#                particles = [ P.Dt1__tilde__, P.t0, P.A1 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_44})

# V_251 = Vertex(name = 'V_251',
#                particles = [ P.Du1__tilde__, P.u0, P.A1 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2 ],
#                couplings = {(0,0):C.GC_44})

# V_252 = Vertex(name = 'V_252',
#                particles = [ P.c0__tilde__, P.c0, P.A0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2, L.FFV3 ],
#                couplings = {(0,0):C.GC_46,(0,1):C.GC_50})

# V_253 = Vertex(name = 'V_253',
#                particles = [ P.t0__tilde__, P.t0, P.A0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2, L.FFV3 ],
#                couplings = {(0,0):C.GC_46,(0,1):C.GC_50})

# V_254 = Vertex(name = 'V_254',
#                particles = [ P.u0__tilde__, P.u0, P.A0 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV2, L.FFV3 ],
#                couplings = {(0,0):C.GC_46,(0,1):C.GC_50})

# V_255 = Vertex(name = 'V_255',
#                particles = [ P.b0__tilde__, P.b0, P.G0 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_2})

# V_256 = Vertex(name = 'V_256',
#                particles = [ P.d0__tilde__, P.d0, P.G0 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_2})

# V_257 = Vertex(name = 'V_257',
#                particles = [ P.s0__tilde__, P.s0, P.G0 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_2})



# V_260 = Vertex(name = 'V_260',
#                particles = [ P.Ds1__tilde__, P.Ds1, P.G0 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_2})

# V_261 = Vertex(name = 'V_261',
#                particles = [ P.Dc1__tilde__, P.Dc1, P.G0 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_2})

# V_262 = Vertex(name = 'V_262',
#                particles = [ P.Dt1__tilde__, P.Dt1, P.G0 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_2})

# V_263 = Vertex(name = 'V_263',
#                particles = [ P.Du1__tilde__, P.Du1, P.G0 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_2})

# V_264 = Vertex(name = 'V_264',
#                particles = [ P.c0__tilde__, P.c0, P.G0 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_2})

# V_265 = Vertex(name = 'V_265',
#                particles = [ P.t0__tilde__, P.t0, P.G0 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_2})

# V_266 = Vertex(name = 'V_266',
#                particles = [ P.u0__tilde__, P.u0, P.G0 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV1 ],
#                couplings = {(0,0):C.GC_2})

# V_279 = Vertex(name = 'V_279',
#                particles = [ P.c0__tilde__, P.Sc1, P.A1 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_49})

# V_280 = Vertex(name = 'V_280',
#                particles = [ P.t0__tilde__, P.St1, P.A1 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_49})

# V_281 = Vertex(name = 'V_281',
#                particles = [ P.u0__tilde__, P.Su1, P.A1 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_49})

# V_282 = Vertex(name = 'V_282',
#                particles = [ P.Sc1__tilde__, P.c0, P.A1 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_49})

# V_283 = Vertex(name = 'V_283',
#                particles = [ P.St1__tilde__, P.t0, P.A1 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_49})

# V_284 = Vertex(name = 'V_284',
#                particles = [ P.Su1__tilde__, P.u0, P.A1 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_49})

# V_285 = Vertex(name = 'V_285',
#                particles = [ P.c0__tilde__, P.Sc1, P.G1 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_2})

# V_286 = Vertex(name = 'V_286',
#                particles = [ P.t0__tilde__, P.St1, P.G1 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_2})

# V_287 = Vertex(name = 'V_287',
#                particles = [ P.u0__tilde__, P.Su1, P.G1 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_2})

# V_288 = Vertex(name = 'V_288',
#                particles = [ P.Sc1__tilde__, P.c0, P.G1 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_2})

# V_289 = Vertex(name = 'V_289',
#                particles = [ P.St1__tilde__, P.t0, P.G1 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_2})

# V_290 = Vertex(name = 'V_290',
#                particles = [ P.Su1__tilde__, P.u0, P.G1 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_2})

# V_291 = Vertex(name = 'V_291',
#                particles = [ P.Sb1__tilde__, P.b0, P.A1 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_52})

# V_292 = Vertex(name = 'V_292',
#                particles = [ P.Sd1__tilde__, P.d0, P.A1 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_52})

# V_293 = Vertex(name = 'V_293',
#                particles = [ P.Ss1__tilde__, P.s0, P.A1 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_52})

# V_294 = Vertex(name = 'V_294',
#                particles = [ P.b0__tilde__, P.Sb1, P.A1 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_52})

# V_295 = Vertex(name = 'V_295',
#                particles = [ P.d0__tilde__, P.Sd1, P.A1 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_52})

# V_296 = Vertex(name = 'V_296',
#                particles = [ P.s0__tilde__, P.Ss1, P.A1 ],
#                color = [ 'Identity(1,2)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_52})

# V_297 = Vertex(name = 'V_297',
#                particles = [ P.Sb1__tilde__, P.b0, P.G1 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_2})

# V_298 = Vertex(name = 'V_298',
#                particles = [ P.Sd1__tilde__, P.d0, P.G1 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_2})

# V_299 = Vertex(name = 'V_299',
#                particles = [ P.Ss1__tilde__, P.s0, P.G1 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_2})

# V_300 = Vertex(name = 'V_300',
#                particles = [ P.b0__tilde__, P.Sb1, P.G1 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_2})

# V_301 = Vertex(name = 'V_301',
#                particles = [ P.d0__tilde__, P.Sd1, P.G1 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_2})

# V_302 = Vertex(name = 'V_302',
#                particles = [ P.s0__tilde__, P.Ss1, P.G1 ],
#                color = [ 'T(3,2,1)' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_2})

# V_303 = Vertex(name = 'V_303',
#                particles = [ P.e1R__plus__, P.e0__minus__, P.A1 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_41})

# V_304 = Vertex(name = 'V_304',
#                particles = [ P.m1R__plus__, P.m0__minus__, P.A1 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_41})

# V_305 = Vertex(name = 'V_305',
#                particles = [ P.t1R__plus__, P.tt0__minus__, P.A1 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_41})

# V_306 = Vertex(name = 'V_306',
#                particles = [ P.e0__plus__, P.e1R__minus__, P.A1 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_41})

# V_307 = Vertex(name = 'V_307',
#                particles = [ P.m0__plus__, P.m1R__minus__, P.A1 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_41})

# V_308 = Vertex(name = 'V_308',
#                particles = [ P.tt0__plus__, P.t1R__minus__, P.A1 ],
#                color = [ '1' ],
#                lorentz = [ L.FFV3 ],
#                couplings = {(0,0):C.GC_41})

V_135 = Vertex(name = 'V_135',
               particles = [ P.e0__plus__, P.e0__minus__, P.A0 ],
               color = [ '1' ],
               lorentz = [ L.FFV],
               couplings = {(0,0):C.GC_41})

V_309 = Vertex(name = 'V_309',
              particles = [P.A0, P.Pion__plus__, P.Pion__minus__],
              color = [ '1' ],
              lorentz = [ L.SSV ],
              couplings = {(0,0):C.GC_41})

V_310 = Vertex(name = 'V_310',
            particles = [P.A0, P.A0, P.Pion__plus__, P.Pion__minus__],
              color = [ '1' ],
              lorentz = [ L.VVSS ],
              couplings = {(0,0):C.GC_68})