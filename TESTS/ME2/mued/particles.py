# This file was automatically created by FeynRules $Revision: 915 $
# Mathematica version: 7.0 for Mac OS X x86 (64-bit) (November 11, 2008)
# Date: Mon 23 Jan 2012 11:06:17


from __future__ import division
from object_library import all_particles, Particle
import parameters as Param

ve0 = Particle(pdg_code = 12,
               name = 've0',
               antiname = 've0~',
               spin = 2,
               color = 1,
               mass = Param.ZERO,
               width = Param.ZERO,
               texname = 've0',
               antitexname = 've0',
               charge = 0,
               GhostNumber = 0,
               LeptonNumber = 1)

ve0__tilde__ = ve0.anti()

vm0 = Particle(pdg_code = 14,
               name = 'vm0',
               antiname = 'vm0~',
               spin = 2,
               color = 1,
               mass = Param.ZERO,
               width = Param.ZERO,
               texname = 'vm0',
               antitexname = 'vm0',
               charge = 0,
               GhostNumber = 0,
               LeptonNumber = 1)

vm0__tilde__ = vm0.anti()

vt0 = Particle(pdg_code = 16,
               name = 'vt0',
               antiname = 'vt0~',
               spin = 2,
               color = 1,
               mass = Param.ZERO,
               width = Param.ZERO,
               texname = 'vt0',
               antitexname = 'vt0',
               charge = 0,
               GhostNumber = 0,
               LeptonNumber = 1)

vt0__tilde__ = vt0.anti()

e0__minus__ = Particle(pdg_code = 11,
                       name = 'e0-',
                       antiname = 'e0+',
                       spin = 2,
                       color = 1,
                       mass = Param.ZERO,
                       width = Param.ZERO,
                       texname = 'e0-',
                       antitexname = 'e0-',
                       charge = -1,
                       GhostNumber = 0,
                       LeptonNumber = 1)

e0__plus__ = e0__minus__.anti()

m0__minus__ = Particle(pdg_code = 13,
                       name = 'm0-',
                       antiname = 'm0+',
                       spin = 2,
                       color = 1,
                       mass = Param.MM,
                       width = Param.ZERO,
                       texname = 'm0-',
                       antitexname = 'm0-',
                       charge = -1,
                       GhostNumber = 0,
                       LeptonNumber = 1)

m0__plus__ = m0__minus__.anti()

tt0__minus__ = Particle(pdg_code = 15,
                        name = 'tt0-',
                        antiname = 'tt0+',
                        spin = 2,
                        color = 1,
                        mass = Param.MTA,
                        width = Param.ZERO,
                        texname = 'tt0-',
                        antitexname = 'tt0-',
                        charge = -1,
                        GhostNumber = 0,
                        LeptonNumber = 1)

tt0__plus__ = tt0__minus__.anti()

u0 = Particle(pdg_code = 2,
              name = 'u0',
              antiname = 'u0~',
              spin = 2,
              color = 3,
              mass = Param.MU,
              width = Param.ZERO,
              texname = 'u0',
              antitexname = 'u0',
              charge = 2/3,
              GhostNumber = 0,
              LeptonNumber = 0)

u0__tilde__ = u0.anti()

c0 = Particle(pdg_code = 4,
              name = 'c0',
              antiname = 'c0~',
              spin = 2,
              color = 3,
              mass = Param.MC,
              width = Param.ZERO,
              texname = 'c0',
              antitexname = 'c0',
              charge = 2/3,
              GhostNumber = 0,
              LeptonNumber = 0)

c0__tilde__ = c0.anti()

t0 = Particle(pdg_code = 6,
              name = 't0',
              antiname = 't0~',
              spin = 2,
              color = 3,
              mass = Param.MT,
              width = Param.WT,
              texname = 't0',
              antitexname = 't0',
              charge = 2/3,
              GhostNumber = 0,
              LeptonNumber = 0)

t0__tilde__ = t0.anti()

d0 = Particle(pdg_code = 1,
              name = 'd0',
              antiname = 'd0~',
              spin = 2,
              color = 3,
              mass = Param.MD,
              width = Param.ZERO,
              texname = 'd0',
              antitexname = 'd0',
              charge = -1/3,
              GhostNumber = 0,
              LeptonNumber = 0)

d0__tilde__ = d0.anti()

s0 = Particle(pdg_code = 3,
              name = 's0',
              antiname = 's0~',
              spin = 2,
              color = 3,
              mass = Param.MS,
              width = Param.ZERO,
              texname = 's0',
              antitexname = 's0',
              charge = -1/3,
              GhostNumber = 0,
              LeptonNumber = 0)

s0__tilde__ = s0.anti()

b0 = Particle(pdg_code = 5,
              name = 'b0',
              antiname = 'b0~',
              spin = 2,
              color = 3,
              mass = Param.MB,
              width = Param.ZERO,
              texname = 'b0',
              antitexname = 'b0',
              charge = -1/3,
              GhostNumber = 0,
              LeptonNumber = 0)

b0__tilde__ = b0.anti()

A0 = Particle(pdg_code = 22,
              name = 'A0',
              antiname = 'A0',
              spin = 3,
              color = 1,
              mass = Param.ZERO,
              width = Param.ZERO,
              texname = 'A0',
              antitexname = 'A0',
              charge = 0,
              GhostNumber = 0,
              LeptonNumber = 0)

Z0 = Particle(pdg_code = 23,
              name = 'Z0',
              antiname = 'Z0',
              spin = 3,
              color = 1,
              mass = Param.MZ,
              width = Param.WZ,
              texname = 'Z0',
              antitexname = 'Z0',
              charge = 0,
              GhostNumber = 0,
              LeptonNumber = 0)

W0__plus__ = Particle(pdg_code = 24,
                      name = 'W0+',
                      antiname = 'W0-',
                      spin = 3,
                      color = 1,
                      mass = Param.MW,
                      width = Param.WW,
                      texname = 'W0+',
                      antitexname = 'W0+',
                      charge = 1,
                      GhostNumber = 0,
                      LeptonNumber = 0)

W0__minus__ = W0__plus__.anti()

G0 = Particle(pdg_code = 21,
              name = 'G0',
              antiname = 'G0',
              spin = 3,
              color = 8,
              mass = Param.ZERO,
              width = Param.ZERO,
              texname = 'G0',
              antitexname = 'G0',
              charge = 0,
              GhostNumber = 0,
              LeptonNumber = 0)

H0 = Particle(pdg_code = 25,
              name = 'H0',
              antiname = 'H0',
              spin = 1,
              color = 1,
              mass = Param.MH,
              width = Param.WH,
              texname = 'H0',
              antitexname = 'H0',
              charge = 0,
              GhostNumber = 0,
              LeptonNumber = 0)

phi0 = Particle(pdg_code = 250,
                name = 'phi0',
                antiname = 'phi0',
                spin = 1,
                color = 1,
                mass = Param.MZ,
                width = Param.ZERO,
                texname = 'phi0',
                antitexname = 'phi0',
                GoldstoneBoson = True,
                charge = 0,
                GhostNumber = 0,
                LeptonNumber = 0)

phi__plus__ = Particle(pdg_code = 251,
                       name = 'phi+',
                       antiname = 'phi-',
                       spin = 1,
                       color = 1,
                       mass = Param.MW,
                       width = Param.ZERO,
                       texname = 'phi+',
                       antitexname = 'phi+',
                       GoldstoneBoson = True,
                       charge = 1,
                       GhostNumber = 0,
                       LeptonNumber = 0)

phi__minus__ = phi__plus__.anti()

ghA = Particle(pdg_code = 9000001,
               name = 'ghA',
               antiname = 'ghA~',
               spin = -1,
               color = 1,
               mass = Param.ZERO,
               width = Param.ZERO,
               texname = 'ghA',
               antitexname = 'ghA',
               charge = 0,
               GhostNumber = 1,
               LeptonNumber = 0)

ghA__tilde__ = ghA.anti()

ghZ = Particle(pdg_code = 9000002,
               name = 'ghZ',
               antiname = 'ghZ~',
               spin = -1,
               color = 1,
               mass = Param.MZ,
               width = Param.ZERO,
               texname = 'ghZ',
               antitexname = 'ghZ',
               charge = 0,
               GhostNumber = 1,
               LeptonNumber = 0)

ghZ__tilde__ = ghZ.anti()

ghWp = Particle(pdg_code = 9000003,
                name = 'ghWp',
                antiname = 'ghWp~',
                spin = -1,
                color = 1,
                mass = Param.MW,
                width = Param.ZERO,
                texname = 'ghWp',
                antitexname = 'ghWp',
                charge = 1,
                GhostNumber = 1,
                LeptonNumber = 0)

ghWp__tilde__ = ghWp.anti()

ghWm = Particle(pdg_code = 9000004,
                name = 'ghWm',
                antiname = 'ghWm~',
                spin = -1,
                color = 1,
                mass = Param.MW,
                width = Param.ZERO,
                texname = 'ghWm',
                antitexname = 'ghWm',
                charge = -1,
                GhostNumber = 1,
                LeptonNumber = 0)

ghWm__tilde__ = ghWm.anti()

ghG = Particle(pdg_code = 9000005,
               name = 'ghG',
               antiname = 'ghG~',
               spin = -1,
               color = 8,
               mass = Param.ZERO,
               width = Param.ZERO,
               texname = 'ghG',
               antitexname = 'ghG',
               charge = 0,
               GhostNumber = 1,
               LeptonNumber = 0)

ghG__tilde__ = ghG.anti()

e1L__minus__ = Particle(pdg_code = 5100011,
                        name = 'e1L-',
                        antiname = 'e1L+',
                        spin = 2,
                        color = 1,
                        mass = Param.DMe,
                        width = Param.DWe,
                        texname = 'e1L-',
                        antitexname = 'e1L-',
                        charge = -1,
                        GhostNumber = 0,
                        LeptonNumber = 1)

e1L__plus__ = e1L__minus__.anti()

m1L__minus__ = Particle(pdg_code = 5100013,
                        name = 'm1L-',
                        antiname = 'm1L+',
                        spin = 2,
                        color = 1,
                        mass = Param.DMmu,
                        width = Param.DWmu,
                        texname = 'm1L-',
                        antitexname = 'm1L-',
                        charge = -1,
                        GhostNumber = 0,
                        LeptonNumber = 1)

m1L__plus__ = m1L__minus__.anti()

t1L__minus__ = Particle(pdg_code = 5100015,
                        name = 't1L-',
                        antiname = 't1L+',
                        spin = 2,
                        color = 1,
                        mass = Param.DMtau,
                        width = Param.DWtau,
                        texname = 't1L-',
                        antitexname = 't1L-',
                        charge = -1,
                        GhostNumber = 0,
                        LeptonNumber = 1)

t1L__plus__ = t1L__minus__.anti()

ve1 = Particle(pdg_code = 5100012,
               name = 've1',
               antiname = 've1~',
               spin = 2,
               color = 1,
               mass = Param.DMn1,
               width = Param.DWn1,
               texname = 've1',
               antitexname = 've1',
               charge = 0,
               GhostNumber = 0,
               LeptonNumber = 1)

ve1__tilde__ = ve1.anti()

vm1 = Particle(pdg_code = 5100014,
               name = 'vm1',
               antiname = 'vm1~',
               spin = 2,
               color = 1,
               mass = Param.DMn2,
               width = Param.DWn2,
               texname = 'vm1',
               antitexname = 'vm1',
               charge = 0,
               GhostNumber = 0,
               LeptonNumber = 1)

vm1__tilde__ = vm1.anti()

vt1 = Particle(pdg_code = 5100016,
               name = 'vt1',
               antiname = 'vt1~',
               spin = 2,
               color = 1,
               mass = Param.DMn3,
               width = Param.DWn3,
               texname = 'vt1',
               antitexname = 'vt1',
               charge = 0,
               GhostNumber = 0,
               LeptonNumber = 1)

vt1__tilde__ = vt1.anti()

e1R__minus__ = Particle(pdg_code = 6100011,
                        name = 'e1R-',
                        antiname = 'e1R+',
                        spin = 2,
                        color = 1,
                        mass = Param.SMe,
                        width = Param.SWe,
                        texname = 'e1R-',
                        antitexname = 'e1R-',
                        charge = -1,
                        GhostNumber = 0,
                        LeptonNumber = 1)

e1R__plus__ = e1R__minus__.anti()

m1R__minus__ = Particle(pdg_code = 6100013,
                        name = 'm1R-',
                        antiname = 'm1R+',
                        spin = 2,
                        color = 1,
                        mass = Param.SMmu,
                        width = Param.SWmu,
                        texname = 'm1R-',
                        antitexname = 'm1R-',
                        charge = -1,
                        GhostNumber = 0,
                        LeptonNumber = 1)

m1R__plus__ = m1R__minus__.anti()

t1R__minus__ = Particle(pdg_code = 6100015,
                        name = 't1R-',
                        antiname = 't1R+',
                        spin = 2,
                        color = 1,
                        mass = Param.SMtau,
                        width = Param.SWtau,
                        texname = 't1R-',
                        antitexname = 't1R-',
                        charge = -1,
                        GhostNumber = 0,
                        LeptonNumber = 1)

t1R__plus__ = t1R__minus__.anti()

Su1 = Particle(pdg_code = 6100002,
               name = 'Su1',
               antiname = 'Su1~',
               spin = 2,
               color = 3,
               mass = Param.SMu,
               width = Param.SWu,
               texname = 'Su1',
               antitexname = 'Su1',
               charge = 2/3,
               GhostNumber = 0,
               LeptonNumber = 0)

Su1__tilde__ = Su1.anti()

Sc1 = Particle(pdg_code = 6100004,
               name = 'Sc1',
               antiname = 'Sc1~',
               spin = 2,
               color = 3,
               mass = Param.SMc,
               width = Param.SWc,
               texname = 'Sc1',
               antitexname = 'Sc1',
               charge = 2/3,
               GhostNumber = 0,
               LeptonNumber = 0)

Sc1__tilde__ = Sc1.anti()

St1 = Particle(pdg_code = 6100006,
               name = 'St1',
               antiname = 'St1~',
               spin = 2,
               color = 3,
               mass = Param.SMtop,
               width = Param.ZERO,
               texname = 'St1',
               antitexname = 'St1',
               charge = 2/3,
               GhostNumber = 0,
               LeptonNumber = 0)

St1__tilde__ = St1.anti()

Sd1 = Particle(pdg_code = 6100001,
               name = 'Sd1',
               antiname = 'Sd1~',
               spin = 2,
               color = 3,
               mass = Param.SMd,
               width = Param.SWd,
               texname = 'Sd1',
               antitexname = 'Sd1',
               charge = -1/3,
               GhostNumber = 0,
               LeptonNumber = 0)

Sd1__tilde__ = Sd1.anti()

Ss1 = Particle(pdg_code = 6100003,
               name = 'Ss1',
               antiname = 'Ss1~',
               spin = 2,
               color = 3,
               mass = Param.SMs,
               width = Param.SWs,
               texname = 'Ss1',
               antitexname = 'Ss1',
               charge = -1/3,
               GhostNumber = 0,
               LeptonNumber = 0)

Ss1__tilde__ = Ss1.anti()

Sb1 = Particle(pdg_code = 6100005,
               name = 'Sb1',
               antiname = 'Sb1~',
               spin = 2,
               color = 3,
               mass = Param.SMb,
               width = Param.SWb,
               texname = 'Sb1',
               antitexname = 'Sb1',
               charge = -1/3,
               GhostNumber = 0,
               LeptonNumber = 0)

Sb1__tilde__ = Sb1.anti()

Du1 = Particle(pdg_code = 5100002,
               name = 'Du1',
               antiname = 'Du1~',
               spin = 2,
               color = 3,
               mass = Param.DMu,
               width = Param.DWu,
               texname = 'Du1',
               antitexname = 'Du1',
               charge = 2/3,
               GhostNumber = 0,
               LeptonNumber = 0)

Du1__tilde__ = Du1.anti()

Dc1 = Particle(pdg_code = 5100004,
               name = 'Dc1',
               antiname = 'Dc1~',
               spin = 2,
               color = 3,
               mass = Param.DMc,
               width = Param.DWc,
               texname = 'Dc1',
               antitexname = 'Dc1',
               charge = 2/3,
               GhostNumber = 0,
               LeptonNumber = 0)

Dc1__tilde__ = Dc1.anti()

Dt1 = Particle(pdg_code = 5100006,
               name = 'Dt1',
               antiname = 'Dt1~',
               spin = 2,
               color = 3,
               mass = Param.DMtop,
               width = Param.DWtop,
               texname = 'Dt1',
               antitexname = 'Dt1',
               charge = 2/3,
               GhostNumber = 0,
               LeptonNumber = 0)

Dt1__tilde__ = Dt1.anti()

Dd1 = Particle(pdg_code = 5100001,
               name = 'Dd1',
               antiname = 'Dd1~',
               spin = 2,
               color = 3,
               mass = Param.DMd,
               width = Param.DWd,
               texname = 'Dd1',
               antitexname = 'Dd1',
               charge = -1/3,
               GhostNumber = 0,
               LeptonNumber = 0)

Dd1__tilde__ = Dd1.anti()

Ds1 = Particle(pdg_code = 5100003,
               name = 'Ds1',
               antiname = 'Ds1~',
               spin = 2,
               color = 3,
               mass = Param.DMs,
               width = Param.DWs,
               texname = 'Ds1',
               antitexname = 'Ds1',
               charge = -1/3,
               GhostNumber = 0,
               LeptonNumber = 0)

Ds1__tilde__ = Ds1.anti()

Db1 = Particle(pdg_code = 5100005,
               name = 'Db1',
               antiname = 'Db1~',
               spin = 2,
               color = 3,
               mass = Param.DMb,
               width = Param.DWb,
               texname = 'Db1',
               antitexname = 'Db1',
               charge = -1/3,
               GhostNumber = 0,
               LeptonNumber = 0)

Db1__tilde__ = Db1.anti()

A1 = Particle(pdg_code = 5100022,
              name = 'A1',
              antiname = 'A1',
              spin = 3,
              color = 1,
              mass = Param.MB1,
              width = Param.ZERO,
              texname = 'A1',
              antitexname = 'A1',
              charge = 0,
              GhostNumber = 0,
              LeptonNumber = 0)

Z1 = Particle(pdg_code = 5100023,
              name = 'Z1',
              antiname = 'Z1',
              spin = 3,
              color = 1,
              mass = Param.MZ1,
              width = Param.WZ1,
              texname = 'Z1',
              antitexname = 'Z1',
              charge = 0,
              GhostNumber = 0,
              LeptonNumber = 0)

W1__plus__ = Particle(pdg_code = 5100024,
                      name = 'W1+',
                      antiname = 'W1-',
                      spin = 3,
                      color = 1,
                      mass = Param.MW1,
                      width = Param.WW1,
                      texname = 'W1+',
                      antitexname = 'W1+',
                      charge = 1,
                      GhostNumber = 0,
                      LeptonNumber = 0)

W1__minus__ = W1__plus__.anti()

G1 = Particle(pdg_code = 5100021,
              name = 'G1',
              antiname = 'G1',
              spin = 3,
              color = 8,
              mass = Param.MG1,
              width = Param.WG1,
              texname = 'G1',
              antitexname = 'G1',
              charge = 0,
              GhostNumber = 0,
              LeptonNumber = 0)

H1 = Particle(pdg_code = 5100030,
              name = 'H1',
              antiname = 'H1',
              spin = 1,
              color = 1,
              mass = Param.MH1,
              width = Param.WH1,
              texname = 'H1',
              antitexname = 'H1',
              charge = 0,
              GhostNumber = 0,
              LeptonNumber = 0)

