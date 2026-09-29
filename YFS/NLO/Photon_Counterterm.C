#include "YFS/NLO/Photon_Counterterm.H"

#include <cmath>
#include <functional>
#include <sstream>
#include <iomanip>
#include <vector>

/*
  The external on-shell photon's charge counterterm, as OpenLoops
  computes it (source: OpenLoops lib_src/openloops/src/renormalisation_ew.F90;
  line numbers below are that file's). With Sherpa's interface OpenLoops runs
  with ew_scheme 2 (the alpha Sherpa passes is taken as alpha(M_Z) - for a
  G_mu card that is alpha_Gmu), ew_renorm_scheme = Sherpa's EW_REN_SCHEME,
  the complex-mass scheme (cms_on 1), massless u, d, s, c, and whatever b mass
  Sherpa passes (0 for a massless b). For a photon registered as pdg 2002
  photon_factors (2711-2764) multiplies the tree by alpha(0)/alpha_in and
  adds to the loop

    loopfactor = [dZe0QEDEWnreg - dZeQEDEW] alpha_in/(4 pi)          (2750)

  in units of the tree, i.e. 2 x loopfactor on V_fin/T. With

    dZe0QEDEWnreg = -dZAAEWnreg/2 - (sw/cw) SiAZ0/MZ2                 (1496)
    dZAAEWnreg    = -(dSiAAheavy0 + PiAAlightZ/Re MZ2 + dAlpha)       (1283)
    dAlpha        = (4 pi/alpha_in)(1 - alpha(0)/alpha(MZ)),
                    alpha(MZ) = alpha_in in ew_scheme 2                (1280)
    dZeGmuQEDEW   = dsw/sw - SiAZ0/(sw cw MZ2)
                    - (6 + (7 - 4 sw2)/(2 sw2) ln cw2)/(2 sw2)
                    + (dZMW2 - SiW0)/(2 MW2)                           (1500-1502)
    dZeZQEDEW     = -(dZAAEWnreg + dAlpha)/2 - (sw/cw) SiAZ0/MZ2       (1505)
    dsw/sw = -(cw2/sw2) dcw/cw,  dcw = (cw/2)(dZMW2/MW2 - dZMZ2/MZ2)  (1487-1488)
    CMS: dZMW2 = SiW + (MW2 - Re MW2) dSiW + 4 (Re MW2 - MW2),
         dZMZ2 = SiZZ - (Re MZ2 - MZ2) dSiZZ                           (1297-1320)

  and the self-energies SiW(Re MW2), dSiW, SiW0 = SiW(0), SiZZ(Re MZ2),
  dSiZZ, SiAZ0 (bosonic 690-749, fermionic 1100-1257; Denner92 in units of
  alpha/4pi, 't Hooft-Feynman gauge), B0 with Delta_UV = 0 at real p^2
  (calcB0: p2 = real(p2)) and, for calcRB0/calcRdB0, the real part above
  the threshold Re p^2 > Re(m1^2 + m2^2). The combination is UV finite and
  mu independent (checked). In the alpha(0) renormalisation scheme the
  photon needs nothing (c = 0); in the alpha(M_Z) scheme
  c = 1 - alpha(0)/alpha_in; in the G_mu scheme it is Delta r with the
  light-fermion Delta alpha replaced by 1 - alpha(0)/alpha_in: at the Z pole
  (alpha_in = 1/131.931) c = 0.0075211, against 0.0075209 from the soft-limit
  calibration and 0.00752 from pyol's 2002 - 22 difference (NOTES-yfsnlo-
  realvirtual-2026-09-28.md, sec. 8).
*/

using namespace YFS;

namespace {

  typedef std::complex<double> C;
  inline double sqr(double x) { return x*x; }

  // adaptive Gauss-Kronrod (7-15) for complex integrands on [a,b]
  const double s_xk[8] = {0.991455371120812639206854697526329,
                          0.949107912342758524526189684047851,
                          0.864864423359769072789712788640926,
                          0.741531185599394439863864773280788,
                          0.586087235467691130294144845693013,
                          0.405845151377397166906606412076961,
                          0.207784955007898467600689403773245,
                          0.000000000000000000000000000000000};
  const double s_wk[8] = {0.022935322010529224963732008058970,
                          0.063092092629978553290700663189204,
                          0.104790010322250183839876322541518,
                          0.140653259715525918745189590510238,
                          0.169004726639267902826583426598550,
                          0.190350578064785409913256402421014,
                          0.204432940075298892414161999234649,
                          0.209482141084727828012999174891714};
  const double s_wg[4] = {0.129484966168869693270611432679082,
                          0.279705391489276667901467771423780,
                          0.381830050505118944950369775488975,
                          0.417959183673469387755102040816327};

  C GK15(const std::function<C(double)> &f, double a, double b, double &err)
  {
    const double h(0.5*(b - a)), m(0.5*(b + a));
    C k(s_wk[7]*f(m)), g(s_wg[3]*f(m));
    for (int j(0); j < 7; ++j) {
      const C fp(f(m + h*s_xk[j])), fm(f(m - h*s_xk[j]));
      k += s_wk[j]*(fp + fm);
      if (j % 2 == 1) g += s_wg[j/2]*(fp + fm);
    }
    err = std::abs(h*(k - g));
    return h*k;
  }

  C Adaptive(const std::function<C(double)> &f, double a, double b,
             double tol, int depth)
  {
    double err(0.);
    const C r(GK15(f, a, b, err));
    if (err < tol || depth > 40) return r;
    const double m(0.5*(a + b));
    return Adaptive(f, a, m, 0.5*tol, depth + 1) + Adaptive(f, m, b, 0.5*tol, depth + 1);
  }

  //! integral over [0,1], pre-split where the integrands peak
  C Integrate01(const std::function<C(double)> &f)
  {
    const double pts[] = {0., 1e-6, 1e-4, 1e-2, 0.5, 1. - 1e-2, 1. - 1e-4, 1. - 1e-6, 1.};
    C r(0.);
    for (int i(0); i < 8; ++i) r += Adaptive(f, pts[i], pts[i+1], 1e-15, 0);
    return r;
  }

  struct Loops {
    double mu2;
    //! calcB0: Delta_UV = 0, p2 -> Re p2; B0(0;0,0) = 0 (scaleless)
    C B0(C p2in, C m12, C m22) const {
      const double p2(p2in.real());
      if (p2 == 0.) {
        if (m12 == 0. && m22 == 0.) return 0.;
        if (m12 == m22) return -std::log(m12/mu2);
        if (m12 == 0.) return 1. - std::log(m22/mu2);
        if (m22 == 0.) return 1. - std::log(m12/mu2);
        return 1. - (m12*std::log(m12/mu2) - m22*std::log(m22/mu2))/(m12 - m22);
      }
      return Integrate01([&](double x) {
        return -std::log((x*m22 + (1.-x)*m12 - x*(1.-x)*p2 - C(0., 1e-300))/mu2); });
    }
    //! calcdB0 = d B0/d p2
    C dB0(C p2in, C m12, C m22) const {
      const double p2(p2in.real());
      if (p2 == 0. && m12 == m22) return 1./(6.*m12);
      return Integrate01([&](double x) {
        return x*(1.-x)/(x*m22 + (1.-x)*m12 - x*(1.-x)*p2 - C(0., 1e-300)); });
    }
    //! calcRB0 / calcRdB0: real part above Re(m1^2 + m2^2)
    static bool Above(C p2, C m12, C m22) { return p2.imag() == 0. && p2.real() > (m12 + m22).real(); }
    C RB0(C p2, C m12, C m22) const { const C r(B0(p2, m12, m22)); return Above(p2, m12, m22) ? C(r.real(), 0.) : r; }
    C RdB0(C p2, C m12, C m22) const { const C r(dB0(p2, m12, m22)); return Above(p2, m12, m22) ? C(r.real(), 0.) : r; }
  };

  C Cplx(double m, double w, bool cms) { return cms ? C(m*m, -m*w) : C(m*m, 0.); }

  struct Lepton { double m2; C B00, B0Z, B000, dB000, B0W0, dB0W0, dB0Z; };
}

double YFS::OnShellPhotonCounterterm(const Photon_CT_Input &in,
                                     photon_ct_scheme scheme, std::string &log)
{
  if (scheme == photon_ct_scheme::alpha0) { log = "alpha(0) scheme: no counterterm"; return 0.; }
  if (scheme == photon_ct_scheme::alphamZ) {
    const double c(1. - in.alpha0/in.alpha_in);
    std::ostringstream o;
    o<<"alpha(M_Z) scheme: c = 1 - alpha(0)/alpha_in = "<<c;
    log = o.str();
    return c;
  }
  const Loops L{sqr(in.MZ)};
  const double nc(3.), pi(M_PI);
  const C MZ2(Cplx(in.MZ, in.WZ, in.cms)), MW2(Cplx(in.MW, in.WW, in.cms));
  const C MH2(Cplx(in.MH, in.WH, in.cms)), MT2(Cplx(in.MT, in.WT, in.cms));
  const double rMZ2(sqr(in.MZ)), rMW2(sqr(in.MW));
  const double MB2(sqr(in.mb));
  const C cw(std::sqrt(MW2)/std::sqrt(MZ2)), cw2(cw*cw), sw2(1. - cw2), sw(std::sqrt(sw2)),
          cw4(cw2*cw2);
  const C gZRH(-sw/cw), gZLH(1./(sw*cw));
  // right (0) and left (1) handed Z couplings, parameters_init.F90:412-415
  const C gZn[2] = {0., gZLH*0.5};
  const C gZl[2] = {-gZRH, gZLH*(-0.5 + sw2)};
  const C gZu[2] = {2.*gZRH/3., gZLH*(0.5 - 2.*sw2/3.)};
  const C gZd[2] = {-gZRH/3., gZLH*(-0.5 + sw2/3.)};
  const C gn2(gZn[0]*gZn[0] + gZn[1]*gZn[1]), gl2(gZl[0]*gZl[0] + gZl[1]*gZl[1]),
          gu2(gZu[0]*gZu[0] + gZu[1]*gZu[1]), gd2(gZd[0]*gZd[0] + gZd[1]*gZd[1]);
  const C Z(0.);
  const C B00WW(L.B0(Z,MW2,MW2)), dB00WW(L.dB0(Z,MW2,MW2)), B00W0(L.B0(Z,MW2,Z)), dB00W0(L.dB0(Z,MW2,Z));
  const C B00HH(L.B0(Z,MH2,MH2)), B00ZZ(L.B0(Z,MZ2,MZ2)), B00ZH(L.B0(Z,MZ2,MH2)), B00WH(L.B0(Z,MW2,MH2)),
          dB00WH(L.dB0(Z,MW2,MH2)), B00WZ(L.B0(Z,MW2,MZ2)), dB00WZ(L.dB0(Z,MW2,MZ2));
  const C B00TT(L.B0(Z,MT2,MT2)), dB00TT(L.dB0(Z,MT2,MT2)), B00BB(L.B0(Z,MB2,MB2)),
          B00TB(L.B0(Z,MT2,MB2)), dB00TB(L.dB0(Z,MT2,MB2));
  const C B0WW0(L.B0(MW2,MW2,Z)), dB0WW0(L.dB0(MW2,MW2,Z)), B0ZZH(L.B0(MZ2,MZ2,MH2)), dB0ZZH(L.dB0(MZ2,MZ2,MH2));
  const C B0WWH(L.B0(MW2,MW2,MH2)), dB0WWH(L.dB0(MW2,MW2,MH2)), B0WWZ(L.B0(MW2,MW2,MZ2)), dB0WWZ(L.dB0(MW2,MW2,MZ2));
  const C B0ZWW(L.RB0(MZ2,MW2,MW2)), dB0ZWW(L.RdB0(MZ2,MW2,MW2));
  const C B0WTB(L.RB0(MW2,MT2,MB2)), dB0WTB(L.RdB0(MW2,MT2,MB2)), B0ZTT(L.RB0(MZ2,MT2,MT2)), dB0ZTT(L.RdB0(MZ2,MT2,MT2));
  const C B0Z00(L.RB0(MZ2,Z,Z)), dB0Z00(L.RdB0(MZ2,Z,Z)), B0W00(L.RB0(MW2,Z,Z)), dB0W00(L.RdB0(MW2,Z,Z));
  const C B0ZBB(L.RB0(MZ2,MB2,MB2)), dB0ZBB(L.RdB0(MZ2,MB2,MB2));
  std::vector<Lepton> lep;
  for (double m : {in.me, in.mmu, in.mtau}) {
    const double m2(m*m);
    lep.push_back(Lepton{m2, L.B0(Z,m2,m2), L.RB0(MZ2,m2,m2), L.B0(Z,Z,m2), L.dB0(Z,Z,m2),
                         L.RB0(MW2,Z,m2), L.RdB0(MW2,Z,m2), L.RdB0(MZ2,m2,m2)});
  }
  // bosonic, 692-749
  C dSiAAheavy0(-3.*B00WW - 4.*MW2*dB00WW);
  const C SiAZ0(2.*MW2*B00WW/(cw*sw));
  C SiZZ(-(((-1. + 4.*cw2)*rMZ2)/3. - (2. - 8.*cw2 + 24.*cw4)*MW2*B00WW
           + ((-10. + 16.*cw2 + 24.*cw4)*MW2 + (-0.5 + 2.*cw2 + 18.*cw4)*rMZ2)*B0ZWW)/(6.*cw2*sw2)
         - ((-2.*rMZ2)/3. - 2.*MH2*B00HH - 2.*MZ2*B00ZZ + (2.*MH2 - 10.*MZ2 - rMZ2)*B0ZZH
            - ((MZ2 - MH2)*(MZ2 - MH2)*(-B00ZH + B0ZZH))/rMZ2)/(12.*cw2*sw2));
  C dSiZZ((2./3. + ((MH2 - MZ2)*(MH2 - MZ2)*(B00ZH - B0ZZH))/(rMZ2*rMZ2) + B0ZZH
           + (2. - 8.*cw2 - 3.*(-1. + 4.*cw2 + 36.*cw4)*B0ZWW
              - 3.*(4.*(-5. + 8.*cw2 + 12.*cw4)*MW2 + (-1. + 4.*cw2 + 36.*cw4)*rMZ2)*dB0ZWW)/3.
           - (2.*MH2 - 10.*MZ2 - rMZ2)*dB0ZZH + ((MH2 - MZ2)*(MH2 - MZ2)*dB0ZZH)/rMZ2)/(12.*cw2*sw2));
  C SiW((-8.*(rMW2/3. - 2.*MW2*B00WW + (MW2*MW2*(B00W0 - B0WW0))/rMW2 + (2.*MW2 + 5.*rMW2)*B0WW0)
         - ((-2.*rMW2)/3. - 2.*MH2*B00HH - 2.*MW2*B00WW + ((MH2 - MW2)*(MH2 - MW2)*(B00WH - B0WWH))/rMW2
            + (2.*MH2 - 10.*MW2 - rMW2)*B0WWH)/sw2
         - ((2.*(-1. + 4.*cw2)*rMW2)/3. - 2.*(1. + 8.*cw2)*(MW2*B00WW + MZ2*B00ZZ)
            + ((1. + 8.*cw2)*(MW2 - MZ2)*(MW2 - MZ2)*(B00WZ - B0WWZ))/rMW2
            + ((54. - 10./cw2 + 16.*cw2)*MW2 + (-1. + 40.*cw2)*rMW2)*B0WWZ)/sw2)/12.);
  C dSiW((-2.*(1./3. + 5.*B0WW0 + (MW2*MW2*(-B00W0 + B0WW0))/(rMW2*rMW2) - (MW2*MW2*dB0WW0)/rMW2
               + (2.*MW2 + 5.*rMW2)*dB0WW0))/3.
         - (-2./3. - B0WWH + ((MW2 - MH2)*(MW2 - MH2)*(-B00WH + B0WWH))/(rMW2*rMW2)
            + (2.*MH2 - 10.*MW2 - rMW2)*dB0WWH - ((MW2 - MH2)*(MW2 - MH2)*dB0WWH)/rMW2)/(12.*sw2)
         - ((2.*(-1. + 4.*cw2))/3. + (-1. + 40.*cw2)*B0WWZ
            + ((1. + 8.*cw2)*(MW2 - MZ2)*(MW2 - MZ2)*(-B00WZ + B0WWZ))/(rMW2*rMW2)
            - ((1. + 8.*cw2)*(MW2 - MZ2)*(MW2 - MZ2)*dB0WWZ)/rMW2
            + ((54. - 10./cw2 + 16.*cw2)*MW2 + (-1. + 40.*cw2)*rMW2)*dB0WWZ)/(12.*sw2));
  C SiW0((-8.*(2.*MW2*B00W0 - 2.*MW2*B00WW - MW2*MW2*dB00W0)
          - (-2.*MH2*B00HH + (2.*MH2 - 10.*MW2)*B00WH - 2.*MW2*B00WW - (MH2 - MW2)*(MH2 - MW2)*dB00WH)/sw2
          - ((54. - 10./cw2 + 16.*cw2)*MW2*B00WZ - 2.*(1. + 8.*cw2)*(MW2*B00WW + MZ2*B00ZZ)
             - (1. + 8.*cw2)*(MW2 - MZ2)*(MW2 - MZ2)*dB00WZ)/sw2)/12.);
  // fermionic, 1104-1257
  C PiAAlightZ(0.);
  for (const Lepton &l : lep)
    PiAAlightZ += -4.*(rMZ2*(1./3. - l.B0Z) - 2.*l.m2*(l.B0Z - l.B00))/3.;
  PiAAlightZ += -2.*20.*nc*rMZ2*(1./3. - B0Z00)/27.;
  PiAAlightZ += -4.*nc*(rMZ2*(1./3. - B0ZBB) - 2.*MB2*(B0ZBB - B00BB))/27.;
  dSiAAheavy0 += -(16.*nc*(1. - 3.*B00TT - 6.*MT2*dB00TT))/81.;
  for (const Lepton &l : lep) {
    SiZZ += 2./3.*(rMZ2*(-1. + 3.*B0Z00)*gn2/3.);
    SiZZ += -2./3.*((3.*l.m2*l.B0Z)/(4.*cw2*sw2) + (rMZ2/3. + 2.*l.m2*l.B00 - (2.*l.m2 + rMZ2)*l.B0Z)*gl2);
    dSiZZ += -2./3.*((3.*l.m2*l.dB0Z)/(4.*cw2*sw2) - ((-1. + 3.*l.B0Z + 3.*(2.*l.m2 + rMZ2)*l.dB0Z)*gl2)/3.);
    dSiZZ += 2./3.*((-1. + 3.*B0Z00 + 3.*rMZ2*dB0Z00)*gn2/3.);
    SiW += -(rMW2/3. + l.m2*l.B00 + ((l.m2 - 2.*rMW2)*l.B0W0)/2. + (l.m2*l.m2*(-l.B000 + l.B0W0))/(2.*rMW2))/(3.*sw2);
    dSiW += -(0.5*(2./3. + l.m2*l.m2*(l.B000 - l.B0W0)/(rMW2*rMW2) - 2.*l.B0W0 + l.m2*l.m2*l.dB0W0/rMW2
                  + (l.m2 - 2.*rMW2)*l.dB0W0))/(3.*sw2);
    SiW0 += -1./(3.*sw2)*(l.m2*l.B00 + l.m2*l.B000/2. + l.m2*l.m2*l.dB000/2.);
  }
  SiZZ += -2./3.*(-2.*nc*rMZ2*(-1. + 3.*B0Z00)*gd2/3.);
  SiZZ += -2./3.*(-2.*nc*rMZ2*(-1. + 3.*B0Z00)*gu2/3.);
  SiZZ += -2./3.*(nc*((3.*MB2*B0ZBB)/(4.*cw2*sw2) + (rMZ2/3. + 2.*MB2*B00BB - (2.*MB2 + rMZ2)*B0ZBB)*gd2));
  SiZZ += -2./3.*(nc*((3.*MT2*B0ZTT)/(4.*cw2*sw2) + (rMZ2/3. + 2.*MT2*B00TT - (2.*MT2 + rMZ2)*B0ZTT)*gu2));
  dSiZZ += -2./3.*(-2.*nc*(-1. + 3.*B0Z00 + 3.*rMZ2*dB0Z00)*gd2/3.);
  dSiZZ += -2./3.*(-2.*nc*(-1. + 3.*B0Z00 + 3.*rMZ2*dB0Z00)*gu2/3.);
  dSiZZ += -2./3.*nc*((3.*MB2*dB0ZBB)/(4.*cw2*sw2) - ((-1. + 3.*B0ZBB + 3.*(2.*MB2 + rMZ2)*dB0ZBB)*gd2)/3.);
  dSiZZ += -2./3.*nc*((3.*MT2*dB0ZTT)/(4.*cw2*sw2) - ((-1. + 3.*B0ZTT + 3.*(2.*MT2 + rMZ2)*dB0ZTT)*gu2)/3.);
  SiW += -(2.*nc*rMW2*(1./3. - B0W00))/(3.*sw2);
  SiW += -(nc*(rMW2/3. + MB2*B00BB + MT2*B00TT + ((MB2 + MT2 - 2.*rMW2)*B0WTB)/2.
               + ((MB2 - MT2)*(MB2 - MT2)*(-B00TB + B0WTB))/(2.*rMW2)))/(3.*sw2);
  dSiW += -(2.*nc*(1./3. - B0W00 - rMW2*dB0W00))/(3.*sw2);
  dSiW += -(nc/2.*(2./3. + ((MB2 - MT2)*(MB2 - MT2)*(B00TB - B0WTB))/(rMW2*rMW2) - 2.*B0WTB
                   + ((MB2 - MT2)*(MB2 - MT2)*dB0WTB)/rMW2 + (MB2 + MT2 - 2.*rMW2)*dB0WTB))/(3.*sw2);
  SiW0 += -nc/(3.*sw2)*(MB2*B00BB + ((MB2 + MT2)*B00TB)/2. + MT2*B00TT + ((MB2 - MT2)*(MB2 - MT2)*dB00TB)/2.);
  // mass RCs (CMS), 1297-1320; weak mixing angle, 1487-1488
  const C dZMW2(in.cms ? SiW + (MW2 - rMW2)*dSiW + 4.*(rMW2 - MW2) : C(SiW.real(), 0.));
  const C dZMZ2(in.cms ? SiZZ - (rMZ2 - MZ2)*dSiZZ : C(SiZZ.real(), 0.));
  const C dcw(cw/2.*(dZMW2/MW2 - dZMZ2/MZ2)), dsw(-cw/sw*dcw);
  // charge RCs, 1280-1283, 1496-1502
  const double dAlpha(4.*pi/in.alpha_in*(1. - in.alpha0/in.alpha_in));
  const C dZAAnreg(-(dSiAAheavy0 + PiAAlightZ/rMZ2 + dAlpha));
  const C dZe0(-0.5*dZAAnreg - sw/cw*SiAZ0/MZ2);
  const C dZeGmu(dsw/sw - 1./sw/cw*SiAZ0/MZ2 - 0.5/sw2*(6. + (7. - 4.*sw2)/(2.*sw2)*std::log(cw2))
                 + 0.5*(dZMW2 - SiW0)/MW2);
  const double c(2.*(dZe0.real() - dZeGmu.real())*in.alpha_in/(4.*pi));
  std::ostringstream o;
  o<<std::setprecision(8)<<"G_mu scheme: c = "<<c<<" (1 - alpha(0)/alpha_in = "
   <<1. - in.alpha0/in.alpha_in<<", the rest of Delta r "<<c - (1. - in.alpha0/in.alpha_in)<<")";
  log = o.str();
  return c;
}
