#include "YFS/Main/YFS_Form_Factor.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#include "ATOOLS/Math/Poincare.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Org/Message.H"
#include "MODEL/Main/Model_Base.H"



using namespace YFS;
using namespace ATOOLS;
using namespace MODEL;

YFS_Form_Factor::YFS_Form_Factor() {
}

YFS_Form_Factor::~YFS_Form_Factor() {

}

double YFS_Form_Factor::BVR_full(double p1p2, double E1, double E2,
                                 double Mas1, double Mas2, double Kmax, double MasPhot, int mode) {
  double alpi = m_alpi;
  double m12 = Mas1 * Mas2;
  double t1, t2, AA4;
  if (p1p2 - m12 < 1e-10) return 0;
  double beta1 = sqrt(1. - sqr(Mas1 / E1));
  double beta2 = sqrt(1. - sqr(Mas2 / E2));
  double rho = sqrt(1 - sqr(m12 / p1p2));
  t1 = (p1p2 * A(p1p2, Mas1, Mas2) - 1) * log(4 * sqr(Kmax / MasPhot));
  if ( mode == 0 ) {
    t2 = p1p2 * A4(p1p2, E1, E2, Mas1, Mas2);
  }
  else {
    // alpi = 1./137.03599976000001/M_PI;
    if (sqr(Mas1) / p1p2 > 1e-10 && sqr(Mas2) / p1p2 > 1e-10  ) {
      AA4 = A4(p1p2, E1, E2, Mas1, Mas2);
      t2 = AA4 * p1p2;
      // t1 = (p1p2*A(p1p2,Mas1,Mas2)-1.)*log(4*(E1*E2/Mas1/Mas2));
    }
    else {
      AA4 = A4(p1p2, E1, E2, Mas1, Mas2);
      t2 = AA4 * p1p2;
    }
  }
  double t3 = Mas1 * Mas1 * A4_eq(E1, Mas1) + Mas2 * Mas2 * A4_eq(E2, Mas2);
  if (IsBad(t1) || IsBad(t2) || IsBad(t3)) {
    msg_Error() << METHOD << "\n"
                << "YFS Form Factor is NaN"
                << "\n T1    = " << t1
                << "\n T2    = " << p1p2*AA4
                << "\n T3    = " << t3 * 0.5
                << "\n E1    = " << E1
                << "\n E2    = " << E2
                << "\n Mass1 = " << Mas1
                << "\n Mass2 = " << Mas2
                << "\n Kmax = " << Kmax
                << "\n MasPhot = " << MasPhot
                << "\n M12 = " << m12
                << "\n A4 = " << AA4
                << "\n p1p2  = " << p1p2
                << "\n form   = " << exp(m_alpi * (t1 + t2 + 0.5 * t3)) << std::endl;
  }
  return m_alpi * (t1 + t2 + 0.5 * t3);
}

double YFS_Form_Factor::BVR_full(Vec4D p1, Vec4D p2,  double Kmax, double MasPhot, int mode) {
  return BVR_full(p1 * p2, p1.E(), p2.E(), p1.Mass(), p2.Mass(), Kmax, MasPhot, mode);
}

double YFS_Form_Factor::BVR_cru(double p1p2, double E1, double E2,
                                double Mas1, double Mas2, double Kmax, double MasPhot) {
  // m_btilcru = m_alpi*(p1p2*BVR_A(p1p2,Mas1,Mas2))
  double t1 = (p1p2 * A(p1p2, Mas1, Mas2)) * log(4.*sqr(Kmax / MasPhot));
  double t2 = p1p2 * A4(p1p2, E1, E2, Mas1, Mas2);
  if (IsBad(t1) || IsBad(t2)) {
    msg_Error() << METHOD << "\n" << "YFS Form Factor is NaN"
                << "\n T1    = " << t1
                << "\n T2    = " << t2
                << "\n E1    = " << E1
                << "\n E2    = " << E2
                << "\n Mass1 = " << Mas1
                << "\n Mass2 = " << Mas2
                << "\n p1p2  = " << p1p2 << std::endl;
  }
  return m_alpi * (t1 + t2);
}


double YFS_Form_Factor::A(double p1p2, double m1, double m2) {
  double m12 = m1 * m2;
  if ((p1p2 - m12) < 1e-10) return 0;
  double xlam = sqrt((p1p2 - m12) * (p1p2 + m12));
  double BVR_A = 1. / xlam * log((p1p2 + xlam) / m12);
  return BVR_A;
}

double YFS_Form_Factor::YijEta(double eta, double y1, double y2, double y3, double y4) {
  double t1 = Zij(eta, y1, y4) + Zij(eta, y2, y1);
  double t2 = Zij(eta, y3, y2) - Zij(eta, y3, y4);
  double t3 = 0.5 * Chi(eta, y1, y2, y3, y4) * Chi(eta, y2, y3, y1, y4);
  double F = t1+t2+t3;
  if(IsBad(F)){
    msg_Error() <<METHOD << "\n "
                << "\n eta    = " << eta
                << "\n Y1    = " << y1
                << "\n Y2    = " << y2
                << "\n Y3    = " << y3
                << "\n T$    = " << y4
                << "\n T1    = " << t1
                << "\n T2    = " << t2
                << "\n T3    = " << t3<< std::endl;
  }
  return F;
}



double YFS_Form_Factor::Chi(double eta, double yi, double yj, double yk, double yl) {
  double nom = (eta - yi) * (eta - yj);
  double den = (eta - yk) * (eta - yl);
  return log(Abs(nom / den));
}

double YFS_Form_Factor::Zij(double eta, double yi, double yj) {
  double x = (yj - yi) / (eta - yi);
  double t1 = 2.*DiLog(x);
  double nom = eta - yi;
  double den = eta - yj;
  double t2 = 0.5 * sqr(log(Abs(nom / den)));
  return t1 + t2;

}


double YFS_Form_Factor::A4(double p1p2, double En1, double En2,
                           double mass1, double mass2) {
  double p1s = En1 * En1 - mass1 * mass1;
  double p2s = En2 * En2 - mass2 * mass2;
  if (p1s < p2s ) {
    double tempM1 = mass1;
    double tempE1 = En1;
    En1 = En2;
    mass1 = mass2;
    En2 = tempE1;
    mass2 = tempM1;
  }
  double Ep = En1 + En2;
  double Em = En1 - En2;
  double sm = mass1 + mass2;
  double dm = mass1 - mass2;
  double Q2 = 2.*p1p2 - mass1 * mass1 - mass2 * mass2;
  double xl = sqrt((Q2 + sm * sm) * (Q2 + dm * dm));
  double xq = sqrt(Q2 + Em * Em);
  double qp = xq + Em;
  double qm = xq - Em;

  double eta0 = sqrt(En2 * En2 - mass2 * mass2);
  if (p1p2 > En1 * En2) eta0 = -eta0;
  double eta1 = sqrt(En1 * En1 - mass1 * mass1) + xq;
  double y1  = 0.5 * ((xq - Ep) + (sm * dm + xl) / qp);
  double y2  = y1 - xl / qp;
  double y3  = 0.5 * ((xq + Ep) + (sm * dm + xl) / qm );
  double y4  = y3 - xl / qm;
  double Eln = 0;
  if (Abs(Em) > 1e-10 ) Eln = log(Abs(qm / qp)) * (Chi(eta1, y1, y4, y2, y3)
                                - Chi(eta0, y1, y4, y2, y3));
  // if(IsBad(Eln)) return 0;
  double V0 = YijEta(eta0, y1, y2, y3, y4);
  double V1 = YijEta(eta1, y1, y2, y3, y4);
  double A = 1. / xl * (Eln + V1 - V0);
  if (IsBad(A)) {
    msg_Error() << METHOD << "\n" << "Xl = " << xl << "\n"
                << "Eln= " << Eln << "\n"
                << "Q2= " << Q2 << "\n"
                << "V1= " << V1 << "\n"
                << "V0= " << V0 << "\n"
                << "qm = " << qm << "\n"
                << "qp = " << qp << "\n"
                << "Ep = " << Ep << "\n"
                << "Em = " << Em << "\n"
                << "En2*En2-mass2*mass2 = " << En2*En2 - mass2*mass2 << "\n"
                << "En1*En1-mass1*mass1 = " << En1*En1 - mass1*mass1 << "\n"
                << "eta1 = " << eta1 << "\n"
                << "eta0 = " << eta0 << "\n";
  }
  return A;
}

double YFS_Form_Factor::A4light(double p1p2, double En1, double En2,
                                double mass1, double mass2) {
  // courtesy of Jadach YFSWW, does not seem to be documented
  double Ep  = En1 + En2;
  double Em  = En1 - En2;
  double Ecm = 4 * En1 * En2;
  double xm12 = mass1 * mass2;
  double Q2 = 2 * p1p2;
  double Del = sqrt(Q2 + Em * Em);
  double DmEm = Del - Em;
  double DpEm = Del + Em;
  double DmEp = Del - Ep;
  double DpEp = Del + Ep;
  double p1p2A4 = - log(Q2 / xm12) * log(Ecm / xm12) + 0.5 * sqr(log(Q2 / xm12))
                  - 0.5 * sqr(log(mass1 / mass2)) + log(En1 / En2) * log(mass1 / mass2)
                  - 0.25 * sqr(log(sqr(DpEm) / Ecm)) - 0.25 * sqr(log(sqr(DmEm) / Ecm))
                  - 0.5 * sqr(log(En1 / En2)) + sqr(M_PI) / 6
                  - DiLog(DpEp / DpEm) - DiLog(DpEp / DmEm)
                  - DiLog(DmEp / DpEm) - DiLog(DmEp / DmEm);
  return p1p2A4 / p1p2;
}


double YFS_Form_Factor::A4_eq(double E, double M) {
  double bet = sqrt(1 - sqr(M / E));
  double b1ln = 2 * log((1. + bet) * E / M);
  return 1 / sqr(M) * b1ln / bet;
}

double YFS_Form_Factor::BVV_full(const ATOOLS::Vec4D p1, const ATOOLS::Vec4D p2, double MasPhot, double Kmax, int mode) {
  double t1, t2, t3;
  double alpi = m_alpha / M_PI;
  double Mas1 = p1.Mass();
  double Mas2 = p2.Mass();
  double m12 = Mas1 * Mas2;
  double E1 = p1.E();
  double E2 = p2.E();
  double p1p2 = p1 * p2;
  double rho = sqrt(1. - sqr(m12 / p1p2));
  double s = (p1 + p2).Abs2();
  double zeta1 = 2 * p1p2 * rho / (sqr(Mas1) + p1p2 * (1. + rho));
  double zeta2 = 2 * p1p2 * rho / (sqr(Mas2) + p1p2 * (1. + rho));
  double beta1 = sqrt(1. - sqr(Mas1 / E1));
  double beta2 = sqrt(1. - sqr(Mas2 / E2));
  double betat = 0.382;
  double beta  = sqrt(1. - 2 * (Mas1 + Mas2) / s + sqr((Mas1 - Mas2) / s));
  // t1 = (1./rho*A(p1p2,Mas1,Mas2)-1.)*2.*log(2.*Kmax/MasPhot);
  if (mode == 0 || mode == 3) {
    t1 = (log(p1p2 * (1. + rho) / m12) / rho - 1) * log(pow(MasPhot, 2) / m12);
    // t1 = (log(sqr(MasPhot)/sqr(250)));
    t2 = p1p2 * rho / s * log(p1p2 * (1. + rho) / m12) + (Mas1 * Mas1 - Mas2 * Mas2) / (2.*s) * log(Mas1 / Mas2) - 1;

    t3 =  -0.5 * log(p1p2 * (1. + rho) / sqr(Mas1)) * log(p1p2 * (1. + rho) / sqr(Mas2))
          - 0.5 * sqr(log((sqr(Mas1) + p1p2 * (1. + rho)) / (sqr(Mas2) + p1p2 * (1. + rho))));
    t3 -= DiLog(zeta1) + DiLog(zeta2);
    t3 += sqr(M_PI);
    t3 /= rho;
  }
  else {
    // deal with interpolation to coulomb if WW
    t1 = (log(p1p2 * (1. + rho) / m12) / rho - 1.) * log(sqr(MasPhot) / m12);
    t2 = p1p2 * rho / s * log(p1p2 * (1. + rho) / m12) + (Mas1 * Mas1 - Mas2 * Mas2) / (2.*s) * log(Mas1 / Mas2);
    t3 = sqr(M_PI) - 0.5 * log(p1p2 * (1. + rho) / sqr(Mas1)) * log(p1p2 * (1. + rho) / sqr(Mas2))
         - 0.5 * sqr(log((sqr(Mas1) + p1p2 * (1. + rho)) / (sqr(Mas2) + p1p2 * (1. + rho))));
    t3 += sqr(M_PI);

    t3 /= rho;
  }
  double virt = m_alpi * (t1 + t2 + t3);
  if (mode == 3) return m_alpi * (t1 + t2 + t3);
  if (mode==4) return m_alpi*t1;
  double real = BVR_full(p1p2, E1, E2, Mas1, Mas2, Kmax, MasPhot, mode);
  if (IsBad(real) || IsBad(virt)) {
    msg_Error() << METHOD << "\n"
                << "p1 = " << p1 << "\n"
                << "p2 = " << p2 << "\n"
                << "p1.Mass = " << p1.Mass() << "\n"
                << "p2.Mass = " << p2.Mass() << "\n"
                << "t1 = " << t1 << "\n"
                << "t2 = " << t2 << "\n"
                << "t3 = " << t3 << "\n"
                << "beta1 = " << beta1 << "\n"
                << "beta2 = " << beta2 << "\n"
                << "zeta1 = " << zeta1 << "\n"
                << "zeta2 = " << zeta2 << "\n"
                << "real = " << real << "\n"
                << "virt = " << virt << "\n"
                << "Mass Photon = " << m_photonMass << "\n";
  }
  return exp(real + virt);
}


DivArrC YFS_Form_Factor::BVV_full_eps(const ATOOLS::Vec4D p1, const ATOOLS::Vec4D p2, double MasPhot, double Kmax, int mode) {
  // for dim-reg
  // DivArrc {UV, IR, IR^2, finite, eps, eps^2, 0}
  double muf = 91.2*91.2;
  double mur = 91.2*91.2;
  double t2, t3;
  DivArrC t1;
  double alpi = m_alpha / M_PI;
  DivArrC massph(1,0,0,log(4*M_PI*mur)-GAMMA_E,1,0);
  double Mas1 = p1.Mass();
  double Mas2 = p2.Mass();
  double m12 = Mas1 * Mas2;
  double E1 = p1.E();
  double E2 = p2.E();
  double p1p2 = p1 * p2;
  double rho = sqrt(1. - sqr(m12 / p1p2));
  double s = (p1 + p2).Abs2();
  double zeta1 = 2 * p1p2 * rho / (sqr(Mas1) + p1p2 * (1. + rho));
  double zeta2 = 2 * p1p2 * rho / (sqr(Mas2) + p1p2 * (1. + rho));
  double beta1 = sqrt(1. - sqr(Mas1 / E1));
  double beta2 = sqrt(1. - sqr(Mas2 / E2));
  double betat = 0.382;
  double beta  = sqrt(1. - 2 * (Mas1 + Mas2) / s + sqr((Mas1 - Mas2) / s));
  // t1 = (1./rho*A(p1p2,Mas1,Mas2)-1.)*2.*log(2.*Kmax/MasPhot);
  t1 = (log(p1p2 * (1. + rho) / m12) / rho - 1) * (massph-log(m12));
  // t1 = (log(sqr(MasPhot)/sqr(250)));
  t2 = p1p2 * rho / s * log(p1p2 * (1. + rho) / m12) + (Mas1 * Mas1 - Mas2 * Mas2) / (2.*s) * log(Mas1 / Mas2) - 1;

  t3 =  -0.5 * log(p1p2 * (1. + rho) / sqr(Mas1)) * log(p1p2 * (1. + rho) / sqr(Mas2))
        - 0.5 * sqr(log((sqr(Mas1) + p1p2 * (1. + rho)) / (sqr(Mas2) + p1p2 * (1. + rho))));
  t3 -= DiLog(zeta1) + DiLog(zeta2);
  t3 += sqr(M_PI);
  t3 /= rho;
  PRINT_VAR(t1.GetEpsilon());
  return m_alpi * (t1 + t2 + t3);
}



double YFS_Form_Factor::WW_t(double t, double m, double M, double k) {
  // t and u virtual form factor
  t = Abs(t);
  double mm = m * M;
  double m2 = m * m;
  double M2 = M * M;
  // double bigL = 0.5*(log(sqrt(t)/m2)+log(sqrt(t)/M2));
  double bigL = log(t / mm);
  double zeta = 1 + M2 / t;
  double t1 = (bigL + log(zeta) - 1) * 2.*log(m_photonMass / m);
  double t2 = 0.5 * zeta * (bigL + log(zeta));
  // double t3 = -0.5*log(t/m2)*log(t/M2)*(bigL+log(zeta)+(zeta-3)/2.)-0.5*log(t/m2)*log(t/M2);
  double t3 = -0.5 * log(t / m2) * log(t / M2) - log(M / m) * (bigL + log(zeta) + 0.5 * (zeta - 3.));
  double t4 = -log(zeta) * (bigL + 0.5 * log(zeta)) + DiLog(1. / zeta) - 1;

  double rel = m_alpi * (t1 + t2 + t3 + t4);
  if (IsBad(rel)) {
    msg_Out() << METHOD << "\n"
              << "(p1-q1)**2 = " << t << "\n"
              << "t1 = " << t1 << "\n"
              << "t2 = " << t2 << "\n"
              << "t3 = " << t3 << "\n"
              << "t4 = " << t4 << "\n"
              << "res = " << rel << "\n"
              << "zeta = " << zeta << "\n"
              << "m = " << m << "\n"
              << "M = " << M << "\n"
              << "alpi = " << m_alpi << "\n";
  }
  return rel;
}

double YFS_Form_Factor::WW_s(Vec4D p1, Vec4D p2) {
  double betat = 0.382;
  double alpi = m_alpha / M_PI;
  double E1 = p1.E();
  double E2 = p2.E();
  double am1s = p1.Abs2();
  double am2s = p2.Abs2();
  double am1  = sqrt(am1s);
  double am2  = sqrt(am2s);
  double am12 = am1 * am2;
  double p1p2 = p1 * p2;
  double s    = 2 * p1p2 + am1s + am2s;
  double beta = sqrt( 1 - 2 * (am1s + am2s) / s + sqr((am1s - am2s) / s) );
  double rho   = sqrt( (1 + am12 / p1p2) * (1 - am12 / p1p2) );
  double pro   = 2 * p1p2 * rho;
  double opr   = p1p2 * (1 + rho);
  double oprm1 = opr + am1s;
  double oprm2 = opr + am2s;
  double Bigl  = log(opr / am12);
  double t1 = (Bigl / rho - 1.) * log(m_photonMass * m_photonMass / am12) + p1p2 * rho / s * log(p1p2 * (1 + rho) / am12) + (am1s - am2s) / (2 * s) * log(am1 / am2);
  double t2 = sqr(M_PI) / 2. - 0.5 * log(p1p2 * (1. + rho) / am1s) * log(p1p2 * (1. + rho) / am2s) - 0.5 * sqr(log((am1s + p1p2 * (1 + rho)) / (am2s + p1p2 * (1 + rho))));
  t2 -= DiLog(2 * p1p2 * rho / (am1s + p1p2 * (1 + rho))) + DiLog(2 * p1p2 * rho / (am2s + p1p2 * (1 + rho)));
  t2 /= rho;
  double t3 = -1;
// ! Interpolation - to match with Coulomb correction
  if (beta > betat && m_useCoulomb) {
    t2 = t2 + M_PI * M_PI / rho;
  }
  else {
    t2 = t2 + M_PI * M_PI * beta / 2.;
  }
  // PRINT_VAR(ReB2pi);
  // if(exp(m_alpi*(t1+t2+t3)) > 50) return 1;
  return exp(m_alpi * (t1 + t2 + t3));
}

double YFS_Form_Factor::BVV_WW(const ATOOLS::Vec4D_Vector born, const ATOOLS::Vec4D_Vector k, const ATOOLS::Vec4D p1, const ATOOLS::Vec4D p2, double MasPhot, double Kmax) {
  // p1 p2 will be the W+W- constructed four momenta
  // born is the born level momenta ! does not have W+W_
  // This is annoying but we have to make sure for now that W- = p2 and W+ = p1
  // todo add flavour map to get this correct?
  // test point from yfsww @ 160
  // p1 =    0.0000000000000000        0.0000000000000000        79.999999998367997        80.000000000000000
  // p2 =    0.0000000000000000        0.0000000000000000       -79.999999998367997        80.000000000000000
  // q1 =   -23.154167635022180       -4.9208126148046487        9.5182228410243219        83.312836026651667
  // q2 =    23.154207795719319        4.9208053213199383       -4.6576849851207687        71.826626117230816
  // double t1 =  -5517.0359043786993;
  // double t2 =  -6169.9920545855257;
  // double u1 =  -7660.4512497937658;
  // double u2 =  -8562.8672134443495;
  // // Virtual correctiion
  // double  Vt1 =  -1.4040641465047601;
  // double  Vt2 =  -1.4135670807177423;
  // double  Vu1 =  -1.4327761837313551;
  // double Vu2 =  -1.4380576392504969;
  m_wm = p2;
  m_wp = p1;
  m_ww_u = m_ww_t = 0;
  // m_photonMass = MasPhot;
  m_beam1 = born[1];//e-
  m_beam2 = born[0];//e+
  // m_beam1 = {80,0,0,79.999999998367997};
  // m_beam2 = {80,0,0,-79.999999998367997};
  // m_wm = {83.312836026651667, -23.154167635022180,-4.9208126148046487,9.5182228410243219};
  // m_wp = {71.826626117230816, 23.154207795719319, 4.9208053213199383, -4.6576849851207687};
  double s = (m_beam1 + m_beam2).Abs2();
  // m_ww_s = BVV_full(m_wm, m_wp, MasPhot, Kmax, 0);
  m_ww_s = exp(BVR_full(p1 * p2, p1.E(), p2.E(), p1.Mass(), p2.Mass(), Kmax, MasPhot, 0));
  m_ww_s *= WW_s(p1, p2);
  m_t1 = (m_beam1 - m_wm).Abs2();
  m_t2 = (m_beam2 - m_wp).Abs2();

  m_u1 = (m_beam1 - m_wp).Abs2();
  m_u2 = (m_beam2 - m_wm).Abs2();

  // PRINT_VAR(m_t1-t1);
  // PRINT_VAR(m_t2-t2);
  // PRINT_VAR(m_u1-u1);
  // PRINT_VAR(m_u2-u2);
  // PRINT_VAR(m_wm.Mass());
  double relt1 = WW_t(m_t1, m_beam1.Mass(), m_wm.Mass(), 1. );
  double relt2 = WW_t(m_t2, m_beam2.Mass(), m_wp.Mass(), 1. );
  double relu1 = WW_t(m_u1, m_beam1.Mass(), m_wp.Mass(), 1. );
  double relu2 = WW_t(m_u2, m_beam2.Mass(), m_wm.Mass(), 1. );
  double rel = relt1 + relt2 - relu1 - relu2;

  m_ww_t += BVR_full(m_beam1, m_wm, Kmax, MasPhot, 1);
  m_ww_t += BVR_full(m_beam2, m_wp, Kmax, MasPhot, 1);
  m_ww_u += BVR_full(m_beam1, m_wp, Kmax, MasPhot, 1);
  m_ww_u += BVR_full(m_beam2, m_wm, Kmax, MasPhot, 1);

  // PRINT_VAR(m_ww_t-Vt1);
  // PRINT_VAR(relt2-Vt2);
  // PRINT_VAR(relu1-Vu1);
  // PRINT_VAR(relu2-Vu2);
  // exit(1);

  double weik = 1;
  double eikii, eikff, eikif, eikfi, eikffff, eikaC, eikbD, eikaD, eikbC;

  double p1p2 = m_beam1 * m_beam2;
  double p1q1 = m_beam1 * m_wm;
  double p2q1 = m_beam2 * m_wm;
  double p1q2 = m_beam1 * m_wp;
  double p2q2 = m_beam2 * m_wp;
  double q1q2 = m_wm * m_wp;
  // calculate eikonals for ee->ww
  if (k.size() != 0) {
    weik = 1;
    for (auto kk : k)
    {
      double p1k = m_beam1 * kk;
      double p2k = m_beam2 * kk;
      double q1k = m_wm * kk;
      double q2k = m_wp * kk;

      // eikif = 2*(p1q1/(p1k*q1k) -p1q2/(p1k*q2k)
      // - p2q1/(p2k*q1k) +p2q2/(p2k*q2k));
      // eikff = 2*q1q2/(q1k*q2k) -sqr(m_wm.Mass()/q1k) -sqr(m_wp.Mass()/q2k);
      // eikii = 2*p1p2/(p1k*p2k) -sqr(m_beam1.Mass()/p1k) -sqr(m_beam2.Mass()/p2k);
      // PRINT_VAR(eikii/(m_beam1/p1k-m_beam2/p2k).Abs2());
      eikii = -(m_beam1 / p1k - m_beam2 / p2k).Abs2();
      eikff = -(m_wm / q1k - m_wp / q2k).Abs2();
      eikif = 2 * (m_beam1 / p1k - m_beam2 / p2k) * (m_wm / q1k - m_wp / q2k);
      weik *= 1 + (eikif) / (eikii + eikff);
    }
  }
  return m_ww_s * exp(rel + m_ww_t - m_ww_u) * weik;
}



double YFS_Form_Factor::BVirtT(const Vec4D &p1, const Vec4D &p2){
  double m1 = p1.Mass();
  double m2 = p2.Mass();
  double p1p2 = p1*p2;
  double t  = (p1-p2).Abs2();
  double ta = abs(t);
  double zeta = 1+ m2*m2/ta;
  double TBvirt, Bv;
  TBvirt = m_alpi*(
       (log(2.0*p1p2/(m1*m2)) -1.0)*log(m_photonMass*m_photonMass/(m1*m2))
       +0.5*zeta*log(ta*zeta/(m1*m2))
       -0.5*log(ta/m1/m1)*log(ta/m2/m2)
       +DiLog(1/zeta) -1.0
       +0.5*(zeta -1.0)*log(m1/m2)
       -log(zeta)*(log(ta/(m1*m2)) +0.5*log(zeta))
       );
  Bv = BVR_full(p1,p2, sqrt(m_s)/2, m_photonMass, 0);
  return exp(TBvirt+Bv);
}

double YFS_Form_Factor::R1(const Vec4D &p1, const Vec4D &p2){
  return exp(R2(p1,p2) + BVR_full(p1,p2,p1.E()*p2.E(), m_photonMass, 0));
}


double YFS_Form_Factor::R2(const Vec4D &p1, const Vec4D &p2){
  double beta1 = (Vec3D(p1).Abs() / p1.E());
  double beta2 = (Vec3D(p2).Abs() / p2.E());
  double logarg =  (1+beta1)*(1+beta2);
  logarg /= (1-beta2)*(1-beta1);

  double biglog = 0.5*(1+beta1*beta2)/(beta1+beta2);
  biglog *= (log(logarg)-2);

  double logp = 0.5*(1+beta1*beta2)/(beta1+beta2);
  logp *= logp;

  double t1 = biglog*2*log(m_photonMass/(p1.E()*p2.E()))+0.25*logp;
  // double t1 = biglog+0.25*logp;
  t1+= -0.5*sqr(log(p1.E()/p2.E()));

  double del = p1.E()-p2.E(); 
  double Delta = sqrt(2*p1*p2-del*del);
  double omega = p1.E()+p2.E();

  double t2 = -0.25*sqr(log(sqr(del+Delta)/(4*p1.E()*p2.E())));
  t2 += -0.25*sqr(log(sqr(del-Delta)/(4*p1.E()*p2.E())));

  t2 += -DiLog((Delta+omega)/(Delta+del)) - DiLog((Delta+omega)/(Delta-del));
  t2 += -DiLog((Delta-omega)/(Delta+del)) - DiLog((Delta-omega)/(Delta-del));
  t2 += M_PI*M_PI/3. - 1.;
  // if(IsNan(t2)) t2 = 0;
  return m_alpi*(t1+t2);
}


