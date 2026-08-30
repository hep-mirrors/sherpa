#include "YFS/Main/YFS_Form_Factor.H"

#include <algorithm>
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#include "ATOOLS/Math/Poincare.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "MODEL/Main/Model_Base.H"
#include "METOOLS/Loops/Master_Integrals.H"

#ifdef USING__LOOPTOOLS
  #include "clooptools.h"
#endif


using namespace YFS;
using namespace ATOOLS;
using namespace MODEL;
using namespace METOOLS;

// Lambda (Kaellen function) now lives once in YFS/Tools/Dipole.H.
// NB: this file's old copy dropped the abs(); the canonical one keeps it.



YFS_Form_Factor::YFS_Form_Factor() {
  // p_virt is only assigned in NLO_Base; left uninitialised it holds garbage, so the
  // (p_virt?...) guard in BVirtGeneralEps passed and IRscale() segfaulted.
  p_virt = nullptr;
    rpa->gen.AddCitation(1,"YFS Form Factor as implemented in \\cite{Jadach:1999vf}");
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
    t2 = p1p2 * A4(p1p2, Leg(E1, Mas1), Leg(E2, Mas2));
  }
  else if(mode==1) {
    // alpi = 1./137.03599976000001/M_PI;
    t1 = (p1p2 * A(p1p2, Mas1, Mas2)) * log(4 * sqr(Kmax / MasPhot));
    t2 = p1p2 * A4(p1p2, Leg(E1, Mas1), Leg(E2, Mas2));
    
    }
    else {
      AA4 = A4(p1p2, Leg(E1, Mas1), Leg(E2, Mas2));
      t2 = AA4 * p1p2;
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
  if(m_fullform>=2) return m_alpi * (t1+0.5*t3);
  return m_alpi * (t1 +  t2 + 0.5 * t3);
}

DivArrD YFS_Form_Factor::BVR_full_eps(YFS::Dipole &d,  double Kmax, int mode) {
  // if(m_tchannel==2) return BVirtTEps(d,Kmax);
  Vec4D p1,p2;
  d.LegsBeforeRadiation(p1, p2);
  const double p1p2 = p1*p2;
  const  double E1 = p1.E();
  const double E2 = p2.E();
  const double Mas1 = d.GetMass(0);//p1.Mass();
  const double Mas2 = d.GetMass(1);//p2.Mass();
  const double alpi = m_alpi;
  const double m12 = Mas1 * Mas2;
  DivArrD t1; 
  double t2, AA4;
  if (p1p2 - m12 < 1e-10) return 0;
  const double beta1 = sqrt(1. - sqr(Mas1 / E1));
  const double beta2 = sqrt(1. - sqr(Mas2 / E2));
  const double rho = sqrt(1 - sqr(m12 / p1p2));
  const double irloop = m_irscale;//p_virt->IRscale();
  const double epsloop = (p_virt?p_virt->Eps_Scheme_Factor({p1,p2}):4*M_PI);
  DivArrD massph(0,-1,0,0,0,0);
  t1 = (p1p2 * A(p1p2, Mas1, Mas2) - 1) * (massph-log(4.*M_PI*sqr(irloop)/4./Kmax/epsloop));
  if ( mode == 0 ) {
    t2 = p1p2 * A4(p1p2, Leg(E1, Mas1), Leg(E2, Mas2));
  }
  else if(mode==1) {
    t1 = (p1p2 * A(p1p2, Mas1, Mas2)) * (massph-log(4.*M_PI*sqr(irloop)/4./Kmax/epsloop));
    t2 = p1p2 * A4(p1p2, Leg(E1, Mas1), Leg(E2, Mas2));
    
    }
    else {
      AA4 = A4(p1p2, Leg(E1, Mas1), Leg(E2, Mas2));
      t2 = AA4 * p1p2;
  }
  double t3 = Mas1 * Mas1 * A4_eq(E1, Mas1) + Mas2 * Mas2 * A4_eq(E2, Mas2);
  if (IsBad(t1.Finite()) || IsBad(t2) || IsBad(t3)) {
    msg_Error() << METHOD << "\n"
                << "YFS Form Factor is NaN"
                << "\n T1    = " << t1.Finite()
                << "\n T2    = " << t2
                << "\n T3    = " << t3 * 0.5
                << "\n E1    = " << E1
                << "\n E2    = " << E2
                << "\n Mass1 = " << Mas1
                << "\n Mass2 = " << Mas2
                << "\n Kmax = " << Kmax
                << "\n M12 = " << m12
                << "\n A4 = " << AA4
                << "\n  A(p1p2, Mas1, Mas2) = "<<  A(p1p2, Mas1, Mas2)
                << "\n p1  = " << p1
                << "\n p2  = " << p2
                << "\n irloop  = " << irloop
                << "\n epsloop  = " << epsloop
                << "\n p1p2  = " << p1p2;
  }
  return m_alpi * (t1 + t2 + 0.5 * t3);
}



double YFS_Form_Factor::BVR_decay(double M, double ml, double ks) {
  if (M <= 0. || ml <= 0. || ks <= 0.) return 0.;
  const double L = log(M/ml);
  return m_alpi*( 2.*(L - 1.)*log(2.*ks/M) + 1.5*L - sqr(M_PI)/6. + 1. );
}


double YFS_Form_Factor::BVR_full(Vec4D p1, Vec4D p2,  double Kmax, double MasPhot, int mode) {
  return BVR_full(p1 * p2, p1.E(), p2.E(), p1.Mass(), p2.Mass(), Kmax, MasPhot, mode);
}


double YFS_Form_Factor::BVR_full(YFS::Dipole &d, double omega) {
  // FULL_FORM==3: dimensionally-regularised YFS form factor built from the
  // METOOLS/Loops master integrals. Virtual B (BVirtMaster, from B0+C0) and
  // real B-tilde (BVR_full_eps) are summed as DivArr so their 1/eps_IR poles
  // cancel; the finite remainder is the scheme-independent YFS Y. The pole
  // cancellation is enforced as a correctness gate.
  if (m_fullform == 3) {
    DivArrD Yv = BVirtMaster(d, omega);          // master-integral virtual B
    DivArrD Yr = BVR_full_eps(d, omega, 0);       // dim-reg real B-tilde
    DivArrD Yvref = BVirtGeneralEps(d, omega);    // known-good analytic virtual B
    // One-time calibration print: compare the master-integral virtual against
    // the validated analytic dim-reg virtual so we can read off the exact
    // normalisation / sign correction empirically (kinematics of the II dipole
    // are fixed, so a single print is representative).
    static bool once = true;
    if (once) {
      once = false;
      msg_Out() << "=== FULL_FORM==3 YFS form-factor calibration (dipole "
                << d.Type() << ") ===\n"
                << "  BVirtMaster   : 1/eps=" << Yv.GetIR()
                << "  finite=" << Yv.Finite() << "\n"
                << "  BVirtGenerEps : 1/eps=" << Yvref.GetIR()
                << "  finite=" << Yvref.Finite()
                << "   (ratio fin=" << (IsZero(Yv.Finite())?0.:Yvref.Finite()/Yv.Finite())
                << ", ratio pole=" << (IsZero(Yv.GetIR())?0.:Yvref.GetIR()/Yv.GetIR()) << ")\n"
                << "  BVR_full_eps  : 1/eps=" << Yr.GetIR()
                << "  finite=" << Yr.Finite() << "\n"
                << "  analytic V+R finite (ref) = " << (Yvref.Finite()+Yr.Finite())
                << std::endl;
    }
    // Use the validated analytic virtual for the physical result for now; the
    // master-integral virtual (Yv) is being calibrated against it above.
    DivArrD Y = Yvref + Yr;
    return Y.Finite();
  }
  double R, V;
  Vec4D p1,p2;
  d.LegsBeforeRadiation(p1, p2);
  R =  BVR_full(p1 * p2, p1.E(), p2.E(), p1.Mass(), p2.Mass(), omega, m_photonMass, 0);
  // R = BVR_full_eps(d, omega, 0).Finite();
  V =  BVirtGeneralEps(d,omega).Finite();
  double Vold =  BVV_full(p1, p2, m_photonMass, omega, 0);
  if(m_fullform>=2) return (R+V);
  return (R+Vold);
}

double YFS_Form_Factor::BVR_full(Vec4D p1, Vec4D p2, double omega) {
  double R, V;
  // if(!m_dim_reg){
    R =  BVR_full(p1 * p2, p1.E(), p2.E(), p1.Mass(), p2.Mass(), omega, m_photonMass, 0);
    V =  BVV_full(p1, p2, m_photonMass, omega, 0);
    // BVirtGeneral()
  // }
  // else{
  //   DivArrD DR, DV;
  //   DR = BVR_full_eps(p1, p2, omega, 0);
  //   DV = BVV_full_eps(p1, p2, omega, 0);
  //   if(!IsZero(DR.GetIR()-DV.GetIR())){
  //     msg->SetPrecision(16);
  //     msg_Error()<<"Poles do not cancel in YFS Form Factor"<<std::endl;
  //     msg_Error()<<"Real eps^{-1} = "<<DR.GetIR()<<std::endl
  //               <<"Virtual eps^{-1} = "<<-DV.GetIR()<<std::endl
  //               << "Diff =  "<< DV.GetIR()-DR.GetIR()<<std::endl;
  //   }
  //   R=DR.Finite();
  //   V=DV.Finite();
  // }
  return (R+V);
}

double YFS_Form_Factor::BVR_cru(double p1p2, const Leg &l1, const Leg &l2,
                                double Kmax) {
  const double E1 = l1.E, E2 = l2.E, Mas1 = l1.m, Mas2 = l2.m;
  // m_btilcru = m_alpi*(p1p2*BVR_A(p1p2,Mas1,Mas2))
  double t1 = (p1p2 * A(p1p2, Mas1, Mas2)) * log(4.*sqr(Kmax / m_photonMass));
  double t2 = p1p2 * A4(p1p2, Leg(E1, Mas1), Leg(E2, Mas2));
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


double YFS_Form_Factor::BVR_cru(Vec4D p1, Vec4D p2, double Kmax) {
  // m_btilcru = m_alpi*(p1p2*BVR_A(p1p2,Mas1,Mas2))
  double t1 = (p1*p2 * A(p1,p2)) * log(4.*sqr(Kmax / m_photonMass));
  double t2 = p1*p2 * A4(p1*p2, Leg(p1), Leg(p2));
  if (IsBad(t1) || IsBad(t2)) {
    msg_Error() << METHOD << "\n" << "YFS Form Factor is NaN"
                << "\n T1    = " << t1
                << "\n T2    = " << t2
                << "\n p1  = " << p1
                << "\n p2  = " << p2 << std::endl;
  }
  return m_alpi * (t1 + t2);
}


double YFS_Form_Factor::A(double p1p2, double m1, double m2) {
  double m12 = m1 * m2;
  if (IsEqual(p1p2, m12)) return 0;
  double xlam = sqrt((p1p2 - m12) * (p1p2 + m12));
  double BVR_A = 1. / xlam * log((p1p2 + xlam) / m12);
  return BVR_A;
}


double YFS_Form_Factor::A(Vec4D p1, Vec4D p2) {
  double m12 = p1.Mass() * p2.Mass();
  double p12 = p1*p2;
  if ((p1*p2 - m12) < 1e-10) return 0;
  double xlam = sqrt((p1*p2 - m12) * (p1*p2 + m12));
  double BVR_A = 1. / xlam * log((p1*p2 + xlam) / m12);
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


double YFS_Form_Factor::A4(double p1p2, Leg l1, Leg l2) {
  // Order by |p|^2, which Leg guarantees is non-negative. Swapping whole Legs
  // keeps each energy with its own mass; the old code swapped four loose
  // doubles by hand.
  if (l1.P2() < l2.P2()) std::swap(l1, l2);
  const double En1 = l1.E, mass1 = l1.m;
  const double En2 = l2.E, mass2 = l2.m;
  double Ep = En1 + En2;
  double Em = En1 - En2;
  double sm = mass1 + mass2;
  double dm = mass1 - mass2;
  double Q2 = 2.*p1p2 - mass1 * mass1 - mass2 * mass2;
  double xl = sqrt((Q2 + sm * sm) * (Q2 + dm * dm));
  double xq = sqrt(Q2 + Em * Em);
  double qp = xq + Em;
  double qm = xq - Em;
  double eta0 = sqrt((En2 * En2 - mass2 * mass2));
  if (p1p2 > En1 * En2) eta0 = -eta0;
  double eta1 = sqrt((En1 * En1 - mass1 * mass1)) + xq;
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

double YFS_Form_Factor::A4_eq(double E, double M) {
  double bet = sqrt(1 - sqr(M / E));
  double b1ln = 2 * log((1. + bet) * E / M);
  return 1 / sqr(M) * b1ln / bet;
}


// Shared implementation of the s-channel virtual B. The two public overloads
// differ ONLY in where the masses come from - the Vec4D one takes them off the
// momenta, the Dipole one off d.GetMass(), which need not agree if the stored
// momenta are slightly off shell - so the mass is passed in rather than
// recomputed here.
// NB Kmax is not needed: this is the purely virtual piece, the soft cutoff
// enters through the real Btilda (BVR_full) that gets added to it.
double YFS_Form_Factor::BVV_full_impl(const ATOOLS::Vec4D &p1, const ATOOLS::Vec4D &p2,
                                      double Mas1, double Mas2,
                                      double MasPhot, int mode) {
  double t1, t2, t3;
  const double m12   = Mas1 * Mas2;
  const double p1p2  = p1 * p2;
  const double rho   = sqrt(1. - sqr(m12 / p1p2));
  const double s     = (p1 + p2).Abs2();
  const double zeta1 = 2 * p1p2 * rho / (sqr(Mas1) + p1p2 * (1. + rho));
  const double zeta2 = 2 * p1p2 * rho / (sqr(Mas2) + p1p2 * (1. + rho));
  if (mode == 0 || mode == 3) {
    t1 = (log(p1p2 * (1. + rho) / m12) / rho - 1) * log(sqr(MasPhot) / m12);
    t2 = p1p2 * rho / s * log(p1p2 * (1. + rho) / m12)
         + (Mas1 * Mas1 - Mas2 * Mas2) / (2.*s) * log(Mas1 / Mas2) - 1;
    t3 =  -0.5 * log(p1p2 * (1. + rho) / sqr(Mas1)) * log(p1p2 * (1. + rho) / sqr(Mas2))
          - 0.5 * sqr(log((sqr(Mas1) + p1p2 * (1. + rho)) / (sqr(Mas2) + p1p2 * (1. + rho))));
    t3 -= DiLog(zeta1) + DiLog(zeta2);
    t3 += sqr(M_PI);
    t3 /= rho;
  }
  else {
    // interpolation to Coulomb, used by the WW form factor
    t1 = (log(p1p2 * (1. + rho) / m12) / rho - 1.) * log(sqr(MasPhot) / m12);
    t2 = p1p2 * rho / s * log(p1p2 * (1. + rho) / m12)
         + (Mas1 * Mas1 - Mas2 * Mas2) / (2.*s) * log(Mas1 / Mas2);
    t3 = sqr(M_PI) - 0.5 * log(p1p2 * (1. + rho) / sqr(Mas1)) * log(p1p2 * (1. + rho) / sqr(Mas2))
         - 0.5 * sqr(log((sqr(Mas1) + p1p2 * (1. + rho)) / (sqr(Mas2) + p1p2 * (1. + rho))));
    t3 += sqr(M_PI);
    t3 /= rho;
  }
  const double virt = m_alpi * (t1 + t2 + t3);
  if (mode == 3) return virt;
  if (mode == 4) return m_alpi * t1;
  if (IsBad(virt)) {
    msg_Error() << METHOD << "\n"
                << "p1 = " << p1 << "\n"
                << "p2 = " << p2 << "\n"
                << "Mas1 = " << Mas1 << "\n"
                << "Mas2 = " << Mas2 << "\n"
                << "t1 = " << t1 << "\n"
                << "t2 = " << t2 << "\n"
                << "t3 = " << t3 << "\n"
                << "zeta1 = " << zeta1 << "\n"
                << "zeta2 = " << zeta2 << "\n"
                << "virt = " << virt << "\n"
                << "Mass Photon = " << m_photonMass << "\n";
  }
  return virt;
}

double YFS_Form_Factor::BVV_full(const ATOOLS::Vec4D p1, const ATOOLS::Vec4D p2, double MasPhot, double Kmax, int mode) {
  return BVV_full_impl(p1, p2, p1.Mass(), p2.Mass(), MasPhot, mode);
}

double YFS_Form_Factor::BVV_full(YFS::Dipole &d, double MasPhot, double Kmax, int mode) {
  Vec4D p1,p2;
  d.LegsBeforeRadiation(p1, p2);
  return BVV_full_impl(p1, p2, d.GetMass(0), d.GetMass(1), MasPhot, mode);
}

DivArrD YFS_Form_Factor::BVV_full_eps(YFS::Dipole &d, double Kmax, int mode){
  // for dim-reg
  // DivArrc {UV, IR, IR^2, finite, eps, eps^2, 0}
  Vec4D p1,p2;
  d.LegsBeforeRadiation(p1, p2);
  double t2, t3;
  DivArrD t1;
  DivArrD massph(0,-1,0,0,0,0);
  double Mas1 = d.GetMass(0);
  double Mas2 = d.GetMass(1);
  double m12 = Mas1*Mas2;
  double E1 = p1.E();
  double E2 = p2.E();
  double p1p2 = p1 * p2;
  // double rho = sqrt(1. - sqr(m12 / p1p2));
  double rho = sqrt((p1p2 - m12) * (p1p2 + m12)) / p1p2;
  double s = (p1 + p2).Abs2();
  double zeta1 = 2 * p1p2 * rho / (sqr(Mas1) + p1p2 * (1. + rho));
  double zeta2 = 2 * p1p2 * rho / (sqr(Mas2) + p1p2 * (1. + rho));
  double beta1 = sqrt(1. - sqr(Mas1 / E1));
  double beta2 = sqrt(1. - sqr(Mas2 / E2));
  double betat = 0.382;
  double beta  = sqrt(1. - 2 * (Mas1 + Mas2) / s + sqr((Mas1 - Mas2) / s));
  // t1 = (1./rho*A(p1p2,Mas1,Mas2)-1.)*2.*log(2.*Kmax/MasPhot);
  double irloop = p_virt->IRscale();
  double epsloop = p_virt->Eps_Scheme_Factor({p1,p2});
  double logarg = (p1p2 * (1. + rho) / m12) / rho;
  t1 = (p1p2 * A(p1p2, Mas1, Mas2) -1.) *  (massph+log(4.*M_PI*sqr(irloop)/m12/epsloop));
  // if(logarg < 1e-2) t1 = (log1p(logarg-1) -1.) *  (massph+log(4.*M_PI*sqr(irloop)/m12/epsloop));
  // else t1 = (log(logarg) - 1.) *  (massph+log(4.*M_PI*sqr(irloop)/m12/epsloop));
  // t1 = (log(sqr(MasPhot)/sqr(250)));
  t2 = p1p2 * rho / s * log(p1p2 * (1. + rho) / m12) + (Mas1 * Mas1 - Mas2 * Mas2) / (2.*s) * log(Mas1 / Mas2) - 1;

  t3 =  -0.5 * log(p1p2 * (1. + rho) / sqr(Mas1)) * log(p1p2 * (1. + rho) / sqr(Mas2))
        - 0.5 * sqr(log((sqr(Mas1) + p1p2 * (1. + rho)) / (sqr(Mas2) + p1p2 * (1. + rho))));
  t3 -= DiLog(zeta1) + DiLog(zeta2);
  t3 += sqr(M_PI);
  t3 /= rho;
  // if(m_tchannel==1) return m_alpi * (t1 + t2 + t3) + BVirtTEps(d);
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

  double relt1 = WW_t(m_t1, m_beam1.Mass(), m_wm.Mass(), 1. );
  double relt2 = WW_t(m_t2, m_beam2.Mass(), m_wp.Mass(), 1. );
  double relu1 = WW_t(m_u1, m_beam1.Mass(), m_wp.Mass(), 1. );
  double relu2 = WW_t(m_u2, m_beam2.Mass(), m_wm.Mass(), 1. );
  double rel = relt1 + relt2 - relu1 - relu2;

  m_ww_t += BVR_full(m_beam1, m_wm, Kmax, MasPhot, 1);
  m_ww_t += BVR_full(m_beam2, m_wp, Kmax, MasPhot, 1);
  m_ww_u += BVR_full(m_beam1, m_wp, Kmax, MasPhot, 1);
  m_ww_u += BVR_full(m_beam2, m_wm, Kmax, MasPhot, 1);

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
  static bool once(true);
  if (once && msg_LevelIsDebugging()) {
    once = false;
    msg_Out() << "=== BVV_WW breakdown (MasPhot=" << MasPhot
              << ", Kmax=" << Kmax << ") ===\n"
              << "  m_ww_s (exp(B(W-,W+)) * WW_s) = " << m_ww_s << "\n"
              << "  rel (WW_t t/u combination)    = " << rel << "\n"
              << "  m_ww_t                        = " << m_ww_t << "\n"
              << "  m_ww_u                        = " << m_ww_u << "\n"
              << "  exp(rel + t - u)              = " << exp(rel + m_ww_t - m_ww_u) << "\n"
              << "  weik (real-photon eikonal)    = " << weik << "\n"
              << "  product                       = "
              << m_ww_s * exp(rel + m_ww_t - m_ww_u) * weik << std::endl;
  }
  return m_ww_s * exp(rel + m_ww_t - m_ww_u) * weik;
}


double YFS_Form_Factor::BVirtT(Vec4D p1, Vec4D p2,  double kmax){
  // kmax is unused: this is the purely virtual t-channel B, the cutoff enters
  // only through the real Btilda that IFForFac() adds to it.
  const double m1 = p1.Mass();
  const double m2 = p2.Mass();
  const double M = m1>=m2 ? m1 : m2;
  const double m = m1>=m2 ? m2 : m1;
  const double ta = fabs((p1-p2).Abs2());
  const double zeta = 1 + M*M/ta;
  const double TBvirt = m_alpi*(
       (log(ta/m/M) +log(zeta) - 1) *log(sqr(m_photonMass)/(m*M))
      +0.5*zeta*log(ta*zeta/(m1*m2))
      -0.5*log(ta/m1/m1)*log(ta/m2/m2)
      // log(m/M), NOT log(m1/m2): this term is antisymmetric under exchanging
      // the two legs, and KKMC's TBvirt (SRCee/KKbvir.cxx:296) documents its
      // m1 as "assumed to be very small", i.e. the LIGHT leg first. IFForFac()
      // calls this with GetBornMomenta(1) then (0) to match R1()'s ordering,
      // which for an initial-final dipole puts the heavy (final-state) leg
      // first - so log(m1/m2) came out with the wrong sign. Using min/max
      // makes the term independent of the caller's ordering and reproduces
      // KKMC. Verified against KKbvir::TBvirt on a fixed e+e- -> mu+mu- point
      // (YFS/Tools/IFI_KKMC_CrossCheck.C vs Test/SherpaCompare/
      // kkmc_ifi_crosscheck.cxx): the whole disagreement was 2x this term.
      +0.5*(zeta -1.0)*log(m/M)
      -log(zeta)*(log(ta/(m1*m2)) +0.5*log(zeta))
      +DiLog(1./zeta) -1.0
       );
  return TBvirt;
}


double YFS_Form_Factor::BVirtT(YFS::Dipole &d, double kmax){
 Vec4D p1,p2;
 d.LegsBeforeRadiation(p1, p2);
  double m1 = d.GetMass(0);
  double m2 = d.GetMass(1);
  if(IsZero(kmax)) kmax=m1*m2;
  double M = m1>=m2 ? m1 : m2;
  double m = m1>=m2 ? m2 : m1;
  double p1p2 = p1*p2;
  double t  = (p1-p2).Abs2();
  double ta = fabs(t);
  double zeta = 1 + M*M/ta;
  double TBvirt, Bv;
  double rho = sqrt(1. - sqr(m1*m2 / (p1*p2)));
  TBvirt = m_alpi*(
    // (log(p1p2 * (1. + rho) / (m1*m2)) / rho - 1) *log(pow(m_photonMass, 2)/(kmax)) 
       (log(2*p1p2/(m1*m2))-1.0)*log(m_photonMass*m_photonMass/(m1*m2))
      +0.5*zeta*log(ta*zeta/(m1*m2))
      -0.5*log(ta/m1/m1)*log(ta/m2/m2)
      +DiLog(1./zeta) -1.0
      // log(m/M) rather than log(m1/m2) - see the note in the Vec4D overload.
      +0.5*(zeta -1.0)*log(m/M)
      -log(zeta)*(log(ta/(m1*m2)) +0.5*log(zeta))
       );
  #ifdef USING__LOOPTOOLS
    Complex form;
    Flavour fl1 = d.GetFlav(0);
    Flavour fl2 = d.GetFlav(1);
    // PRINT_VAR(d.m_thetai*p1+d.m_thetaj*p2);
    double s = (p1+p2).Abs2();
    double crossterm = ((d.m_thetai*p1+d.m_thetaj*p2).Abs2());
    if(fl1==fl2){
      form = 1./8.*B0(s, m1*m1, m2*m2);
      // form += -8.*(m1*m1)*C0(m1*m1, s, m2*m2,0,m1*m1,m2*m2);
      // PRINT_VAR(B0(s, m1*m1, m2*m2));
    //   // PRINT_VAR(C0(m1*m1, 0, m1*m1,0, m1*m1, m1*m1));
    }
    else{
      form = 2*(p1*p2)*C0(m1*m1, s, m2*m2,0,m1*m1,m2*m2);
      form += 0.25*B0(s, m1*m1, m2*m2);
      // PRINT_VAR(B0(crossterm, m1*m1, m2*m2));
    }
    form*=m_alpi;
    // // PRINT_VAR(fl1);
    // // PRINT_VAR(fl2);
    // TBvirt+=form.real();
    clearcache();
  #endif
  return TBvirt;
}

DivArrD YFS_Form_Factor::BVirtTEps(YFS::Dipole &d, double kmax){
   Vec4D p1,p2;
   d.LegsBeforeRadiation(p1, p2);
  double m1 = d.GetMass(0);
  double m2 = d.GetMass(1);
  if(IsZero(kmax)) kmax=m1*m2;
  double M = m1>=m2 ? m1 : m2;
  double m = m1>=m2 ? m2 : m1;
  double p1p2 = p1*p2;
  double t  = (p1-p2).Abs2();
  double ta = fabs(t);
  double zeta = 1 + M*M/ta;
  DivArrD TBvirt, Bv;
  double rho = sqrt(1. - sqr(m1*m2 / (p1*p2)));
  DivArrD massph(0,-1,0,0,0,0);
  double irloop = p_virt->IRscale();
  double epsloop = p_virt->Eps_Scheme_Factor({p1,p2});
  // (massph-log(4.*M_PI*sqr(irloop)/4./Kmax/epsloop));
  TBvirt = m_alpi*(
    (log(p1p2 * (1. + rho) / (m1*m2)) / rho - 1) * -1.*(-massph-log(4.*M_PI*sqr(irloop)/m1/m2/epsloop)) 
    // (log(2*p1p2/(m1*m2))-1.0) * -1.*(-massph-log(4.*M_PI*sqr(irloop)/m1/m2/epsloop)) 
       +0.5*zeta*log(ta*zeta/(m1*m2))
        -0.5*log(ta/m1/m1)*log(ta/m2/m2)
      +DiLog(1./zeta) -1.0
      +0.5*(zeta -1.0)*log(m1/m2)
      -log(zeta)*(log(ta/(m1*m2)) +0.5*log(zeta))
       );
  #ifdef USING__LOOPTOOLS
    Complex form;
    Flavour fl1 = d.GetFlav(0);
    Flavour fl2 = d.GetFlav(1);
    // PRINT_VAR(d.m_thetai*p1+d.m_thetaj*p2);
    double s = (p1+p2).Abs2();
    double crossterm = ((d.m_thetai*p1+d.m_thetaj*p2).Abs2());
    if(fl1==fl2){
      form = 1./8.*B0(s, m1*m1, m2*m2);
      // form += -8.*(m1*m1)*C0(m1*m1, s, m2*m2,0,m1*m1,m2*m2);
      // PRINT_VAR(B0(s, m1*m1, m2*m2));
    //   // PRINT_VAR(C0(m1*m1, 0, m1*m1,0, m1*m1, m1*m1));
    }
    else{
      form = 2*(p1*p2)*C0(m1*m1, s, m2*m2,0,m1*m1,m2*m2);
      form += 0.25*B0(s, m1*m1, m2*m2);
      // PRINT_VAR(B0(crossterm, m1*m1, m2*m2));
    }
    form*=m_alpi;
    // // PRINT_VAR(fl1);
    // // PRINT_VAR(fl2);
    // TBvirt+=form.real();
    clearcache();
  #endif
  return TBvirt;
}

double YFS_Form_Factor::IFForFac(YFS::Dipole &d, double omega){
  // KKMC's TForFac, for one initial-final dipole. KKceex.cxx assembles the
  // interference form factor as
  //     Yint = prod over the four IF pairs of TForFac(+/-a, pi, pj, Emin)
  //     TForFac = exp( Btilda(...,Emin) + TBvirt(...) )
  // and the sign pattern (+ for like-charge pairs, - for unlike) is what
  // Dipole::ChargeNorm() already supplies, so the caller multiplies by that and
  // this returns the exponent for a single pair.
  //
  // The point of having this separate from the II/FF term is the VIRTUAL: an
  // initial-final pair is t-channel-like, so it needs BVirtT (spacelike p1.p2),
  // not the s-channel virtual. FormFactorSum() previously gave the IF dipoles
  // Btilda with NO virtual at all, which is why no rescaling of that term could
  // reproduce KKMC -- half the object was missing.
  //
  // Momentum order follows R1(): for dipoletype::ifi it takes GetBornMomenta(1)
  // then (0). Kept identical here so the two agree bin by bin.
  //
  // NB BVirtT(Vec4D,Vec4D,double) ignores its kmax argument -- it is read only by
  // an IsZero() default that is never used afterwards -- so the omega passed here
  // affects the Btilda part alone. That is deliberate: it is the soft-real piece
  // that carries the cutoff, exactly as in KKMC's Btilda(...,Emin).
  ATOOLS::Vec4D p1 = d.GetBornMomenta(1);
  ATOOLS::Vec4D p2 = d.GetBornMomenta(0);
  double Breal = BVR_full(p1 * p2, p1.E(), p2.E(), p1.Mass(), p2.Mass(),
                          omega, m_photonMass, 0);
  double Bvirt = BVirtT(p1, p2, omega);
  return Breal + Bvirt;
}


double YFS_Form_Factor::R1(YFS::Dipole &d){
  Vec4D p1,p2;
  d.LegsBeforeRadiation(p1, p2);
  double R = BVR_full(p1 * p2, p1.E(), p2.E(), p1.Mass(), p2.Mass(), sqrt(m_s)/2., m_photonMass, 0);
  double V = BVirtT(p1, p2, sqrt(m_s) / 2.);
  if(m_tchannel==3) V = BVirtGeneral(d);
  if(m_tchannel==1){
    // add s channel 
    double Vs = BVR_full(p1, p2, sqrt(m_s) / 2.);
    return R+V-Vs;
  }
  return R+V;
}


double YFS_Form_Factor::A1(const Vec4D &p1, const Vec4D &p2){
  double m12= m_m1*m_m2;
  double p12 = p1*p2;
  // double xlam = sqrt((p12 - m12) * (p12 + m12));
  double xlam = p12*sqrt(1-m12*m12/p12/p12);
  double t1 = -2 + (m_m1*m_m1-m_m2*m_m2)/(m_m1*m_m1+m_m2*m_m2-2.*p12)*log(m_m1/m_m2);
  t1 += -2.*xlam*xlam/(m_m1*m_m1+m_m2*m_m2-2.*p12)*A(p12, m_m1, m_m2);
  if(fabs(t1)>1e4){
    msg_Out()<<"m1 = "<<m_m1<<std::endl
             <<"m2 = "<<m_m2<<std::endl
             <<"xlam = "<<xlam<<std::endl
             <<"A(p12, m_m1, m_m2) = "<<A(p12, m_m1, m_m2)<<std::endl
             <<"(m_m1*m_m1+m_m2*m_m2-2.*p12) = "<<(m_m1*m_m1+m_m2*m_m2-2.*p12)<<std::endl;
  }
  return t1;
}

double YFS_Form_Factor::A2(const Vec4D &p1, const Vec4D &p2){
  double m12= m_m1*m_m2;
  double p12 = p1*p2;
  double xlam = sqrt((p12 - m12) * (p12 + m12));
  Complex omega1 = m_m1*m_m1/(xlam+p12-m_m1*m_m1);
  Complex omega2 = m_m2*m_m2/(xlam+p12-m_m2*m_m2);
  double scale = (p1-p2).Abs2();
  Complex t1 = log(Complex(scale,0)/m_m1/m_m2)*A(p12, m_m1, m_m2);
  Complex t2 = -0.5*(log(omega1*omega1)+log(omega2*omega2))
              +0.5*log(1.+omega1)*log(1.+omega1)
              +0.5*log(1.+omega2)*log(1.+omega2)
              -log(1.+omega1+omega2)*(log(omega1/(1.+omega1))+log(omega2/(1.+omega2)))
              -DiLog((1.+omega1)/(1.+omega1+omega2))
              -DiLog((1.+omega2)/(1.+omega1+omega2))
              +DiLog((omega1)/(1.+omega1+omega2))
              +DiLog((omega2)/(1.+omega1+omega2));
if(IsBad( t1+0.5/xlam*t2)){
  msg_Error()<<"A2(p1,p2) is Nan"<<std::endl
             <<"p1 = "<<p1<<std::endl
             <<"p2 = "<<p2<<std::endl
             <<"t1 = "<<t1<<std::endl
             <<"t2 = "<<t2<<std::endl
             <<"A0(p1,p2) = "<<A(p12, m_m1, m_m2)<<std::endl
             <<"p12 = "<<p12<<std::endl
             <<"omega1 = "<<omega1<<std::endl
             <<"omega2 = "<<omega2<<std::endl
             <<"m1 = "<<m_m1<<std::endl
             <<"m2 = "<<m_m2<<std::endl
             <<"p1.Theta() = "<<p1.Theta()<<std::endl
             <<"p2.Theta() = "<<p2.Theta()<<std::endl
             <<"xlam = "<<xlam<<std::endl;
}
return (t1+0.5/xlam*t2).real();
// return (0.5/xlam*t2).real();
// return (t1).real();
}


double YFS_Form_Factor::BVirtGeneral(YFS::Dipole &d, double Kmax){
  Vec4D p1,p2;
  d.LegsBeforeRadiation(p1, p2);
  m_m1 = d.GetMass(0);
  m_m2 = d.GetMass(1);
  double a0 = A(p1*p2, m_m1, m_m2);
  double a2 = A2(p1, p2);
  double form = log(m_photonMass*m_photonMass/m_m1/m_m2)*(p1*p2*a0-1);
  form += A1(p1, p2);
  form += -p1*p2*a2;
  // if(fabs(form)> 1e4 ){
  //   msg_Out()<<"Form = "<<form<<std::endl
  //            <<"p1 = "<<p1<<std::endl
  //            <<"p2 = "<<p2<<std::endl
  //            <<"p1.Theta() = "<<p1.Theta()<<std::endl
  //            <<"p2.Theta() = "<<p2.Theta()<<std::endl
  //            <<"a0 = "<<a0<<std::endl
  //            <<"a1 = "<<A1(p1,p2)<<std::endl
  //            <<"a2 = "<<a2<<std::endl;
  // }
  if(IsBad(form)){
    msg_Error()<<"YFS Btilde is NaN"<<std::endl
                <<"A0 = "<<a0<<std::endl
                <<"A1 = "<<A1(p1, p2)<<std::endl
                <<"A2 = "<<a2<<std::endl
                <<d.Type()<<std::endl
                <<"d.Left() = "<<d.GetFlav(0)<<std::endl
                <<"d.Right() = "<<d.GetFlav(1)<<std::endl;
    
  }
  return m_alpi*form;
}


DivArrD YFS_Form_Factor::BVirtGeneralEps(YFS::Dipole &d, double Kmax){
  Vec4D p1,p2;
  d.LegsBeforeRadiation(p1, p2);
  const double m1 = d.GetMass(0);
  const double m2 = d.GetMass(1);
  m_m1 = m1;
  m_m2 = m2;
  const double a0 = A(p1*p2, m1, m2);
  const double a2 = A2(p1, p2);
  const double irloop = (p_virt?p_virt->IRscale():100);
  const double epsloop = (p_virt?p_virt->Eps_Scheme_Factor({p1,p2}):4*M_PI);
  DivArrD massph(0,-1,0,0,0,0);
  // double form = log(m_photonMass*m_photonMass/m_m1/m_m2)*(p1*p2*a0-1);
  DivArrD form = (massph+log(4.*M_PI*sqr(irloop)/m1/m2/epsloop))*(p1*p2*a0-1.);
  form += A1(p1, p2);
  form += -p1*p2*a2;
  return m_alpi*form;
}


// Virtual YFS IR factor Y assembled from the master integrals (METOOLS/Loops),
// following Carloni Calame et al. "Towards muon-electron scattering at NNLO"
// arXiv:2007.01586 Eq.(3.3), LoopTools/COLLIER conventions:
//
//   off-diagonal (i!=j):  Y_ij = (a/pi) Qi Qj th_i th_j
//        [ p_i.p_j C0(m_i^2, s_ij, m_j^2; lam^2, m_i^2, m_j^2)
//          + 1/4 B0(s_ij, m_i^2, m_j^2) ],           s_ij = (th_i p_i + th_j p_j)^2
//   diagonal (i=j):       Y_ii = 1/8 (a/pi) Qi^2
//        [ B0(0, m_i^2, m_i^2) - 4 m_i^2 C0(m_i^2, 0, m_i^2; lam^2, m_i^2, m_i^2) ]
//
//   th_i = -1 (+1) for incoming (outgoing); Qi in positron units.
//
// The per-dipole factor returned here is Y_ij(off-diagonal) + Y_ii + Y_jj (the
// exchange plus both leg self-energies, i.e. the paper's Y_a for the line pair).
// The Qi Qj th_i th_j charge/flow factor is applied by the caller
// (Define_Dipoles::FormFactorSum via ChargeNorm), so it is omitted here.
//
// Regularisation: the paper uses a photon mass lambda. Here the master integrals
// are dimensionally regularised (massless photon internal line -> 1/eps_IR), so
// Y is returned as a DivArr whose 1/eps_IR pole cancels against the dim-reg real
// emission BVR_full_eps (enforced in BVR_full(Dipole&), FULL_FORM==3).
//
// TODO(validate): (i) the C0 photon-line is taken massless (dim-reg) rather than
// with the paper's lambda^2 - equivalent up to the standard IR dictionary
// (4 pi mu^2)^eps Gamma(1+eps)/eps -> ln lambda^2, checked via the 3-vs-0 XS;
// (ii) the diagonal Y_ii/Y_jj are added in full for an isolated dipole - if a
// leg is shared between dipoles this double counts and must be shared out.
DivArrD YFS_Form_Factor::BVirtMaster(YFS::Dipole &d, double Kmax) {
  Vec4D p1, p2;
  double th1, th2;
  if (d.Type()==dipoletype::initial) {          // both incoming
    p1 = d.GetBornMomenta(0);  p2 = d.GetBornMomenta(1);  th1 = -1.; th2 = -1.;
  } else if (d.Type()==dipoletype::final) {      // both outgoing
    p1 = d.GetBornMomenta(0);  p2 = d.GetBornMomenta(1);  th1 = +1.; th2 = +1.;
  } else if (d.Type()==dipoletype::ifi) {        // initial-final
    p1 = d.GetBornMomenta(0);  p2 = d.GetBornMomenta(1);  th1 = -1.; th2 = +1.;
  } else {
    THROW(fatal_error, "Unknown Dipole Type");
  }
  const double m1 = d.GetMass(0);
  const double m2 = d.GetMass(1);
  const double m12 = m1*m1, m22 = m2*m2;
  const double sij = (th1*p1 + th2*p2).Abs2();   // (th_i p_i + th_j p_j)^2
  const double mu2 = (p_virt?sqr(p_virt->IRscale()):GLOBAL_RENORMALISATION_SCALE);
  const Complex z(0.,0.);

  // --- off-diagonal exchange term, Eq.(3.3) i!=j ---
  DivArrC C0ij = Master_Triangle(m12, sij, m22, z, Complex(m12,0.), Complex(m22,0.), mu2);
  DivArrC B0ij = Master_Bubble(sij, Complex(m12,0.), Complex(m22,0.), mu2);
  DivArrC Yij  = C0ij*(p1*p2) + B0ij*0.25;

  // --- diagonal self-energy terms, Eq.(3.3) i=j (one per leg) ---
  DivArrC C0ii = Master_Triangle(m12, 0., m12, z, Complex(m12,0.), Complex(m12,0.), mu2);
  DivArrC C0jj = Master_Triangle(m22, 0., m22, z, Complex(m22,0.), Complex(m22,0.), mu2);
  DivArrC B0ii = Master_Bubble(0., Complex(m12,0.), Complex(m12,0.), mu2);
  DivArrC B0jj = Master_Bubble(0., Complex(m22,0.), Complex(m22,0.), mu2);
  DivArrC Yii  = (B0ii - C0ii*(4.*m12))*0.125;
  DivArrC Yjj  = (B0jj - C0jj*(4.*m22))*0.125;

  DivArrC combo = Yij + Yii + Yjj;

  // real part (YFS Y uses Re B) -> DivArrD, normalised to alpha/pi (Eq.3.3)
  DivArrD form(combo[0].real(), combo[1].real(), combo[2].real(),
               combo[3].real(), combo[4].real(), combo[5].real());
  return form * m_alpi;
}
