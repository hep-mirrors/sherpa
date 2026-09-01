#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/FormFactor.H"
#include "METOOLS/Loops/Master_Integrals.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "clooptools.h"
#include <iostream>

using namespace PHASIC;
using namespace ATOOLS;
using namespace METOOLS;
using namespace MODEL;

namespace EXTRAXS {
  class PionPionVirtual: public PHASIC::Virtual_ME2_Base {
    double m_eps2, m_eps, m_fin, m_s, m_t;
    double ME2, MP2, m_alpha;

  public:

    PionPionVirtual(const Process_Info& pi, const Flavour_Vector& flavs,
                    const double& ep2, const double& ep);

    ~PionPionVirtual() {}

    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    Complex BornTriangle();
    Complex BornBox();
    Complex BornSelf();
    Complex Full();
    double Kappa();
    double Kappat(const double t, const double u);
    Complex ISR();
    Complex FSR();
    Complex IFI();
    Complex C0e(const double x, const double y, const double mass);
    Complex C0e(const double x, const double y, const double mass1,
                const double mass2);
    Complex C0e(const double x, const double mass);
    Complex D0e(const double z, const double x, const double y,
                const double mass1, const double mass2);
    inline double Den(const double& a, const double& b)
    {
      return (a == b ? 0. : 1. / (a - b));
    }
    std::unique_ptr<FormFactor> p_formfactor;

    // Complex MBub(double a, double b, double c);
    // Complex MTri(double a, double b, double c,
    //             double aa, double bb, double cc);
    // Complex MBox(double a, double b, double c);
    // inline Complex Conjugate(Complex a) {return conj(a);}
    // inline double Re(Complex a) {return a.real();}
  };
} // namespace EXTRAXS

using namespace EXTRAXS;

PionPionVirtual::PionPionVirtual(const Process_Info& pi,
                                 const Flavour_Vector& flavs, const double& ep2,
                                 const double& ep)
    : Virtual_ME2_Base(pi, flavs), m_eps2(ep2), m_eps(ep)
{
  // if (!s_loader->LoadLibrary("ooptools")) THROW(fatal_error, "Failed to load
  // libooptools.");
  ltini();
  Setlambda(0.);
  ME = Flavour(kf_e).Mass();
  MM = (Flavour(kf_mu).Mass());
  ML = (Flavour(kf_tau).Mass());
  MP = (Flavour(kf_pi).Mass());
  ME2 = (Flavour(kf_e).Mass() * Flavour(kf_e).Mass());
  MM2 = (Flavour(kf_mu).Mass() * Flavour(kf_mu).Mass());
  ML2 = (Flavour(kf_tau).Mass() * Flavour(kf_tau).Mass());
  MP2 = (Flavour(kf_pi).Mass() * Flavour(kf_pi).Mass());
  MZ = Flavour(kf_Z).Mass();
  MZ2 = MZ * MZ;
  double MW = Flavour(kf_Wplus).Mass();
  double MW2 = MW * MW;

  MP2 = (Flavour(kf_pi).Mass() * Flavour(kf_pi).Mass());
  SW2 = std::abs(MODEL::s_model->ComplexConstant(string("csin2_thetaW")));
  CW2 = 1. - SW * SW;
  CW = sqrt(CW2);
  SW = sqrt(SW2);
  m_mode = 2;
  m_IRscale = 1.;
  m_UVscale = 1.;
  p_formfactor = std::unique_ptr<FormFactor>(new FormFactor());
  Scoped_Settings s{Settings::GetMainSettings()["YFS"]};
  m_photonmass = s["PHOTON_MASS"].Get<double>();
  int dimreg = s["Dim_Reg"].Get<int>();
  if (dimreg) m_photonmass = 0;
  setcmpbits(64);
  // setmudim(m_IRscale);
  // setdelta(4.*M_PI);
  // Setminmass(0.);
}

// ---- numbers against the Mathematica values printed below.
void PionPionVirtual::MICrossCheck() const
{
  const double s = 1.0, me2 = ME2, mp2 = MP2, mu2 = 1.0;
  msg_Out() << "B0 (s,me2,me2)   = " << B0_(s, me2, me2, mu2) << std::endl;
  msg_Out() << "C0 (me2,me2,s)   = \
" << C0_(me2, me2, s, 0., me2, me2, mu2)
            << std::endl;
  msg_Out() << "C00(me2,me2,s)   = \
" << C00_(me2, me2, s, 0., me2, me2, mu2)
            << std::endl;
  msg_Out() << "C1 (me2,me2,s)   = \
" << C1_(me2, me2, s, 0., me2, me2, mu2)
            << std::endl;
}


void PionPionVirtual::Calc(const Vec4D_Vector& momenta)
{
  double factor(1.);
  if (m_stype & sbt::qcd)
    factor = 2 * M_PI / AlphaQCD();
  else if (m_stype & sbt::qed)
    factor = 2 * M_PI / AlphaQED();
  else
    THROW(fatal_error, "Unknown coupling.");
  // 1/epsIR
  // m_res.IR()=m_eps*factor;
  // // 1/epsIR2
  // m_res.IR2()=m_eps2*factor;
  // // finite
  // if(m_count==1000){
  // m_count=0;
  clearcache(); // Since in general s will be different we gain nothing from
                // cache
  m_s = (momenta[0] + momenta[1]).Abs2();
  m_alpha = (*aqed)(m_s);
  m_res.Finite() = 1. / m_s / m_s * 2 * (BornBox() + BornTriangle()).real();
  // FORTRAN(ltexi)();
}

DivArrC PionPionVirtual::AmpISR() const
{
  return -(Power(Pi, 2) *
           Global_B0(Power(ME 2), 0,
                     Power(ME 2)) *
           Power(Sqrt(4.*M_PI*m_alpha), 6)) +
         (D * Power(Pi, 2) *
          Global_B0(Power(ME 2), 0,
                    Power(ME 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             2. +
         (Power(Pi, 2) * Power(t, 2) *
          Global_B0(Power(ME 2), 0,
                    Power(ME 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             Power(s, 2) -
         (D * Power(Pi, 2) * Power(t, 2) *
          Global_B0(Power(ME 2), 0,
                    Power(ME 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             (2. * Power(s, 2)) -
         (2 * Power(Pi, 2) * t * u *
          Global_B0(Power(ME 2), 0,
                    Power(ME 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             Power(s, 2) +
         (D * Power(Pi, 2) * t * u *
          Global_B0(Power(ME 2), 0,
                    Power(ME 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             Power(s, 2) +
         (Power(Pi, 2) * Power(u, 2) *
          Global_B0(Power(ME 2), 0,
                    Power(ME 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             Power(s, 2) -
         (D * Power(Pi, 2) * Power(u, 2) *
          Global_B0(Power(ME 2), 0,
                    Power(ME 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             (2. * Power(s, 2)) -
         Power(Pi, 2) * s *
             Global_PaVe(1,
                         List(s, Power(ME 2),
                              Power(ME 2)),
                         List(Power(ME 2),
                              Power(ME 2), 0),
                         Rule(Global_PaVeAutoOrder, True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) +
         (D * Power(Pi, 2) * s *
          Global_PaVe(
              1,
              List(s, Power(ME 2), Power(ME 2)),
              List(Power(ME 2), Power(ME 2), 0),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             2. +
         (Power(Pi, 2) * Power(t, 2) *
          Global_PaVe(
              1,
              List(s, Power(ME 2), Power(ME 2)),
              List(Power(ME 2), Power(ME 2), 0),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             s -
         (D * Power(Pi, 2) * Power(t, 2) *
          Global_PaVe(
              1,
              List(s, Power(ME 2), Power(ME 2)),
              List(Power(ME 2), Power(ME 2), 0),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             (2. * s) -
         (2 * Power(Pi, 2) * t * u *
          Global_PaVe(
              1,
              List(s, Power(ME 2), Power(ME 2)),
              List(Power(ME 2), Power(ME 2), 0),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             s +
         (D * Power(Pi, 2) * t * u *
          Global_PaVe(
              1,
              List(s, Power(ME 2), Power(ME 2)),
              List(Power(ME 2), Power(ME 2), 0),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             s +
         (Power(Pi, 2) * Power(u, 2) *
          Global_PaVe(
              1,
              List(s, Power(ME 2), Power(ME 2)),
              List(Power(ME 2), Power(ME 2), 0),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             s -
         (D * Power(Pi, 2) * Power(u, 2) *
          Global_PaVe(
              1,
              List(s, Power(ME 2), Power(ME 2)),
              List(Power(ME 2), Power(ME 2), 0),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             (2. * s) +
         Power(Pi, 2) * s *
             Global_PaVe(1,
                         List(Power(ME 2),
                              Power(ME 2), s),
                         List(Power(ME 2), 0,
                              Power(ME 2)),
                         Rule(Global_PaVeAutoOrder, True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) -
         (Power(Pi, 2) * Power(t, 2) *
          Global_PaVe(
              1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             s +
         (2 * Power(Pi, 2) * t * u *
          Global_PaVe(
              1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             s -
         (Power(Pi, 2) * Power(u, 2) *
          Global_PaVe(
              1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             s +
         2 * Power(Pi, 2) *
             Global_PaVe(0, 0,
                         List(s, Power(ME 2),
                              Power(ME 2)),
                         List(Power(ME 2),
                              Power(ME 2), 0),
                         Rule(Global_PaVeAutoOrder, True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) -
         D * Power(Pi, 2) *
             Global_PaVe(0, 0,
                         List(s, Power(ME 2),
                              Power(ME 2)),
                         List(Power(ME 2),
                              Power(ME 2), 0),
                         Rule(Global_PaVeAutoOrder, True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) -
         (2 * Power(Pi, 2) * Power(t, 2) *
          Global_PaVe(
              0, 0,
              List(s, Power(ME 2), Power(ME 2)),
              List(Power(ME 2), Power(ME 2), 0),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             Power(s, 2) +
         (D * Power(Pi, 2) * Power(t, 2) *
          Global_PaVe(
              0, 0,
              List(s, Power(ME 2), Power(ME 2)),
              List(Power(ME 2), Power(ME 2), 0),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             Power(s, 2) +
         (4 * Power(Pi, 2) * t * u *
          Global_PaVe(
              0, 0,
              List(s, Power(ME 2), Power(ME 2)),
              List(Power(ME 2), Power(ME 2), 0),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             Power(s, 2) -
         (2 * D * Power(Pi, 2) * t * u *
          Global_PaVe(
              0, 0,
              List(s, Power(ME 2), Power(ME 2)),
              List(Power(ME 2), Power(ME 2), 0),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             Power(s, 2) -
         (2 * Power(Pi, 2) * Power(u, 2) *
          Global_PaVe(
              0, 0,
              List(s, Power(ME 2), Power(ME 2)),
              List(Power(ME 2), Power(ME 2), 0),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             Power(s, 2) +
         (D * Power(Pi, 2) * Power(u, 2) *
          Global_PaVe(
              0, 0,
              List(s, Power(ME 2), Power(ME 2)),
              List(Power(ME 2), Power(ME 2), 0),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             Power(s, 2) +
         (4 * Power(Pi, 2) * t *
          Global_C0(s, Power(ME 2), Power(ME 2),
                    Power(ME 2), Power(ME 2),
                    0) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s -
         (D * Power(Pi, 2) * t *
          Global_C0(s, Power(ME 2), Power(ME 2),
                    Power(ME 2), Power(ME 2),
                    0) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s +
         (4 * Power(Pi, 2) * u *
          Global_C0(s, Power(ME 2), Power(ME 2),
                    Power(ME 2), Power(ME 2),
                    0) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s -
         (D * Power(Pi, 2) * u *
          Global_C0(s, Power(ME 2), Power(ME 2),
                    Power(ME 2), Power(ME 2),
                    0) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s +
         (6 * Power(Pi, 2) * t *
          Global_PaVe(
              1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s +
         (2 * D * Power(Pi, 2) * Power(t, 2) *
          Global_PaVe(
              1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             Power(s, 2) +
         (6 * Power(Pi, 2) * u *
          Global_PaVe(
              1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s -
         (2 * D * Power(Pi, 2) * u *
          Global_PaVe(
              1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s -
         (2 * D * Power(Pi, 2) * t * u *
          Global_PaVe(
              1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             Power(s, 2) -
         (2 * Power(Pi, 2) * t *
          Global_PaVe(
              1, 1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s +
         (D * Power(Pi, 2) * t *
          Global_PaVe(
              1, 1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s -
         (4 * Power(Pi, 2) * Power(t, 2) *
          Global_PaVe(
              1, 1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             Power(s, 2) +
         (2 * D * Power(Pi, 2) * Power(t, 2) *
          Global_PaVe(
              1, 1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             Power(s, 2) +
         (2 * Power(Pi, 2) * u *
          Global_PaVe(
              1, 1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s -
         (D * Power(Pi, 2) * u *
          Global_PaVe(
              1, 1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s +
         (4 * Power(Pi, 2) * t * u *
          Global_PaVe(
              1, 1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             Power(s, 2) -
         (2 * D * Power(Pi, 2) * t * u *
          Global_PaVe(
              1, 1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             Power(s, 2) -
         (8 * Power(Pi, 2) *
          Global_C0(s, Power(ME 2), Power(ME 2),
                    Power(ME 2), Power(ME 2),
                    0) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 4)) /
             s +
         (2 * D * Power(Pi, 2) *
          Global_C0(s, Power(ME 2), Power(ME 2),
                    Power(ME 2), Power(ME 2),
                    0) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 4)) /
             s -
         (12 * Power(Pi, 2) *
          Global_PaVe(
              1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 4)) /
             s +
         (2 * D * Power(Pi, 2) *
          Global_PaVe(
              1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 4)) /
             s -
         (2 * D * Power(Pi, 2) * t *
          Global_PaVe(
              1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 4)) /
             Power(s, 2) +
         (2 * D * Power(Pi, 2) * u *
          Global_PaVe(
              1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 4)) /
             Power(s, 2) +
         (4 * Power(Pi, 2) * t *
          Global_PaVe(
              1, 1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 4)) /
             Power(s, 2) -
         (2 * D * Power(Pi, 2) * t *
          Global_PaVe(
              1, 1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 4)) /
             Power(s, 2) -
         (4 * Power(Pi, 2) * u *
          Global_PaVe(
              1, 1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 4)) /
             Power(s, 2) +
         (2 * D * Power(Pi, 2) * u *
          Global_PaVe(
              1, 1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 4)) /
             Power(s, 2) +
         (4 * Power(Pi, 2) *
          Global_B0(Power(ME 2), 0,
                    Power(ME 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2)) /
             s -
         (2 * D * Power(Pi, 2) *
          Global_B0(Power(ME 2), 0,
                    Power(ME 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2)) /
             s +
         4 * Power(Pi, 2) *
             Global_PaVe(1,
                         List(s, Power(ME 2),
                              Power(ME 2)),
                         List(Power(ME 2),
                              Power(ME 2), 0),
                         Rule(Global_PaVeAutoOrder, True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2) -
         2 * D * Power(Pi, 2) *
             Global_PaVe(1,
                         List(s, Power(ME 2),
                              Power(ME 2)),
                         List(Power(ME 2),
                              Power(ME 2), 0),
                         Rule(Global_PaVeAutoOrder, True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2) -
         4 * Power(Pi, 2) *
             Global_PaVe(1,
                         List(Power(ME 2),
                              Power(ME 2), s),
                         List(Power(ME 2), 0,
                              Power(ME 2)),
                         Rule(Global_PaVeAutoOrder, True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2) -
         (8 * Power(Pi, 2) *
          Global_PaVe(
              0, 0,
              List(s, Power(ME 2), Power(ME 2)),
              List(Power(ME 2), Power(ME 2), 0),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2)) /
             s +
         (4 * D * Power(Pi, 2) *
          Global_PaVe(
              0, 0,
              List(s, Power(ME 2), Power(ME 2)),
              List(Power(ME 2), Power(ME 2), 0),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2)) /
             s -
         (8 * Power(Pi, 2) *
          Global_C0(s, Power(ME 2), Power(ME 2),
                    Power(ME 2), Power(ME 2),
                    0) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
          Power(Global_SMP("m_pi"), 2)) /
             s +
         (2 * D * Power(Pi, 2) *
          Global_C0(s, Power(ME 2), Power(ME 2),
                    Power(ME 2), Power(ME 2),
                    0) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
          Power(Global_SMP("m_pi"), 2)) /
             s -
         (16 * Power(Pi, 2) * t *
          Global_C0(s, Power(ME 2), Power(ME 2),
                    Power(ME 2), Power(ME 2),
                    0) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
          Power(Global_SMP("m_pi"), 2)) /
             Power(s, 2) +
         (4 * D * Power(Pi, 2) * t *
          Global_C0(s, Power(ME 2), Power(ME 2),
                    Power(ME 2), Power(ME 2),
                    0) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
          Power(Global_SMP("m_pi"), 2)) /
             Power(s, 2) -
         (16 * Power(Pi, 2) * u *
          Global_C0(s, Power(ME 2), Power(ME 2),
                    Power(ME 2), Power(ME 2),
                    0) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
          Power(Global_SMP("m_pi"), 2)) /
             Power(s, 2) +
         (4 * D * Power(Pi, 2) * u *
          Global_C0(s, Power(ME 2), Power(ME 2),
                    Power(ME 2), Power(ME 2),
                    0) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
          Power(Global_SMP("m_pi"), 2)) /
             Power(s, 2) -
         (12 * Power(Pi, 2) *
          Global_PaVe(
              1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
          Power(Global_SMP("m_pi"), 2)) /
             s +
         (2 * D * Power(Pi, 2) *
          Global_PaVe(
              1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
          Power(Global_SMP("m_pi"), 2)) /
             s -
         (24 * Power(Pi, 2) * t *
          Global_PaVe(
              1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
          Power(Global_SMP("m_pi"), 2)) /
             Power(s, 2) +
         (2 * D * Power(Pi, 2) * t *
          Global_PaVe(
              1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
          Power(Global_SMP("m_pi"), 2)) /
             Power(s, 2) -
         (24 * Power(Pi, 2) * u *
          Global_PaVe(
              1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
          Power(Global_SMP("m_pi"), 2)) /
             Power(s, 2) +
         (6 * D * Power(Pi, 2) * u *
          Global_PaVe(
              1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
          Power(Global_SMP("m_pi"), 2)) /
             Power(s, 2) +
         (4 * Power(Pi, 2) * t *
          Global_PaVe(
              1, 1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
          Power(Global_SMP("m_pi"), 2)) /
             Power(s, 2) -
         (2 * D * Power(Pi, 2) * t *
          Global_PaVe(
              1, 1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
          Power(Global_SMP("m_pi"), 2)) /
             Power(s, 2) -
         (4 * Power(Pi, 2) * u *
          Global_PaVe(
              1, 1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
          Power(Global_SMP("m_pi"), 2)) /
             Power(s, 2) +
         (2 * D * Power(Pi, 2) * u *
          Global_PaVe(
              1, 1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
          Power(Global_SMP("m_pi"), 2)) /
             Power(s, 2) +
         (32 * Power(Pi, 2) *
          Global_C0(s, Power(ME 2), Power(ME 2),
                    Power(ME 2), Power(ME 2),
                    0) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 4) *
          Power(Global_SMP("m_pi"), 2)) /
             Power(s, 2) -
         (8 * D * Power(Pi, 2) *
          Global_C0(s, Power(ME 2), Power(ME 2),
                    Power(ME 2), Power(ME 2),
                    0) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 4) *
          Power(Global_SMP("m_pi"), 2)) /
             Power(s, 2) +
         (48 * Power(Pi, 2) *
          Global_PaVe(
              1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 4) *
          Power(Global_SMP("m_pi"), 2)) /
             Power(s, 2) -
         (8 * D * Power(Pi, 2) *
          Global_PaVe(
              1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 4) *
          Power(Global_SMP("m_pi"), 2)) /
             Power(s, 2) +
         (32 * Power(Pi, 2) *
          Global_C0(s, Power(ME 2), Power(ME 2),
                    Power(ME 2), Power(ME 2),
                    0) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
          Power(Global_SMP("m_pi"), 4)) /
             Power(s, 2) -
         (8 * D * Power(Pi, 2) *
          Global_C0(s, Power(ME 2), Power(ME 2),
                    Power(ME 2), Power(ME 2),
                    0) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
          Power(Global_SMP("m_pi"), 4)) /
             Power(s, 2) +
         (48 * Power(Pi, 2) *
          Global_PaVe(
              1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
          Power(Global_SMP("m_pi"), 4)) /
             Power(s, 2) -
         (8 * D * Power(Pi, 2) *
          Global_PaVe(
              1,
              List(Power(ME 2), Power(ME 2), s),
              List(Power(ME 2), 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
          Power(Global_SMP("m_pi"), 4)) /
             Power(s, 2);
}

DivArrC PionPionVirtual::AmpFSR() const
{
  return -(Power(Pi, 2) *
           Global_B0(Power(Global_SMP("m_pi"), 2), 0,
                     Power(Global_SMP("m_pi"), 2)) *
           Power(Sqrt(4.*M_PI*m_alpha), 6)) +
         (Power(Pi, 2) * Power(t, 2) *
          Global_B0(Power(Global_SMP("m_pi"), 2), 0,
                    Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             Power(s, 2) -
         (2 * Power(Pi, 2) * t * u *
          Global_B0(Power(Global_SMP("m_pi"), 2), 0,
                    Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             Power(s, 2) +
         (Power(Pi, 2) * Power(u, 2) *
          Global_B0(Power(Global_SMP("m_pi"), 2), 0,
                    Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             Power(s, 2) -
         Power(Pi, 2) * s *
             Global_C0(s, Power(Global_SMP("m_pi"), 2),
                       Power(Global_SMP("m_pi"), 2),
                       Power(Global_SMP("m_pi"), 2),
                       Power(Global_SMP("m_pi"), 2), 0) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) +
         (Power(Pi, 2) * Power(t, 2) *
          Global_C0(s, Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), 0) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             s -
         (2 * Power(Pi, 2) * t * u *
          Global_C0(s, Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), 0) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             s +
         (Power(Pi, 2) * Power(u, 2) *
          Global_C0(s, Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), 0) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             s -
         2 * Power(Pi, 2) * s *
             Global_PaVe(1,
                         List(s, Power(Global_SMP("m_pi"), 2),
                              Power(Global_SMP("m_pi"), 2)),
                         List(Power(Global_SMP("m_pi"), 2),
                              Power(Global_SMP("m_pi"), 2), 0),
                         Rule(Global_PaVeAutoOrder, True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) +
         (2 * Power(Pi, 2) * Power(t, 2) *
          Global_PaVe(1,
                      List(s, Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2)),
                      List(Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2), 0),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             s -
         (4 * Power(Pi, 2) * t * u *
          Global_PaVe(1,
                      List(s, Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2)),
                      List(Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2), 0),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             s +
         (2 * Power(Pi, 2) * Power(u, 2) *
          Global_PaVe(1,
                      List(s, Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2)),
                      List(Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2), 0),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             s +
         (4 * Power(Pi, 2) *
          Global_B0(Power(Global_SMP("m_pi"), 2), 0,
                    Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2)) /
             s +
         6 * Power(Pi, 2) *
             Global_C0(s, Power(Global_SMP("m_pi"), 2),
                       Power(Global_SMP("m_pi"), 2),
                       Power(Global_SMP("m_pi"), 2),
                       Power(Global_SMP("m_pi"), 2), 0) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2) -
         (2 * Power(Pi, 2) * Power(t, 2) *
          Global_C0(s, Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), 0) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2)) /
             Power(s, 2) +
         (4 * Power(Pi, 2) * t * u *
          Global_C0(s, Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), 0) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2)) /
             Power(s, 2) -
         (2 * Power(Pi, 2) * Power(u, 2) *
          Global_C0(s, Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), 0) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2)) /
             Power(s, 2) +
         12 * Power(Pi, 2) *
             Global_PaVe(1,
                         List(s, Power(Global_SMP("m_pi"), 2),
                              Power(Global_SMP("m_pi"), 2)),
                         List(Power(Global_SMP("m_pi"), 2),
                              Power(Global_SMP("m_pi"), 2), 0),
                         Rule(Global_PaVeAutoOrder, True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2) -
         (4 * Power(Pi, 2) * Power(t, 2) *
          Global_PaVe(1,
                      List(s, Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2)),
                      List(Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2), 0),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2)) /
             Power(s, 2) +
         (8 * Power(Pi, 2) * t * u *
          Global_PaVe(1,
                      List(s, Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2)),
                      List(Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2), 0),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2)) /
             Power(s, 2) -
         (4 * Power(Pi, 2) * Power(u, 2) *
          Global_PaVe(1,
                      List(s, Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2)),
                      List(Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2), 0),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2)) /
             Power(s, 2) -
         (8 * Power(Pi, 2) *
          Global_C0(s, Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), 0) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 4)) /
             s -
         (16 * Power(Pi, 2) *
          Global_PaVe(1,
                      List(s, Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2)),
                      List(Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2), 0),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 4)) /
             s;
}

DivArrC PionPionVirtual::AmpIFI() const
{
  return -0.5 * (Power(Pi, 2) * s * t *
                 Global_D0(s, Power(ME 2), t,
                           Power(Global_SMP("m_pi"), 2),
                           Power(ME 2),
                           Power(Global_SMP("m_\
pi"),
                                 2),
                           0, 0, Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)) *
                 Power(Sqrt(4.*M_PI*m_alpha), 6)) +
         (Power(Pi, 2) * Power(t, 3) *
          Global_D0(s, Power(ME 2), t,
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_\
SMP("m_pi"),
                          2),
                    0, 0, Power(ME 2),
                    Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             (2. * s) -
         (Power(Pi, 2) * Power(t, 2) * u *
          Global_D0(s, Power(ME 2), t,
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_\
SMP("m_pi"),
                          2),
                    0, 0, Power(ME 2),
                    Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             s +
         (Power(Pi, 2) * t * Power(u, 2) *
          Global_D0(s, Power(ME 2), t,
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_\
SMP("m_pi"),
                          2),
                    0, 0, Power(ME 2),
                    Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             (2. * s) -
         (Power(Pi, 2) * Power(s, 2) *
          Global_D0(s, Power(ME 2),
                    -s - t + 2 * Power(ME 2) +
                        2 * Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_SMP("m_pi"), 2), 0, 0,
                    Power(ME 2), Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             2. -
         (Power(Pi, 2) * s * t *
          Global_D0(s, Power(ME 2),
                    -s - t + 2 * Power(ME 2) +
                        2 * Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_SMP("m_pi"), 2), 0, 0,
                    Power(ME 2), Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             2. +
         (Power(Pi, 2) * Power(t, 2) *
          Global_D0(s, Power(ME 2),
                    -s - t + 2 * Power(ME 2) +
                        2 * Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_SMP("m_pi"), 2), 0, 0,
                    Power(ME 2), Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             2. +
         (Power(Pi, 2) * Power(t, 3) *
          Global_D0(s, Power(ME 2),
                    -s - t + 2 * Power(ME 2) +
                        2 * Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_SMP("m_pi"), 2), 0, 0,
                    Power(ME 2), Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             (2. * s) -
         Power(Pi, 2) * t * u *
             Global_D0(
                 s, Power(ME 2),
                 -s - t + 2 * Power(ME 2) +
                     2 * Power(Global_SMP("m_pi"), 2),
                 Power(Global_SMP("m_pi"), 2), Power(ME 2),
                 Power(Global_SMP("m_pi"), 2), 0, 0,
                 Power(ME 2), Power(Global_SMP("m_pi"), 2)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) -
         (Power(Pi, 2) * Power(t, 2) * u *
          Global_D0(s, Power(ME 2),
                    -s - t + 2 * Power(ME 2) +
                        2 * Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_SMP("m_pi"), 2), 0, 0,
                    Power(ME 2), Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             s +
         (Power(Pi, 2) * Power(u, 2) *
          Global_D0(s, Power(ME 2),
                    -s - t + 2 * Power(ME 2) +
                        2 * Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_SMP("m_pi"), 2), 0, 0,
                    Power(ME 2), Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             2. +
         (Power(Pi, 2) * t * Power(u, 2) *
          Global_D0(s, Power(ME 2),
                    -s - t + 2 * Power(ME 2) +
                        2 * Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_SMP("m_pi"), 2), 0, 0,
                    Power(ME 2), Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             (2. * s) +
         Power(Pi, 2) * s * t *
             Global_PaVe(1,
                         List(s, Power(ME 2), t,
                              Power(Global_SMP("m_pi"), 2),
                              Power(ME 2),
                              Power(Global_\
SMP("m_pi"),
                                    2)),
                         List(0, 0, Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         Rule(Global_PaVeAutoOrder, True),
                         Rule(Global_\
PaVeAutoReduce,
                              True)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) -
         (Power(Pi, 2) * Power(t, 3) *
          Global_PaVe(1,
                      List(s, Power(ME 2), t,
                           Power(Global_SMP("m_pi"), 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      List(0, 0, Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(Global_\
PaVeAutoReduce,
                           True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             s +
         (2 * Power(Pi, 2) * Power(t, 2) * u *
          Global_PaVe(1,
                      List(s,
                           Power(Global_SMP("m_\
e"),
                                 2),
                           t, Power(Global_SMP("m_pi"), 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      List(0, 0, Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(Global_\
PaVeAutoReduce,
                           True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             s -
         (Power(Pi, 2) * t * Power(u, 2) *
          Global_PaVe(1,
                      List(s, Power(ME 2), t,
                           Power(Global_SMP("m_pi"), 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      List(0, 0, Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(Global_\
PaVeAutoReduce,
                           True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6)) /
             s +
         Power(Pi, 2) * Power(s, 2) *
             Global_PaVe(1,
                         List(s, Power(ME 2),
                              -s - t + 2 * Power(ME 2) +
                                  2 * Power(Global_SMP("m_pi"), 2),
                              Power(Global_SMP("m_pi"), 2),
                              Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         List(0, 0, Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         Rule(Global_\
PaVeAutoOrder,
                              True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Global_\
SMP("e"),
                   6) +
         Power(Pi, 2) * s * t *
             Global_PaVe(1,
                         List(s, Power(ME 2),
                              -s - t + 2 * Power(ME 2) +
                                  2 * Power(Global_SMP("m_pi"), 2),
                              Power(Global_SMP("m_pi"), 2),
                              Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         List(0, 0, Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         Rule(Global_\
PaVeAutoOrder,
                              True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Global_\
SMP("e"),
                   6) -
         Power(Pi, 2) * Power(t, 2) *
             Global_PaVe(1,
                         List(s, Power(ME 2),
                              -s - t + 2 * Power(ME 2) +
                                  2 * Power(Global_SMP("m_pi"), 2),
                              Power(Global_SMP("m_pi"), 2),
                              Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         List(0, 0, Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         Rule(Global_\
PaVeAutoOrder,
                              True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Global_\
SMP("e"),
                   6) -
         (Power(Pi, 2) * Power(t, 3) *
          Global_PaVe(1,
                      List(s, Power(ME 2),
                           -s - t + 2 * Power(ME 2) +
                               2 * Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      List(0, 0, Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_\
PaVeAutoOrder,
                           True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Global_\
SMP("e"),
                6)) /
             s +
         2 * Power(Pi, 2) * t * u *
             Global_PaVe(1,
                         List(s, Power(ME 2),
                              -s - t + 2 * Power(ME 2) +
                                  2 * Power(Global_SMP("m_pi"), 2),
                              Power(Global_SMP("m_pi"), 2),
                              Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         List(0, 0, Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         Rule(Global_\
PaVeAutoOrder,
                              True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Global_\
SMP("e"),
                   6) +
         (2 * Power(Pi, 2) * Power(t, 2) * u *
          Global_PaVe(1,
                      List(s,
                           Power(Global_SMP("m_\
e"),
                                 2),
                           -s - t + 2 * Power(ME 2) +
                               2 * Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      List(0, 0, Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_\
PaVeAutoOrder,
                           True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Global_\
SMP("e"),
                6)) /
             s -
         Power(Pi, 2) * Power(u, 2) *
             Global_PaVe(1,
                         List(s, Power(ME 2),
                              -s - t + 2 * Power(ME 2) +
                                  2 * Power(Global_SMP("m_pi"), 2),
                              Power(Global_SMP("m_pi"), 2),
                              Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         List(0, 0, Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         Rule(Global_\
PaVeAutoOrder,
                              True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Global_\
SMP("e"),
                   6) -
         (Power(Pi, 2) * t * Power(u, 2) *
          Global_PaVe(1,
                      List(s, Power(ME 2),
                           -s - t + 2 * Power(ME 2) +
                               2 * Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      List(0, 0, Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_\
PaVeAutoOrder,
                           True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Global_\
SMP("e"),
                6)) /
             s -
         (2 * D * Power(Pi, 2) * t *
          Global_C0(s, Power(ME 2), Power(ME 2),
                    0, 0, Power(ME 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s +
         (2 * D * Power(Pi, 2) * u *
          Global_C0(s, Power(ME 2), Power(ME 2),
                    0, 0, Power(ME 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s +
         (Power(Pi, 2) * s *
          Global_D0(s, Power(ME 2), t,
                    Power(Global_\
SMP("m_pi"),
                          2),
                    Power(ME 2), Power(Global_SMP("m_pi"), 2),
                    0, 0, Power(ME 2),
                    Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             2. -
         (Power(Pi, 2) * Power(t, 2) *
          Global_D0(s, Power(ME 2), t,
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_\
SMP("m_pi"),
                          2),
                    0, 0, Power(ME 2),
                    Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             (2. * s) +
         (Power(Pi, 2) * t * u *
          Global_D0(s, Power(ME 2), t,
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_SMP("m_\
pi"),
                          2),
                    0, 0, Power(ME 2),
                    Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s -
         (Power(Pi, 2) * Power(u, 2) *
          Global_D0(s, Power(ME 2), t,
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_\
SMP("m_pi"),
                          2),
                    0, 0, Power(ME 2),
                    Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             (2. * s) +
         (Power(Pi, 2) * s *
          Global_D0(s, Power(ME 2),
                    -s - t + 2 * Power(ME 2) +
                        2 * Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_SMP("m_pi"), 2), 0, 0,
                    Power(ME 2), Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             2. -
         (Power(Pi, 2) * Power(t, 2) *
          Global_D0(s, Power(ME 2),
                    -s - t + 2 * Power(ME 2) +
                        2 * Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_SMP("m_pi"), 2), 0, 0,
                    Power(ME 2), Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             (2. * s) +
         (Power(Pi, 2) * t * u *
          Global_D0(s, Power(ME 2),
                    -s - t + 2 * Power(ME 2) +
                        2 * Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_SMP("m_pi"), 2), 0, 0,
                    Power(ME 2), Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s -
         (Power(Pi, 2) * Power(u, 2) *
          Global_D0(s, Power(ME 2),
                    -s - t + 2 * Power(ME 2) +
                        2 * Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_SMP("m_pi"), 2), 0, 0,
                    Power(ME 2), Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             (2. * s) +
         (8 * Power(Pi, 2) * t *
          Global_PaVe(
              1,
              List(s, Power(ME 2), Power(ME 2)),
              List(0, 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s -
         (4 * D * Power(Pi, 2) * t *
          Global_PaVe(
              1,
              List(s, Power(ME 2), Power(ME 2)),
              List(0, 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s -
         (8 * Power(Pi, 2) * u *
          Global_PaVe(
              1,
              List(s, Power(ME 2), Power(ME 2)),
              List(0, 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s +
         (4 * D * Power(Pi, 2) * u *
          Global_PaVe(
              1,
              List(s, Power(ME 2), Power(ME 2)),
              List(0, 0, Power(ME 2)),
              Rule(Global_PaVeAutoOrder, True),
              Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s -
         Power(Pi, 2) * s *
             Global_PaVe(1,
                         List(s, Power(ME 2), t,
                              Power(Global_SMP("m_pi"), 2),
                              Power(ME 2),
                              Power(
                                  Global_SMP("m_\
pi"),
                                  2)),
                         List(0, 0, Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         Rule(Global_PaVeAutoOrder, True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) +
         (Power(Pi, 2) * Power(t, 2) *
          Global_PaVe(1,
                      List(s, Power(ME 2), t,
                           Power(Global_SMP("m_pi"), 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      List(0, 0, Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(
                          Global_\
PaVeAutoReduce,
                          True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s -
         (2 * Power(Pi, 2) * t * u *
          Global_PaVe(1,
                      List(s, Power(ME 2), t,
                           Power(Global_SMP("m_pi"), 2),
                           Power(ME 2),
                           Power(
                               Global_\
SMP("m_pi"),
                               2)),
                      List(0, 0, Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(
                          Global_\
PaVeAutoReduce,
                          True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s +
         (Power(Pi, 2) * Power(u, 2) *
          Global_PaVe(1,
                      List(s, Power(ME 2), t,
                           Power(Global_SMP("m_pi"), 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      List(0, 0, Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(
                          Global_\
PaVeAutoReduce,
                          True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s -
         Power(Pi, 2) * s *
             Global_PaVe(1,
                         List(s, Power(ME 2),
                              -s - t + 2 * Power(ME 2) +
                                  2 * Power(Global_SMP("m_pi"), 2),
                              Power(Global_SMP("m_pi"), 2),
                              Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         List(0, 0, Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         Rule(Global_\
PaVeAutoOrder,
                              True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Global_\
SMP("e"),
                   6) *
             Power(ME 2) +
         (Power(Pi, 2) * Power(t, 2) *
          Global_PaVe(1,
                      List(s, Power(ME 2),
                           -s - t + 2 * Power(ME 2) +
                               2 * Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      List(0, 0, Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_\
PaVeAutoOrder,
                           True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Global_\
SMP("e"),
                6) *
          Power(ME 2)) /
             s -
         (2 * Power(Pi, 2) * t * u *
          Global_PaVe(1,
                      List(s, Power(ME 2),
                           -s - t + 2 * Power(ME 2) +
                               2 * Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      List(0, 0, Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_\
PaVeAutoOrder,
                           True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Global_\
SMP("e"),
                6) *
          Power(ME 2)) /
             s +
         (Power(Pi, 2) * Power(u, 2) *
          Global_PaVe(1,
                      List(s, Power(ME 2),
                           -s - t + 2 * Power(ME 2) +
                               2 * Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      List(0, 0, Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_\
PaVeAutoOrder,
                           True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Global_\
SMP("e"),
                6) *
          Power(ME 2)) /
             s +
         2 * Power(Pi, 2) * t *
             Global_PaVe(1,
                         List(t, Power(ME 2), s,
                              Power(Global_SMP("m_pi"), 2),
                              Power(Global_SMP("m_pi"), 2),
                              Power(
                                  Global_\
SMP("m_e"),
                                  2)),
                         List(Power(Global_SMP("m_pi"), 2),
                              Power(Global_SMP("m_\
e"),
                                    2),
                              0, 0),
                         Rule(Global_PaVeAutoOrder, True),
                         Rule(
                             Global_\
PaVeAutoReduce,
                             True)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) -
         (2 * Power(Pi, 2) * Power(t, 2) *
          Global_PaVe(1,
                      List(t, Power(ME 2), s,
                           Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(ME 2)),
                      List(Power(Global_SMP("m_pi"), 2),
                           Power(ME 2), 0, 0),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(
                          Global_\
PaVeAutoReduce,
                          True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s +
         (2 * Power(Pi, 2) * t * u *
          Global_PaVe(1,
                      List(t, Power(ME 2), s,
                           Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(
                               Global_\
SMP("m_e"),
                               2)),
                      List(Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_\
e"),
                                 2),
                           0, 0),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(
                          Global_\
PaVeAutoReduce,
                          True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2)) /
             s +
         2 * Power(Pi, 2) * t *
             Global_PaVe(1,
                         List(Power(ME 2),
                              Power(ME 2),
                              Power(Global_SMP("m_pi"), 2),
                              Power(
                                  Global_SMP("m_\
pi"),
                                  2),
                              s,
                              -s - t + 2 * Power(ME 2) +
                                  2 * Power(Global_SMP("m_pi"), 2)),
                         List(0, Power(ME 2), 0,
                              Power(Global_SMP("m_pi"), 2)),
                         Rule(Global_PaVeAutoOrder, True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) *
             Power(
                 Global_\
SMP("m_e"),
                 2) +
         (2 * Power(Pi, 2) * Power(t, 2) *
          Global_PaVe(1,
                      List(Power(ME 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2), s,
                           -s - t + 2 * Power(ME 2) +
                               2 * Power(Global_SMP("m_pi"), 2)),
                      List(0, Power(ME 2), 0,
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) *
          Power(
              Global_\
SMP("m_e"),
              2)) /
             s -
         4 * Power(Pi, 2) * u *
             Global_PaVe(1,
                         List(Power(ME 2),
                              Power(ME 2),
                              Power(Global_SMP("m_pi"), 2),
                              Power(
                                  Global_SMP("m_\
pi"),
                                  2),
                              s,
                              -s - t + 2 * Power(ME 2) +
                                  2 * Power(Global_SMP("m_pi"), 2)),
                         List(0, Power(ME 2), 0,
                              Power(Global_SMP("m_pi"), 2)),
                         Rule(Global_PaVeAutoOrder, True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) *
             Power(
                 Global_\
SMP("m_e"),
                 2) -
         (2 * Power(Pi, 2) * t * u *
          Global_PaVe(1,
                      List(Power(ME 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(
                               Global_\
SMP("m_pi"),
                               2),
                           s,
                           -s - t + 2 * Power(ME 2) +
                               2 * Power(Global_SMP("m_pi"), 2)),
                      List(0, Power(ME 2), 0,
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) *
          Power(
              Global_\
SMP("m_e"),
              2)) /
             s -
         2 * Power(Pi, 2) *
             Global_PaVe(1,
                         List(t, Power(ME 2), s,
                              Power(Global_SMP("m_pi"), 2),
                              Power(Global_SMP("m_pi"), 2),
                              Power(Global_SMP("m_\
e"),
                                    2)),
                         List(Power(Global_SMP("m_pi"), 2),
                              Power(ME 2), 0, 0),
                         Rule(Global_PaVeAutoOrder, True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 4) +
         (2 * Power(Pi, 2) * t *
          Global_PaVe(1,
                      List(t, Power(ME 2), s,
                           Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(
                               Global_\
SMP("m_e"),
                               2)),
                      List(Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_\
e"),
                                 2),
                           0, 0),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(
                          Global_\
PaVeAutoReduce,
                          True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 4)) /
             s -
         (2 * Power(Pi, 2) * u *
          Global_PaVe(1,
                      List(t, Power(ME 2), s,
                           Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(
                               Global_\
SMP("m_e"),
                               2)),
                      List(Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_\
e"),
                                 2),
                           0, 0),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(
                          Global_\
PaVeAutoReduce,
                          True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 4)) /
             s +
         2 * Power(Pi, 2) *
             Global_PaVe(1,
                         List(Power(ME 2),
                              Power(ME 2),
                              Power(Global_SMP("m_pi"), 2),
                              Power(
                                  Global_SMP("m_\
pi"),
                                  2),
                              s,
                              -s - t + 2 * Power(ME 2) +
                                  2 * Power(Global_SMP("m_pi"), 2)),
                         List(0, Power(ME 2), 0,
                              Power(Global_SMP("m_pi"), 2)),
                         Rule(Global_PaVeAutoOrder, True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) *
             Power(
                 Global_\
SMP("m_e"),
                 4) -
         (2 * Power(Pi, 2) * t *
          Global_PaVe(1,
                      List(Power(ME 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(
                               Global_SMP("m_\
pi"),
                               2),
                           s,
                           -s - t + 2 * Power(ME 2) +
                               2 * Power(Global_SMP("m_pi"), 2)),
                      List(0, Power(ME 2), 0,
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) *
          Power(
              Global_\
SMP("m_e"),
              4)) /
             s +
         (2 * Power(Pi, 2) * u *
          Global_PaVe(1,
                      List(Power(ME 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(
                               Global_SMP("m_\
pi"),
                               2),
                           s,
                           -s - t + 2 * Power(ME 2) +
                               2 * Power(Global_SMP("m_pi"), 2)),
                      List(0, Power(ME 2), 0,
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) *
          Power(
              Global_\
SMP("m_e"),
              4)) /
             s +
         (Power(Pi, 2) * s *
          Global_D0(s, Power(ME 2), t,
                    Power(Global_\
SMP("m_pi"),
                          2),
                    Power(ME 2), Power(Global_SMP("m_pi"), 2),
                    0, 0, Power(ME 2),
                    Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2)) /
             2. +
         2 * Power(Pi, 2) * t *
             Global_D0(s, Power(ME 2), t,
                       Power(Global_\
SMP("m_pi"),
                             2),
                       Power(ME 2),
                       Power(Global_SMP("m_pi"), 2), 0, 0,
                       Power(ME 2),
                       Power(Global_SMP("m_pi"), 2)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2) -
         (Power(Pi, 2) * Power(t, 2) *
          Global_D0(s, Power(ME 2), t,
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_\
SMP("m_pi"),
                          2),
                    0, 0, Power(ME 2),
                    Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2)) /
             (2. * s) +
         (Power(Pi, 2) * t * u *
          Global_D0(s, Power(ME 2), t,
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_SMP("m_\
pi"),
                          2),
                    0, 0, Power(ME 2),
                    Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2)) /
             s -
         (Power(Pi, 2) * Power(u, 2) *
          Global_D0(s, Power(ME 2), t,
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_\
SMP("m_pi"),
                          2),
                    0, 0, Power(ME 2),
                    Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2)) /
             (2. * s) +
         (5 * Power(Pi, 2) * s *
          Global_D0(s, Power(ME 2),
                    -s - t + 2 * Power(ME 2) +
                        2 * Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_SMP("m_pi"), 2), 0, 0,
                    Power(ME 2), Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2)) /
             2. +
         2 * Power(Pi, 2) * t *
             Global_D0(
                 s, Power(ME 2),
                 -s - t + 2 * Power(ME 2) +
                     2 * Power(Global_SMP("m_pi"), 2),
                 Power(Global_SMP("m_pi"), 2), Power(ME 2),
                 Power(Global_SMP("m_pi"), 2), 0, 0,
                 Power(ME 2), Power(Global_SMP("m_pi"), 2)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2) -
         (Power(Pi, 2) * Power(t, 2) *
          Global_D0(s, Power(ME 2),
                    -s - t + 2 * Power(ME 2) +
                        2 * Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_SMP("m_pi"), 2), 0, 0,
                    Power(ME 2), Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2)) /
             (2. * s) +
         (Power(Pi, 2) * t * u *
          Global_D0(s, Power(ME 2),
                    -s - t + 2 * Power(ME 2) +
                        2 * Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_SMP("m_pi"), 2), 0, 0,
                    Power(ME 2), Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2)) /
             s -
         (Power(Pi, 2) * Power(u, 2) *
          Global_D0(s, Power(ME 2),
                    -s - t + 2 * Power(ME 2) +
                        2 * Power(Global_SMP("m_pi"), 2),
                    Power(Global_SMP("m_pi"), 2), Power(ME 2),
                    Power(Global_SMP("m_pi"), 2), 0, 0,
                    Power(ME 2), Power(Global_SMP("m_pi"), 2)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2)) /
             (2. * s) -
         Power(Pi, 2) * s *
             Global_PaVe(1,
                         List(s, Power(ME 2), t,
                              Power(Global_SMP("m_pi"), 2),
                              Power(ME 2),
                              Power(
                                  Global_SMP("m_\
pi"),
                                  2)),
                         List(0, 0, Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         Rule(Global_PaVeAutoOrder, True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 2) -
         4 * Power(Pi, 2) * t *
             Global_PaVe(1,
                         List(s, Power(ME 2), t,
                              Power(Global_SMP("m_pi"), 2),
                              Power(ME 2),
                              Power(
                                  Global_\
SMP("m_pi"),
                                  2)),
                         List(0, 0, Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         Rule(Global_PaVeAutoOrder, True),
                         Rule(
                             Global_\
PaVeAutoReduce,
                             True)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) *
             Power(
                 Global_SMP("m_\
pi"),
                 2) +
         (Power(Pi, 2) * Power(t, 2) *
          Global_PaVe(1,
                      List(s, Power(ME 2), t,
                           Power(Global_SMP("m_pi"), 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      List(0, 0, Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(
                          Global_\
PaVeAutoReduce,
                          True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) *
          Power(
              Global_SMP("m_\
pi"),
              2)) /
             s -
         (2 * Power(Pi, 2) * t * u *
          Global_PaVe(1,
                      List(s, Power(ME 2), t,
                           Power(Global_SMP("m_pi"), 2),
                           Power(ME 2),
                           Power(
                               Global_\
SMP("m_pi"),
                               2)),
                      List(0, 0, Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(
                          Global_\
PaVeAutoReduce,
                          True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) *
          Power(
              Global_SMP("m_\
pi"),
              2)) /
             s +
         (Power(Pi, 2) * Power(u, 2) *
          Global_PaVe(1,
                      List(s, Power(ME 2), t,
                           Power(Global_SMP("m_pi"), 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      List(0, 0, Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(
                          Global_\
PaVeAutoReduce,
                          True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) *
          Power(
              Global_SMP("m_\
pi"),
              2)) /
             s -
         5 * Power(Pi, 2) * s *
             Global_PaVe(1,
                         List(s, Power(ME 2),
                              -s - t + 2 * Power(ME 2) +
                                  2 * Power(Global_SMP("m_pi"), 2),
                              Power(Global_SMP("m_pi"), 2),
                              Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         List(0, 0, Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         Rule(Global_\
PaVeAutoOrder,
                              True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Global_\
SMP("e"),
                   6) *
             Power(Global_SMP("m_pi"), 2) -
         4 * Power(Pi, 2) * t *
             Global_PaVe(1,
                         List(s, Power(ME 2),
                              -s - t + 2 * Power(ME 2) +
                                  2 * Power(Global_SMP("m_pi"), 2),
                              Power(Global_SMP("m_pi"), 2),
                              Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         List(0, 0, Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         Rule(Global_\
PaVeAutoOrder,
                              True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Global_\
SMP("e"),
                   6) *
             Power(Global_SMP("m_pi"), 2) +
         (Power(Pi, 2) * Power(t, 2) *
          Global_PaVe(1,
                      List(s, Power(ME 2),
                           -s - t + 2 * Power(ME 2) +
                               2 * Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      List(0, 0, Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_\
PaVeAutoOrder,
                           True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Global_\
SMP("e"),
                6) *
          Power(Global_SMP("m_pi"), 2)) /
             s -
         (2 * Power(Pi, 2) * t * u *
          Global_PaVe(1,
                      List(s, Power(ME 2),
                           -s - t + 2 * Power(ME 2) +
                               2 * Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      List(0, 0, Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_\
PaVeAutoOrder,
                           True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Global_\
SMP("e"),
                6) *
          Power(Global_SMP("m_pi"), 2)) /
             s +
         (Power(Pi, 2) * Power(u, 2) *
          Global_PaVe(1,
                      List(s, Power(ME 2),
                           -s - t + 2 * Power(ME 2) +
                               2 * Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      List(0, 0, Power(ME 2),
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_\
PaVeAutoOrder,
                           True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Global_\
SMP("e"),
                6) *
          Power(Global_SMP("m_pi"), 2)) /
             s -
         2 * Power(Pi, 2) *
             Global_D0(s, Power(ME 2), t,
                       Power(Global_\
SMP("m_pi"),
                             2),
                       Power(ME 2),
                       Power(Global_SMP("m_pi"), 2), 0, 0,
                       Power(ME 2),
                       Power(Global_SMP("m_pi"), 2)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
             Power(Global_SMP("m_pi"), 2) -
         2 * Power(Pi, 2) *
             Global_D0(
                 s, Power(ME 2),
                 -s - t + 2 * Power(ME 2) +
                     2 * Power(Global_SMP("m_pi"), 2),
                 Power(Global_SMP("m_pi"), 2), Power(ME 2),
                 Power(Global_SMP("m_pi"), 2), 0, 0,
                 Power(ME 2), Power(Global_SMP("m_pi"), 2)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
             Power(Global_SMP("m_pi"), 2) +
         4 * Power(Pi, 2) *
             Global_PaVe(1,
                         List(s, Power(ME 2), t,
                              Power(Global_SMP("m_pi"), 2),
                              Power(ME 2),
                              Power(
                                  Global_SMP("m_\
pi"),
                                  2)),
                         List(0, 0, Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         Rule(Global_PaVeAutoOrder, True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
             Power(Global_SMP("m_pi"), 2) +
         4 * Power(Pi, 2) *
             Global_PaVe(1,
                         List(s, Power(ME 2),
                              -s - t + 2 * Power(ME 2) +
                                  2 * Power(Global_SMP("m_pi"), 2),
                              Power(Global_SMP("m_pi"), 2),
                              Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         List(0, 0, Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         Rule(Global_\
PaVeAutoOrder,
                              True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Global_\
SMP("e"),
                   6) *
             Power(ME 2) * Power(Global_SMP("m_pi"), 2) -
         2 * Power(Pi, 2) *
             Global_PaVe(1,
                         List(t, Power(ME 2), s,
                              Power(Global_SMP("m_pi"), 2),
                              Power(Global_SMP("m_pi"), 2),
                              Power(Global_SMP("m_\
e"),
                                    2)),
                         List(Power(Global_SMP("m_pi"), 2),
                              Power(ME 2), 0, 0),
                         Rule(Global_PaVeAutoOrder, True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
             Power(Global_SMP("m_pi"), 2) -
         (6 * Power(Pi, 2) * t *
          Global_PaVe(1,
                      List(t, Power(ME 2), s,
                           Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(
                               Global_\
SMP("m_e"),
                               2)),
                      List(Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_\
e"),
                                 2),
                           0, 0),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(
                          Global_\
PaVeAutoReduce,
                          True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
          Power(Global_SMP("m_pi"), 2)) /
             s -
         (2 * Power(Pi, 2) * u *
          Global_PaVe(1,
                      List(t, Power(ME 2), s,
                           Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(
                               Global_\
SMP("m_e"),
                               2)),
                      List(Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_\
e"),
                                 2),
                           0, 0),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(
                          Global_\
PaVeAutoReduce,
                          True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
          Power(Global_SMP("m_pi"), 2)) /
             s +
         2 * Power(Pi, 2) *
             Global_PaVe(1,
                         List(Power(ME 2),
                              Power(ME 2),
                              Power(Global_SMP("m_pi"), 2),
                              Power(
                                  Global_SMP("m_\
pi"),
                                  2),
                              s,
                              -s - t + 2 * Power(ME 2) +
                                  2 * Power(Global_SMP("m_pi"), 2)),
                         List(0, Power(ME 2), 0,
                              Power(Global_SMP("m_pi"), 2)),
                         Rule(Global_PaVeAutoOrder, True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) *
             Power(
                 Global_\
SMP("m_e"),
                 2) *
             Power(Global_SMP("m_pi"), 2) -
         (2 * Power(Pi, 2) * t *
          Global_PaVe(1,
                      List(Power(ME 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(
                               Global_SMP("m_\
pi"),
                               2),
                           s,
                           -s - t + 2 * Power(ME 2) +
                               2 * Power(Global_SMP("m_pi"), 2)),
                      List(0, Power(ME 2), 0,
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) *
          Power(
              Global_\
SMP("m_e"),
              2) *
          Power(Global_SMP("m_pi"), 2)) /
             s +
         (10 * Power(Pi, 2) * u *
          Global_PaVe(1,
                      List(Power(ME 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(
                               Global_SMP("m_\
pi"),
                               2),
                           s,
                           -s - t + 2 * Power(ME 2) +
                               2 * Power(Global_SMP("m_pi"), 2)),
                      List(0, Power(ME 2), 0,
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) *
          Power(
              Global_\
SMP("m_e"),
              2) *
          Power(Global_SMP("m_pi"), 2)) /
             s +
         (8 * Power(Pi, 2) *
          Global_PaVe(1,
                      List(t, Power(ME 2), s,
                           Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(
                               Global_\
SMP("m_e"),
                               2)),
                      List(Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_\
e"),
                                 2),
                           0, 0),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(
                          Global_\
PaVeAutoReduce,
                          True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 4) *
          Power(Global_SMP("m_pi"), 2)) /
             s -
         (8 * Power(Pi, 2) *
          Global_PaVe(1,
                      List(Power(ME 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(
                               Global_SMP("m_\
pi"),
                               2),
                           s,
                           -s - t + 2 * Power(ME 2) +
                               2 * Power(Global_SMP("m_pi"), 2)),
                      List(0, Power(ME 2), 0,
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) *
          Power(
              Global_\
SMP("m_e"),
              4) *
          Power(Global_SMP("m_pi"), 2)) /
             s -
         2 * Power(Pi, 2) *
             Global_D0(s, Power(ME 2), t,
                       Power(Global_\
SMP("m_pi"),
                             2),
                       Power(ME 2),
                       Power(Global_SMP("m_pi"), 2), 0, 0,
                       Power(ME 2),
                       Power(Global_SMP("m_pi"), 2)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 4) -
         2 * Power(Pi, 2) *
             Global_D0(
                 s, Power(ME 2),
                 -s - t + 2 * Power(ME 2) +
                     2 * Power(Global_SMP("m_pi"), 2),
                 Power(Global_SMP("m_pi"), 2), Power(ME 2),
                 Power(Global_SMP("m_pi"), 2), 0, 0,
                 Power(ME 2), Power(Global_SMP("m_pi"), 2)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 4) +
         4 * Power(Pi, 2) *
             Global_PaVe(1,
                         List(s, Power(ME 2), t,
                              Power(Global_SMP("m_pi"), 2),
                              Power(ME 2),
                              Power(
                                  Global_SMP("m_\
pi"),
                                  2)),
                         List(0, 0, Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         Rule(Global_PaVeAutoOrder, True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(Global_SMP("m_pi"), 4) +
         4 * Power(Pi, 2) *
             Global_PaVe(1,
                         List(s, Power(ME 2),
                              -s - t + 2 * Power(ME 2) +
                                  2 * Power(Global_SMP("m_pi"), 2),
                              Power(Global_SMP("m_pi"), 2),
                              Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         List(0, 0, Power(ME 2),
                              Power(Global_SMP("m_pi"), 2)),
                         Rule(Global_\
PaVeAutoOrder,
                              True),
                         Rule(Global_PaVeAutoReduce, True)) *
             Power(Global_\
SMP("e"),
                   6) *
             Power(Global_SMP("m_pi"), 4) +
         (8 * Power(Pi, 2) *
          Global_PaVe(1,
                      List(t, Power(ME 2), s,
                           Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(
                               Global_\
SMP("m_e"),
                               2)),
                      List(Power(Global_SMP("m_pi"), 2),
                           Power(Global_SMP("m_\
e"),
                                 2),
                           0, 0),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(
                          Global_\
PaVeAutoReduce,
                          True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) * Power(ME 2) *
          Power(Global_SMP("m_pi"), 4)) /
             s -
         (8 * Power(Pi, 2) *
          Global_PaVe(1,
                      List(Power(ME 2),
                           Power(ME 2),
                           Power(Global_SMP("m_pi"), 2),
                           Power(
                               Global_SMP("m_\
pi"),
                               2),
                           s,
                           -s - t + 2 * Power(ME 2) +
                               2 * Power(Global_SMP("m_pi"), 2)),
                      List(0, Power(ME 2), 0,
                           Power(Global_SMP("m_pi"), 2)),
                      Rule(Global_PaVeAutoOrder, True),
                      Rule(Global_PaVeAutoReduce, True)) *
          Power(Sqrt(4.*M_PI*m_alpha), 6) *
          Power(
              Global_\
SMP("m_e"),
              2) *
          Power(Global_SMP("m_pi"), 4)) /
             s;
}


DECLARE_VIRTUALME2_GETTER(EXTRAXS::PionPionVirtual, "PionPionVirtual")
Virtual_ME2_Base*
ATOOLS::Getter<PHASIC::Virtual_ME2_Base, PHASIC::Process_Info,
               EXTRAXS::PionPionVirtual>::operator()(const Process_Info& pi)
    const
{
  if (pi.m_loopgenerator.find("Internal") != 0) return NULL;
  // if (pi.m_fi.m_nlotype==nlo_type::loop) {
  Flavour_Vector fl(pi.ExtractFlavours());
  if (fl.size() != 4) return NULL;
  if (fl[0] == Flavour(kf_e) && fl[1] == fl[0].Bar() &&
      (fl[2].Kfcode() == kf_pi_plus || fl[2].Kfcode() == -kf_pi_plus) &&
      fl[3] == fl[2].Bar()) {
    return new PionPionVirtual(pi, fl, 0.0, 0);
  }
  return NULL;
}
