#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "METOOLS/Loops/Master_Integrals.H"
#include "clooptools.h"
#include <iostream>
#include "ATOOLS/Phys/FormFactor.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace METOOLS;
using namespace MODEL;

namespace EXTRAXS {
  class PionPionVirtual : public PHASIC::Virtual_ME2_Base {
    double m_eps2, m_eps, m_fin, m_s, m_t;
    double ME2, MP2, m_alpha;
  public:
    PionPionVirtual(const Process_Info& pi, const Flavour_Vector& flavs,
                  const double& ep2, const double& ep);

    ~PionPionVirtual() {
    }

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
    Complex C0e(const double x, const double y, const double mass1, const double mass2);
    Complex C0e(const double x, const double mass);
    Complex D0e(const double z, const double x, const double y, const double mass1, const double mass2);
    inline double Den(const double &a, const double &b) {return (a==b?0.:1./(a-b));}
    std::unique_ptr<FormFactor> p_formfactor;

    // Complex MBub(double a, double b, double c);
    // Complex MTri(double a, double b, double c,
    //             double aa, double bb, double cc);
    // Complex MBox(double a, double b, double c);
    // inline Complex Conjugate(Complex a) {return conj(a);}
    // inline double Re(Complex a) {return a.real();}

  };
}

using namespace EXTRAXS;

PionPionVirtual::PionPionVirtual(const Process_Info& pi, const Flavour_Vector& flavs,
                  const double& ep2, const double& ep) :
      Virtual_ME2_Base(pi, flavs),
      m_eps2(ep2), m_eps(ep)
   {
      // if (!s_loader->LoadLibrary("ooptools")) THROW(fatal_error, "Failed to load libooptools.");
      ltini();
      Setlambda(0.);
      ME = Flavour(kf_e).Mass();
      MM = (Flavour(kf_mu).Mass());
      ML = (Flavour(kf_tau).Mass());
      MP = (Flavour(kf_pi).Mass());
      ME2 = (Flavour(kf_e).Mass()*Flavour(kf_e).Mass());
      MM2 = (Flavour(kf_mu).Mass()*Flavour(kf_mu).Mass());
      ML2 = (Flavour(kf_tau).Mass()*Flavour(kf_tau).Mass());
      MP2 = (Flavour(kf_pi).Mass()*Flavour(kf_pi).Mass());
      MZ = Flavour(kf_Z).Mass();
      MZ2 = MZ*MZ;
      double  MW= Flavour(kf_Wplus).Mass();
      double MW2 = MW*MW;

      MP2 = (Flavour(kf_pi).Mass()*Flavour(kf_pi).Mass());
      SW2 = std::abs(MODEL::s_model->ComplexConstant(string("csin2_thetaW")));
      CW2 = 1.-SW*SW;
      CW = sqrt(CW2);
      SW = sqrt(SW2);
      m_mode=2;
      m_IRscale=1.;
      m_UVscale=1.;
      p_formfactor = std::unique_ptr<FormFactor>(new FormFactor());
      Scoped_Settings s{ Settings::GetMainSettings()["YFS"] };
      m_photonmass = s["PHOTON_MASS"].Get<double>();
      int dimreg = s["Dim_Reg"].Get<int>();
      if(dimreg) m_photonmass = 0;
      setcmpbits(64);
      // setmudim(m_IRscale);
      // setdelta(4.*M_PI);
      // Setminmass(0.);
   }


void PionPionVirtual::Calc(const Vec4D_Vector& momenta) {
  double factor(1.);
  if      (m_stype&sbt::qcd) factor=2*M_PI/AlphaQCD();
  else if (m_stype&sbt::qed) factor=2*M_PI/AlphaQED();
  else THROW(fatal_error,"Unknown coupling.");
  // 1/epsIR
  // m_res.IR()=m_eps*factor;
  // // 1/epsIR2
  // m_res.IR2()=m_eps2*factor;
  // // finite
  // if(m_count==1000){
    // m_count=0;
  clearcache();// Since in general s will be different we gain nothing from cache
  m_s = (momenta[0]+momenta[1]).Abs2();
  m_alpha = (*aqed)(m_s);
  m_res.Finite()= 1./m_s/m_s*2*(BornBox()+BornTriangle()).real();
  // FORTRAN(ltexi)();
}
  
Complex PionPionVirtual::BornBox(){
return 
-16*m_alpha*M_PI*(ME2*(-2*MP2 + m_s + 2*m_t)*
   Conjugate(C0i(cc0, m_s, ME2, ME2, 0, 0, ME2)) + 
  ME2*(2*MP2 - m_s - 2*m_t)*Conjugate(C0i(cc2, m_s, ME2, ME2, 0, 0, ME2)) + 
  ME2*ME2*ME2*Conjugate(D0i(dd0, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, ME2, 
     MP2)) - 3*ME2*ME2*MP2*Conjugate(D0i(dd0, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) + 
  ME2*MP2*MP2*Conjugate(D0i(dd0, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) + MP2*MP2*MP2*Conjugate(D0i(dd0, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) + 
  2*ME2*ME2*m_s*Conjugate(D0i(dd0, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) - 3*ME2*MP2*m_s*Conjugate(D0i(dd0, m_s, ME2, 2*MP2 - m_s - m_t, MP2, 
     ME2, MP2, 0, 0, ME2, MP2)) - 
  MP2*MP2*m_s*Conjugate(D0i(dd0, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, ME2, 
     MP2)) + ME2*m_s*m_s*Conjugate(D0i(dd0, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)) + ME2*ME2*m_t*Conjugate(D0i(dd0, m_s, ME2, 2*MP2 - m_s - m_t, 
     MP2, ME2, MP2, 0, 0, ME2, MP2)) - 
  3*MP2*MP2*m_t*Conjugate(D0i(dd0, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) + 3*MP2*m_s*m_t*Conjugate(D0i(dd0, m_s, ME2, 2*MP2 - m_s - m_t, MP2, 
     ME2, MP2, 0, 0, ME2, MP2)) - 
  m_s*m_s*m_t*Conjugate(D0i(dd0, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, ME2, 
     MP2)) - ME2*m_t*m_t*Conjugate(D0i(dd0, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)) + 3*MP2*m_t*m_t*Conjugate(D0i(dd0, m_s, ME2, 2*MP2 - m_s - m_t, 
     MP2, ME2, MP2, 0, 0, ME2, MP2)) - 
  2*m_s*m_t*m_t*Conjugate(D0i(dd0, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, ME2, 
     MP2)) - m_t*m_t*m_t*Conjugate(D0i(dd0, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 
     0, ME2, MP2)) - ME2*ME2*ME2*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) + ME2*ME2*MP2*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 0, 
     0, ME2, MP2)) + 3*ME2*MP2*MP2*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)) + MP2*MP2*MP2*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 0, 
     0, ME2, MP2)) - ME2*ME2*m_s*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 0, 
     0, ME2, MP2)) - ME2*MP2*m_s*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)) + ME2*ME2*m_t*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)) - 4*ME2*MP2*m_t*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) - 
  3*MP2*MP2*m_t*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) + 
  2*ME2*m_s*m_t*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) + 
  MP2*m_s*m_t*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) + 
  ME2*m_t*m_t*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) + 
  3*MP2*m_t*m_t*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) - 
  m_s*m_t*m_t*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) - 
  m_t*m_t*m_t*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) - 
  ME2*ME2*ME2*Conjugate(D0i(dd2, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, ME2, 
     MP2)) + 3*ME2*ME2*MP2*Conjugate(D0i(dd2, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) - 
  ME2*MP2*MP2*Conjugate(D0i(dd2, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) + MP2*MP2*MP2*Conjugate(D0i(dd2, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) - 
  ME2*ME2*m_s*Conjugate(D0i(dd2, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, ME2, 
     MP2)) + ME2*MP2*m_s*Conjugate(D0i(dd2, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) - 
  MP2*MP2*m_s*Conjugate(D0i(dd2, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, ME2, 
     MP2)) - ME2*ME2*m_t*Conjugate(D0i(dd2, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)) + 2*ME2*MP2*m_t*Conjugate(D0i(dd2, m_s, ME2, 2*MP2 - m_s - m_t, 
     MP2, ME2, MP2, 0, 0, ME2, MP2)) - 
  3*MP2*MP2*m_t*Conjugate(D0i(dd2, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) - ME2*m_s*m_t*Conjugate(D0i(dd2, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) + 
  3*MP2*m_s*m_t*Conjugate(D0i(dd2, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) - m_s*m_s*m_t*Conjugate(D0i(dd2, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) - 
  ME2*m_t*m_t*Conjugate(D0i(dd2, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, ME2, 
     MP2)) + 3*MP2*m_t*m_t*Conjugate(D0i(dd2, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) - 
  2*m_s*m_t*m_t*Conjugate(D0i(dd2, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, ME2, 
     MP2)) - m_t*m_t*m_t*Conjugate(D0i(dd2, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 
     0, ME2, MP2)) + ME2*ME2*ME2*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) - ME2*ME2*MP2*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 0, 
     0, ME2, MP2)) + ME2*MP2*MP2*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)) + MP2*MP2*MP2*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 0, 
     0, ME2, MP2)) - ME2*MP2*m_s*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)) - ME2*ME2*m_t*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)) - 2*ME2*MP2*m_t*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) - 
  3*MP2*MP2*m_t*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) + 
  ME2*m_s*m_t*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) + 
  MP2*m_s*m_t*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) + 
  ME2*m_t*m_t*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) + 
  3*MP2*m_t*m_t*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) - 
  m_s*m_t*m_t*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) - 
  m_t*m_t*m_t*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) + 
  ME2*ME2*MP2*Conjugate(D0i(dd3, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) - 2*ME2*MP2*MP2*Conjugate(D0i(dd3, m_s, ME2, 2*MP2 - m_s - m_t, MP2, 
     ME2, MP2, 0, 0, ME2, MP2)) - 
  MP2*MP2*MP2*Conjugate(D0i(dd3, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, ME2, 
     MP2)) + ME2*MP2*m_s*Conjugate(D0i(dd3, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) + 
  2*MP2*MP2*m_t*Conjugate(D0i(dd3, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) - MP2*m_s*m_t*Conjugate(D0i(dd3, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) - 
  MP2*m_t*m_t*Conjugate(D0i(dd3, m_s, ME2, 2*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, ME2, 
     MP2)) - ME2*ME2*MP2*Conjugate(D0i(dd3, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) + 2*ME2*MP2*MP2*Conjugate(D0i(dd3, m_s, ME2, m_t, MP2, ME2, MP2, 0, 
     0, ME2, MP2)) + MP2*MP2*MP2*Conjugate(D0i(dd3, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) - ME2*MP2*m_s*Conjugate(D0i(dd3, m_s, ME2, m_t, MP2, ME2, MP2, 0, 
     0, ME2, MP2)) - 2*MP2*MP2*m_t*Conjugate(D0i(dd3, m_s, ME2, m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)) + MP2*m_s*m_t*Conjugate(D0i(dd3, m_s, ME2, m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)) + MP2*m_t*m_t*Conjugate(D0i(dd3, m_s, ME2, m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)))/m_s;
}


Complex PionPionVirtual::BornTriangle(){
return 8*m_alpha*M_PI/m_s/m_s*
 (-4*(ME2*ME2 - MP2*MP2 + ME2*(-2*MP2 + m_s) + 2*MP2*m_t - m_t*(m_s + m_t))*
   Conjugate(B0i(bb0, MP2, 0, MP2)) + 
  2*(ME2*ME2 - MP2*MP2 + ME2*(-2*MP2 + m_s) + 2*MP2*m_t - m_t*(m_s + m_t))*
   Conjugate(B0i(bb0, m_s, ME2, ME2)) - 
  2*ME2*ME2*Conjugate(B0i(bb1, MP2, 0, MP2)) + 
  4*ME2*MP2*Conjugate(B0i(bb1, MP2, 0, MP2)) + 
  2*MP2*MP2*Conjugate(B0i(bb1, MP2, 0, MP2)) - 
  2*ME2*m_s*Conjugate(B0i(bb1, MP2, 0, MP2)) - 
  4*MP2*m_t*Conjugate(B0i(bb1, MP2, 0, MP2)) + 
  2*m_s*m_t*Conjugate(B0i(bb1, MP2, 0, MP2)) + 
  2*m_t*m_t*Conjugate(B0i(bb1, MP2, 0, MP2)) + 
  4*ME2*ME2*ME2*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) - 
  8*ME2*ME2*MP2*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4*ME2*MP2*MP2*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2*ME2*ME2*m_s*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4*ME2*MP2*m_s*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2*MP2*MP2*m_s*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) - 
  2*ME2*m_s*m_s*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) + 
  8*ME2*MP2*m_t*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4*ME2*m_s*m_t*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4*MP2*m_s*m_t*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2*m_s*m_s*m_t*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4*ME2*m_t*m_t*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2*m_s*m_t*m_t*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4*ME2*ME2*MP2*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8*ME2*MP2*MP2*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) - 
  4*MP2*MP2*MP2*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) - 
  2*ME2*ME2*m_s*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) + 
  8*ME2*MP2*m_s*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) + 
  2*MP2*MP2*m_s*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) - 
  2*ME2*m_s*m_s*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) + 
  8*MP2*MP2*m_t*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8*MP2*m_s*m_t*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) + 
  2*m_s*m_s*m_t*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) - 
  4*MP2*m_t*m_t*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) + 
  2*m_s*m_t*m_t*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) - 
  4*ME2*ME2*Conjugate(C0i(cc00, ME2, m_s, ME2, 0, ME2, ME2)) + 
  8*ME2*MP2*Conjugate(C0i(cc00, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4*MP2*MP2*Conjugate(C0i(cc00, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4*ME2*m_s*Conjugate(C0i(cc00, ME2, m_s, ME2, 0, ME2, ME2)) - 
  8*MP2*m_t*Conjugate(C0i(cc00, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4*m_s*m_t*Conjugate(C0i(cc00, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4*m_t*m_t*Conjugate(C0i(cc00, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4*ME2*ME2*Conjugate(C0i(cc00, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8*ME2*MP2*Conjugate(C0i(cc00, MP2, m_s, MP2, 0, MP2, MP2)) - 
  4*MP2*MP2*Conjugate(C0i(cc00, MP2, m_s, MP2, 0, MP2, MP2)) + 
  4*ME2*m_s*Conjugate(C0i(cc00, MP2, m_s, MP2, 0, MP2, MP2)) + 
  8*MP2*m_t*Conjugate(C0i(cc00, MP2, m_s, MP2, 0, MP2, MP2)) - 
  4*m_s*m_t*Conjugate(C0i(cc00, MP2, m_s, MP2, 0, MP2, MP2)) - 
  4*m_t*m_t*Conjugate(C0i(cc00, MP2, m_s, MP2, 0, MP2, MP2)) + 
  8*ME2*ME2*ME2*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) - 
  16*ME2*ME2*MP2*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4*ME2*MP2*MP2*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) + 
  6*ME2*ME2*m_s*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2*MP2*MP2*m_s*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) - 
  ME2*m_s*m_s*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) + 
  8*ME2*MP2*m_t*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4*ME2*m_s*m_t*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4*MP2*m_s*m_t*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2*m_s*m_s*m_t*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4*ME2*m_t*m_t*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2*m_s*m_t*m_t*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) + 
  8*ME2*ME2*MP2*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) - 
  16*ME2*MP2*MP2*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8*MP2*MP2*MP2*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) - 
  3*ME2*ME2*m_s*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) + 
  14*ME2*MP2*m_s*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) + 
  3*MP2*MP2*m_s*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) - 
  3*ME2*m_s*m_s*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) + 
  16*MP2*MP2*m_t*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) - 
  14*MP2*m_s*m_t*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) + 
  3*m_s*m_s*m_t*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8*MP2*m_t*m_t*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) + 
  3*m_s*m_t*m_t*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) + 
  4*ME2*MP2*MP2*Conjugate(C0i(cc11, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4*ME2*MP2*m_s*Conjugate(C0i(cc11, ME2, m_s, ME2, 0, ME2, ME2)) + 
  ME2*m_s*m_s*Conjugate(C0i(cc11, ME2, m_s, ME2, 0, ME2, ME2)) - 
  8*ME2*MP2*m_t*Conjugate(C0i(cc11, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4*ME2*m_s*m_t*Conjugate(C0i(cc11, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4*ME2*m_t*m_t*Conjugate(C0i(cc11, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4*ME2*ME2*MP2*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8*ME2*MP2*MP2*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) - 
  4*MP2*MP2*MP2*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) - 
  ME2*ME2*m_s*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) + 
  6*ME2*MP2*m_s*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) + 
  MP2*MP2*m_s*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) - 
  ME2*m_s*m_s*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) + 
  8*MP2*MP2*m_t*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) - 
  6*MP2*m_s*m_t*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) + 
  m_s*m_s*m_t*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) - 
  4*MP2*m_t*m_t*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) + 
  m_s*m_t*m_t*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) + 
  8*ME2*MP2*MP2*Conjugate(C0i(cc12, ME2, m_s, ME2, 0, ME2, ME2)) - 
  8*ME2*MP2*m_s*Conjugate(C0i(cc12, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2*ME2*m_s*m_s*Conjugate(C0i(cc12, ME2, m_s, ME2, 0, ME2, ME2)) - 
  16*ME2*MP2*m_t*Conjugate(C0i(cc12, ME2, m_s, ME2, 0, ME2, ME2)) + 
  8*ME2*m_s*m_t*Conjugate(C0i(cc12, ME2, m_s, ME2, 0, ME2, ME2)) + 
  8*ME2*m_t*m_t*Conjugate(C0i(cc12, ME2, m_s, ME2, 0, ME2, ME2)) + 
  8*ME2*ME2*MP2*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) - 
  16*ME2*MP2*MP2*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8*MP2*MP2*MP2*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) - 
  2*ME2*ME2*m_s*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) + 
  12*ME2*MP2*m_s*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) + 
  2*MP2*MP2*m_s*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) - 
  2*ME2*m_s*m_s*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) + 
  16*MP2*MP2*m_t*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) - 
  12*MP2*m_s*m_t*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) + 
  2*m_s*m_s*m_t*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8*MP2*m_t*m_t*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) + 
  2*m_s*m_t*m_t*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) + 
  4*ME2*ME2*ME2*Conjugate(C0i(cc2, ME2, m_s, ME2, 0, ME2, ME2)) - 
  8*ME2*ME2*MP2*Conjugate(C0i(cc2, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2*ME2*ME2*m_s*Conjugate(C0i(cc2, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2*MP2*MP2*m_s*Conjugate(C0i(cc2, ME2, m_s, ME2, 0, ME2, ME2)) - 
  ME2*m_s*m_s*Conjugate(C0i(cc2, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4*MP2*m_s*m_t*Conjugate(C0i(cc2, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2*m_s*m_s*m_t*Conjugate(C0i(cc2, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2*m_s*m_t*m_t*Conjugate(C0i(cc2, ME2, m_s, ME2, 0, ME2, ME2)) + 
  8*ME2*ME2*MP2*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) - 
  16*ME2*MP2*MP2*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8*MP2*MP2*MP2*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) - 
  3*ME2*ME2*m_s*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) + 
  14*ME2*MP2*m_s*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) + 
  3*MP2*MP2*m_s*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) - 
  3*ME2*m_s*m_s*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) + 
  16*MP2*MP2*m_t*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) - 
  14*MP2*m_s*m_t*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) + 
  3*m_s*m_s*m_t*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8*MP2*m_t*m_t*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) + 
  3*m_s*m_t*m_t*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) + 
  4*ME2*MP2*MP2*Conjugate(C0i(cc22, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4*ME2*MP2*m_s*Conjugate(C0i(cc22, ME2, m_s, ME2, 0, ME2, ME2)) + 
  ME2*m_s*m_s*Conjugate(C0i(cc22, ME2, m_s, ME2, 0, ME2, ME2)) - 
  8*ME2*MP2*m_t*Conjugate(C0i(cc22, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4*ME2*m_s*m_t*Conjugate(C0i(cc22, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4*ME2*m_t*m_t*Conjugate(C0i(cc22, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4*ME2*ME2*MP2*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8*ME2*MP2*MP2*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) - 
  4*MP2*MP2*MP2*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) - 
  ME2*ME2*m_s*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) + 
  6*ME2*MP2*m_s*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) + 
  MP2*MP2*m_s*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) - 
  ME2*m_s*m_s*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) + 
  8*MP2*MP2*m_t*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) - 
  6*MP2*m_s*m_t*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) + 
  m_s*m_s*m_t*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) - 
  4*MP2*m_t*m_t*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) + 
  m_s*m_t*m_t*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) - 
  2*ME2*ME2*Re(B0i(bb0, ME2, 0, ME2)) + 4*ME2*MP2*Re(B0i(bb0, ME2, 0, ME2)) + 
  2*MP2*MP2*Re(B0i(bb0, ME2, 0, ME2)) - 2*ME2*m_s*Re(B0i(bb0, ME2, 0, ME2)) - 
  4*MP2*m_t*Re(B0i(bb0, ME2, 0, ME2)) + 2*m_s*m_t*Re(B0i(bb0, ME2, 0, ME2)) + 
  2*m_t*m_t*Re(B0i(bb0, ME2, 0, ME2)) + 3*ME2*ME2*Re(B0i(bb0, MP2, 0, MP2)) - 
  6*ME2*MP2*Re(B0i(bb0, MP2, 0, MP2)) - 3*MP2*MP2*Re(B0i(bb0, MP2, 0, MP2)) + 
  3*ME2*m_s*Re(B0i(bb0, MP2, 0, MP2)) + 6*MP2*m_t*Re(B0i(bb0, MP2, 0, MP2)) - 
  3*m_s*m_t*Re(B0i(bb0, MP2, 0, MP2)) - 3*m_t*m_t*Re(B0i(bb0, MP2, 0, MP2)) - 
  2*ME2*ME2*Re(B0i(bb1, ME2, 0, ME2)) + 4*ME2*MP2*Re(B0i(bb1, ME2, 0, ME2)) + 
  2*MP2*MP2*Re(B0i(bb1, ME2, 0, ME2)) - 2*ME2*m_s*Re(B0i(bb1, ME2, 0, ME2)) - 
  4*MP2*m_t*Re(B0i(bb1, ME2, 0, ME2)) + 2*m_s*m_t*Re(B0i(bb1, ME2, 0, ME2)) + 
  2*m_t*m_t*Re(B0i(bb1, ME2, 0, ME2)) + 2*ME2*ME2*Re(B0i(bb1, MP2, 0, MP2)) - 
  4*ME2*MP2*Re(B0i(bb1, MP2, 0, MP2)) - 2*MP2*MP2*Re(B0i(bb1, MP2, 0, MP2)) + 
  2*ME2*m_s*Re(B0i(bb1, MP2, 0, MP2)) + 4*MP2*m_t*Re(B0i(bb1, MP2, 0, MP2)) - 
  2*m_s*m_t*Re(B0i(bb1, MP2, 0, MP2)) - 2*m_t*m_t*Re(B0i(bb1, MP2, 0, MP2)) + 
  4*ME2*ME2*ME2*Re(B0i(dbb0, ME2, 0, ME2)) - 
  8*ME2*ME2*MP2*Re(B0i(dbb0, ME2, 0, ME2)) - 
  4*ME2*MP2*MP2*Re(B0i(dbb0, ME2, 0, ME2)) + 
  4*ME2*ME2*m_s*Re(B0i(dbb0, ME2, 0, ME2)) + 
  8*ME2*MP2*m_t*Re(B0i(dbb0, ME2, 0, ME2)) - 
  4*ME2*m_s*m_t*Re(B0i(dbb0, ME2, 0, ME2)) - 
  4*ME2*m_t*m_t*Re(B0i(dbb0, ME2, 0, ME2)) + 
  4*ME2*ME2*MP2*Re(B0i(dbb0, MP2, 0, MP2)) - 
  8*ME2*MP2*MP2*Re(B0i(dbb0, MP2, 0, MP2)) - 
  4*MP2*MP2*MP2*Re(B0i(dbb0, MP2, 0, MP2)) + 
  4*ME2*MP2*m_s*Re(B0i(dbb0, MP2, 0, MP2)) + 
  8*MP2*MP2*m_t*Re(B0i(dbb0, MP2, 0, MP2)) - 
  4*MP2*m_s*m_t*Re(B0i(dbb0, MP2, 0, MP2)) - 
  4*MP2*m_t*m_t*Re(B0i(dbb0, MP2, 0, MP2)) - 4*ME2*ME2*ME2*Re(B0i(dbb1, ME2, 0, ME2)) + 
  8*ME2*ME2*MP2*Re(B0i(dbb1, ME2, 0, ME2)) + 
  4*ME2*MP2*MP2*Re(B0i(dbb1, ME2, 0, ME2)) - 
  4*ME2*ME2*m_s*Re(B0i(dbb1, ME2, 0, ME2)) - 
  8*ME2*MP2*m_t*Re(B0i(dbb1, ME2, 0, ME2)) + 
  4*ME2*m_s*m_t*Re(B0i(dbb1, ME2, 0, ME2)) + 
  4*ME2*m_t*m_t*Re(B0i(dbb1, ME2, 0, ME2)) + 
  2*ME2*ME2*MP2*Re(B0i(dbb1, MP2, 0, MP2)) - 
  4*ME2*MP2*MP2*Re(B0i(dbb1, MP2, 0, MP2)) - 
  2*MP2*MP2*MP2*Re(B0i(dbb1, MP2, 0, MP2)) + 
  2*ME2*MP2*m_s*Re(B0i(dbb1, MP2, 0, MP2)) + 
  4*MP2*MP2*m_t*Re(B0i(dbb1, MP2, 0, MP2)) - 
  2*MP2*m_s*m_t*Re(B0i(dbb1, MP2, 0, MP2)) - 2*MP2*m_t*m_t*Re(B0i(dbb1, MP2, 0, MP2)));
}


DECLARE_VIRTUALME2_GETTER(EXTRAXS::PionPionVirtual,"PionPionVirtual")
Virtual_ME2_Base *ATOOLS::Getter
<PHASIC::Virtual_ME2_Base,PHASIC::Process_Info,EXTRAXS::PionPionVirtual>::
operator()(const Process_Info &pi) const
{
  if (pi.m_loopgenerator.find("Internal")!=0) return NULL;
  // if (pi.m_fi.m_nlotype==nlo_type::loop) {
  Flavour_Vector fl(pi.ExtractFlavours());
  if (fl.size()!=4) return NULL;
  if (fl[0]==Flavour(kf_e) && fl[1]==fl[0].Bar() &&
      (fl[2].Kfcode()==kf_pi_plus || fl[2].Kfcode()==-kf_pi_plus) && fl[3]==fl[2].Bar())
  {
    return new PionPionVirtual(pi,fl,0.0,0);
  }
  return NULL;
}
