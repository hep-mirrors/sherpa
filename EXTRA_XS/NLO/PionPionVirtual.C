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

using namespace PHASIC;
using namespace ATOOLS;
using namespace METOOLS;
using namespace MODEL;

namespace EXTRAXS {
  class PionPionVirtual : public PHASIC::Virtual_ME2_Base {
    double m_eps2, m_eps, m_fin, m_s, m_t;
    double ME2, MP2, m_alpha, ML2, MM2;
    double S, T, U;
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
    inline double Den(const double &a, const double &b) {return (a==b?0.:1./(a-b));}
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
      m_eps2(ep2), m_eps(ep), ME2(flavs[0].Mass()*flavs[0].Mass()),
      MP2(flavs[3].Mass()*flavs[3].Mass()){
     FORTRAN(ltini)();
     Setlambda(0.);
      MM2 = (Flavour(kf_mu).Mass()*Flavour(kf_mu).Mass());
      ML2 = (Flavour(kf_tau).Mass()*Flavour(kf_tau).Mass());
      m_mode=1;
      m_IRscale=1.;
      m_UVscale=1.;
      setmudim(m_IRscale);
      setdelta(m_IRscale);
      }


void PionPionVirtual::Calc(const Vec4D_Vector& momenta) {
  double factor(1.);
  if      (m_stype&sbt::qcd) factor=2.*M_PI/AlphaQCD();
  else if (m_stype&sbt::qed) factor=2.*M_PI/AlphaQED();
  else THROW(fatal_error,"Unknown coupling.");
  Setlambda(0.);
  m_s = (momenta[0]+momenta[1]).Abs2();
  S=m_s;
  T = (momenta[0]-momenta[2]).Abs2();
  m_alpha = (*aqed)(m_s);
  Complex ffull=BornTriangle()+BornBox()+BornSelf();
  m_res.Finite()= m_alpha*m_alpha*m_alpha*2.*(Full()).real();
  Setlambda(-1.);
  m_res.IR() = m_alpha*m_alpha*2.*Full().real();
  Setlambda(-2.);
  m_res.IR2() = m_alpha*m_alpha*2.*Full().real();
  clearcache();// Since in general s will be different we gain nothing from cache
}
  
Complex PionPionVirtual::Full(){
   // return BornTriangle();
   return 4.*M_PI*Den(S,0)*(ME2*(-T + U)*
     (4.*C0i(cc0,S,ME2,ME2,0,0,ME2) - 4.*C0i(cc2,S,ME2,ME2,0,0,ME2) + 
       4.*(-ME2 + T)*D0i(dd2,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2) + 
       4.*(-ME2 + U)*D0i(dd2,S,ME2,U,MP2,ME2,MP2,0,0,ME2,MP2) + 
       2.*(-T + U)*(C0i(cc1,ME2,S,ME2,0,ME2,ME2) + 
          C0i(cc11,ME2,S,ME2,0,ME2,ME2) + 2.*C0i(cc12,ME2,S,ME2,0,ME2,ME2) + 
          C0i(cc2,ME2,S,ME2,0,ME2,ME2) + C0i(cc22,ME2,S,ME2,0,ME2,ME2))*
        Den(S,0)) + (ME2*ME2 - MP2*MP2 + T*U - ME2*(T + U))*
     (4.*((ME2 + MP2 - T)*D0i(dd0,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2) + 
          (ME2 - MP2 + U)*D0i(dd2,S,ME2,U,MP2,ME2,MP2,0,0,ME2,MP2) + 
          MP2*D0i(dd3,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2)) - 
       4.*((ME2 + MP2 - U)*D0i(dd0,S,ME2,U,MP2,ME2,MP2,0,0,ME2,MP2) + 
          (ME2 - MP2 + T)*D0i(dd2,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2) + 
          MP2*D0i(dd3,S,ME2,U,MP2,ME2,MP2,0,0,ME2,MP2)) + 
       2.*(-1. - 4.*C0i(cc00,ME2,S,ME2,0,ME2,ME2) + 
          4.*(B0i(bb1,S,ME2,ME2) + B0i(bb1,S,ML2,ML2) + 
             B0i(bb1,S,MM2,MM2) + C0i(cc00,MP2,S,MP2,0,MP2,MP2)) + 
          (8.*MP2 - 3*S)*C0i(cc1,MP2,S,MP2,0,MP2,MP2) + 
          2.*(B0i(bb0,S,ME2,ME2) + 
             (2.*ME2 - S)*C0i(cc0,ME2,S,ME2,0,ME2,ME2) + 
             (2.*MP2 - S)*C0i(cc0,MP2,S,MP2,0,MP2,MP2) + 
             (4.*ME2 - S)*C0i(cc1,ME2,S,ME2,0,ME2,ME2) + 
             (4.*MP2 - S)*C0i(cc12,MP2,S,MP2,0,MP2,MP2)) - 
          2.*(B0i(bb0,MP2,0,MP2) + B0i(bb1,MP2,0,MP2) - 
             (2.*ME2 - 2.*MP2 + T + U)*C0i(cc2,ME2,S,ME2,0,ME2,ME2)) + 
          (8.*MP2 - 3*S)*C0i(cc2,MP2,S,MP2,0,MP2,MP2) + 
          (4.*MP2 - S)*(C0i(cc11,MP2,S,MP2,0,MP2,MP2) + 
             C0i(cc22,MP2,S,MP2,0,MP2,MP2)) + 
          4.*ME2*(-1. + 2.*B0i(bb0,ME2,0,ME2) - 2.*B0i(bb1,ME2,0,ME2))*
           Den(ME2,ME2) + 4.*MP2*(2.*B0i(bb0,MP2,0,MP2) + B0i(bb1,MP2,0,MP2))*
           Den(MP2,MP2))*Den(S,0) - 
       4.*(A0(MP2) + 4.*(B0i(bb00,S,ME2,ME2) + B0i(bb00,S,ML2,ML2) + 
             B0i(bb00,S,MM2,MM2)) - 
          2.*(A0(ME2) + A0(ML2) + A0(MM2) + B0i(bb00,S,MP2,MP2)))*
        pow(Den(S,0),2)));
}
//    return Den(S,0)*(2.*ME2*(ME2 + MP2 - T)*
//      (2.*Conjugate(C0i(cc0,S,ME2,ME2,0,0,ME2)) - 
//        2.*(Conjugate(C0i(cc2,S,ME2,ME2,0,0,ME2)) + 
//           (ME2 - 2.*MP2 + S + T)*
//            Conjugate(D0i(dd2,S,ME2,2.*MP2 - S - T,MP2,ME2,MP2,0,0,ME2,MP2)) 
// + (ME2 - T)*Conjugate(D0i(dd2,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2))) + 
//        (2.*MP2 - S - 2.*T)*(Conjugate(C0i(cc1,ME2,S,ME2,0,ME2,ME2)) + 
//           Conjugate(C0i(cc11,ME2,S,ME2,0,ME2,ME2)) + 
//           2.*Conjugate(C0i(cc12,ME2,S,ME2,0,ME2,ME2)) + 
//           Conjugate(C0i(cc2,ME2,S,ME2,0,ME2,ME2)) + 
//           Conjugate(C0i(cc22,ME2,S,ME2,0,ME2,ME2)))*Den(S,0)) - 
//     2.*ME2*(ME2 - MP2 + S + T)*(2.*Conjugate(C0i(cc0,S,ME2,ME2,0,0,ME2)) - 
//        2.*(Conjugate(C0i(cc2,S,ME2,ME2,0,0,ME2)) + 
//           (ME2 - 2.*MP2 + S + T)*
//            Conjugate(D0i(dd2,S,ME2,2.*MP2 - S - T,MP2,ME2,MP2,0,0,ME2,MP2)) 
// + (ME2 - T)*Conjugate(D0i(dd2,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2))) + 
//        (2.*MP2 - S - 2.*T)*(Conjugate(C0i(cc1,ME2,S,ME2,0,ME2,ME2)) + 
//           Conjugate(C0i(cc11,ME2,S,ME2,0,ME2,ME2)) + 
//           2.*Conjugate(C0i(cc12,ME2,S,ME2,0,ME2,ME2)) + 
//           Conjugate(C0i(cc2,ME2,S,ME2,0,ME2,ME2)) + 
//           Conjugate(C0i(cc22,ME2,S,ME2,0,ME2,ME2)))*Den(S,0)) + 
//     4.*ME2*MP2*(2.*((ME2 - MP2 + S + T)*
//            Conjugate(D0i(dd0,S,ME2,2.*MP2 - S - T,MP2,ME2,MP2,0,0,ME2,MP2)) 
// - (ME2 + MP2 - T)*Conjugate(D0i(dd0,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2)) + 
//           (-ME2 - MP2 + S + T)*
//            Conjugate(D0i(dd2,S,ME2,2.*MP2 - S - T,MP2,ME2,MP2,0,0,ME2,MP2)) 
// + (ME2 - MP2 + T)*Conjugate(D0i(dd2,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2)) + 
//           MP2*(Conjugate(D0i(dd3,S,ME2,2.*MP2 - S - T,MP2,ME2,MP2,0,0,ME2,
//                MP2)) - Conjugate(D0i(dd3,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2)))) 
// + (1. - 2.*Conjugate(B0i(bb0,S,ME2,ME2)) - 4.*Conjugate(B0i(bb1,S,ME2,ME2)) - 
//           4.*Conjugate(B0i(bb1,S,ML2,ML2)) - 
//           4.*Conjugate(B0i(bb1,S,MM2,MM2)) + 
//           2.*(-2.*ME2 + S)*Conjugate(C0i(cc0,ME2,S,ME2,0,ME2,ME2)) + 
//           4.*Conjugate(C0i(cc00,ME2,S,ME2,0,ME2,ME2)) - 
//           4.*Conjugate(C0i(cc00,MP2,S,MP2,0,MP2,MP2)) + 
//           S*(2.*Conjugate(C0i(cc0,MP2,S,MP2,0,MP2,MP2)) + 
//              2.*Conjugate(C0i(cc1,ME2,S,ME2,0,ME2,ME2)) + 
//              3.*Conjugate(C0i(cc1,MP2,S,MP2,0,MP2,MP2)) + 
//              Conjugate(C0i(cc11,MP2,S,MP2,0,MP2,MP2)) + 
//              2.*Conjugate(C0i(cc12,MP2,S,MP2,0,MP2,MP2)) + 
//              2.*Conjugate(C0i(cc2,ME2,S,ME2,0,ME2,ME2)) + 
//              3.*Conjugate(C0i(cc2,MP2,S,MP2,0,MP2,MP2)) + 
//              Conjugate(C0i(cc22,MP2,S,MP2,0,MP2,MP2))) - 
//           4.*MP2*(Conjugate(C0i(cc0,MP2,S,MP2,0,MP2,MP2)) + 
//              2.*Conjugate(C0i(cc1,MP2,S,MP2,0,MP2,MP2)) + 
//              Conjugate(C0i(cc11,MP2,S,MP2,0,MP2,MP2)) + 
//              2.*(Conjugate(C0i(cc12,MP2,S,MP2,0,MP2,MP2)) + 
//                 Conjugate(C0i(cc2,MP2,S,MP2,0,MP2,MP2))) + 
//              Conjugate(C0i(cc22,MP2,S,MP2,0,MP2,MP2))) - 
//           4.*ME2*(2.*Conjugate(C0i(cc1,ME2,S,ME2,0,ME2,ME2)) + 
//              Conjugate(C0i(cc2,ME2,S,ME2,0,ME2,ME2)) + 
//              (-1. + 2.*Conjugate(B0i(bb0,ME2,0,ME2)) - 
//                 2.*Conjugate(B0i(bb1,ME2,0,ME2)))*Den(ME2,ME2)) + 
//           Conjugate(B0i(bb0,MP2,0,MP2))*(2 - 8.*MP2*Den(MP2,MP2)) + 
//           2.*Conjugate(B0i(bb1,MP2,0,MP2))*(1. - 2.*MP2*Den(MP2,MP2)))*Den(S,0) 
// + 2.*(-2.*Conjugate(A0(ME2)) - 2.*Conjugate(A0(ML2)) - 2.*Conjugate(A0(MM2)) + 
//           Conjugate(A0(MP2)) + 
//           4.*(Conjugate(B0i(bb00,S,ME2,ME2)) + 
//              Conjugate(B0i(bb00,S,ML2,ML2)) + 
//              Conjugate(B0i(bb00,S,MM2,MM2))) - 
//           2.*Conjugate(B0i(bb00,S,MP2,MP2)))*pow(Den(S,0),2)) + 
//     (MP2*(MP2 - 2.*T) - (ME2 - T)*(ME2 + S + T))*
//      (2.*((ME2 - MP2 + S + T)*Conjugate(D0i(dd0,S,ME2,2.*MP2 - S - T,MP2,ME2,
//              MP2,0,0,ME2,MP2)) - 
//           (ME2 + MP2 - T)*Conjugate(D0i(dd0,S,ME2,T,MP2,ME2,MP2,0,0,ME2,
//              MP2)) + (-ME2 - MP2 + S + T)*
//            Conjugate(D0i(dd2,S,ME2,2.*MP2 - S - T,MP2,ME2,MP2,0,0,ME2,MP2)) 
// + (ME2 - MP2 + T)*Conjugate(D0i(dd2,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2)) + 
//           MP2*(Conjugate(D0i(dd3,S,ME2,2.*MP2 - S - T,MP2,ME2,MP2,0,0,ME2,
//                MP2)) - Conjugate(D0i(dd3,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2)))) 
// + (1. - 2.*Conjugate(B0i(bb0,S,ME2,ME2)) - 4.*Conjugate(B0i(bb1,S,ME2,ME2)) - 
//           4.*Conjugate(B0i(bb1,S,ML2,ML2)) - 
//           4.*Conjugate(B0i(bb1,S,MM2,MM2)) + 
//           2.*(-2.*ME2 + S)*Conjugate(C0i(cc0,ME2,S,ME2,0,ME2,ME2)) + 
//           4.*Conjugate(C0i(cc00,ME2,S,ME2,0,ME2,ME2)) - 
//           4.*Conjugate(C0i(cc00,MP2,S,MP2,0,MP2,MP2)) + 
//           S*(2.*Conjugate(C0i(cc0,MP2,S,MP2,0,MP2,MP2)) + 
//              2.*Conjugate(C0i(cc1,ME2,S,ME2,0,ME2,ME2)) + 
//              3.*Conjugate(C0i(cc1,MP2,S,MP2,0,MP2,MP2)) + 
//              Conjugate(C0i(cc11,MP2,S,MP2,0,MP2,MP2)) + 
//              2.*Conjugate(C0i(cc12,MP2,S,MP2,0,MP2,MP2)) + 
//              2.*Conjugate(C0i(cc2,ME2,S,ME2,0,ME2,ME2)) + 
//              3.*Conjugate(C0i(cc2,MP2,S,MP2,0,MP2,MP2)) + 
//              Conjugate(C0i(cc22,MP2,S,MP2,0,MP2,MP2))) - 
//           4.*MP2*(Conjugate(C0i(cc0,MP2,S,MP2,0,MP2,MP2)) + 
//              2.*Conjugate(C0i(cc1,MP2,S,MP2,0,MP2,MP2)) + 
//              Conjugate(C0i(cc11,MP2,S,MP2,0,MP2,MP2)) + 
//              2.*(Conjugate(C0i(cc12,MP2,S,MP2,0,MP2,MP2)) + 
//                 Conjugate(C0i(cc2,MP2,S,MP2,0,MP2,MP2))) + 
//              Conjugate(C0i(cc22,MP2,S,MP2,0,MP2,MP2))) - 
//           4.*ME2*(2.*Conjugate(C0i(cc1,ME2,S,ME2,0,ME2,ME2)) + 
//              Conjugate(C0i(cc2,ME2,S,ME2,0,ME2,ME2)) + 
//              (-1. + 2.*Conjugate(B0i(bb0,ME2,0,ME2)) - 
//                 2.*Conjugate(B0i(bb1,ME2,0,ME2)))*Den(ME2,ME2)) + 
//           Conjugate(B0i(bb0,MP2,0,MP2))*(2 - 8.*MP2*Den(MP2,MP2)) + 
//           2.*Conjugate(B0i(bb1,MP2,0,MP2))*(1. - 2.*MP2*Den(MP2,MP2)))*Den(S,0) 
// + 2.*(-2.*Conjugate(A0(ME2)) - 2.*Conjugate(A0(ML2)) - 2.*Conjugate(A0(MM2)) + 
//           Conjugate(A0(MP2)) + 
//           4.*(Conjugate(B0i(bb00,S,ME2,ME2)) + 
//              Conjugate(B0i(bb00,S,ML2,ML2)) + 
//              Conjugate(B0i(bb00,S,MM2,MM2))) - 
//           2.*Conjugate(B0i(bb00,S,MP2,MP2)))*pow(Den(S,0),2)) - 
//     (-(MP2*(MP2 - 2.*T)) + (ME2 - T)*(ME2 + S + T))*
//      (2.*((ME2 - MP2 + S + T)*Conjugate(D0i(dd0,S,ME2,2.*MP2 - S - T,MP2,ME2,
//              MP2,0,0,ME2,MP2)) - 
//           (ME2 + MP2 - T)*Conjugate(D0i(dd0,S,ME2,T,MP2,ME2,MP2,0,0,ME2,
//              MP2)) + (-ME2 - MP2 + S + T)*
//            Conjugate(D0i(dd2,S,ME2,2.*MP2 - S - T,MP2,ME2,MP2,0,0,ME2,MP2)) + 
//           (ME2 - MP2 + T)*Conjugate(D0i(dd2,S,ME2,T,MP2,ME2,MP2,0,0,ME2,
//              MP2)) + MP2*(Conjugate(D0i(dd3,S,ME2,2.*MP2 - S - T,MP2,ME2,MP2,
//                0,0,ME2,MP2)) - 
//              Conjugate(D0i(dd3,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2)))) + 
//        (1. - 2.*Conjugate(B0i(bb0,S,ME2,ME2)) - 
//           4.*Conjugate(B0i(bb1,S,ME2,ME2)) - 
//           4.*Conjugate(B0i(bb1,S,ML2,ML2)) - 
//           4.*Conjugate(B0i(bb1,S,MM2,MM2)) + 
//           2.*(-2.*ME2 + S)*Conjugate(C0i(cc0,ME2,S,ME2,0,ME2,ME2)) + 
//           4.*Conjugate(C0i(cc00,ME2,S,ME2,0,ME2,ME2)) - 
//           4.*Conjugate(C0i(cc00,MP2,S,MP2,0,MP2,MP2)) + 
//           S*(2.*Conjugate(C0i(cc0,MP2,S,MP2,0,MP2,MP2)) + 
//              2.*Conjugate(C0i(cc1,ME2,S,ME2,0,ME2,ME2)) + 
//              3.*Conjugate(C0i(cc1,MP2,S,MP2,0,MP2,MP2)) + 
//              Conjugate(C0i(cc11,MP2,S,MP2,0,MP2,MP2)) + 
//              2.*Conjugate(C0i(cc12,MP2,S,MP2,0,MP2,MP2)) + 
//              2.*Conjugate(C0i(cc2,ME2,S,ME2,0,ME2,ME2)) + 
//              3.*Conjugate(C0i(cc2,MP2,S,MP2,0,MP2,MP2)) + 
//              Conjugate(C0i(cc22,MP2,S,MP2,0,MP2,MP2))) - 
//           4.*MP2*(Conjugate(C0i(cc0,MP2,S,MP2,0,MP2,MP2)) + 
//              2.*Conjugate(C0i(cc1,MP2,S,MP2,0,MP2,MP2)) + 
//              Conjugate(C0i(cc11,MP2,S,MP2,0,MP2,MP2)) + 
//              2.*(Conjugate(C0i(cc12,MP2,S,MP2,0,MP2,MP2)) + 
//                 Conjugate(C0i(cc2,MP2,S,MP2,0,MP2,MP2))) + 
//              Conjugate(C0i(cc22,MP2,S,MP2,0,MP2,MP2))) - 
//           4.*ME2*(2.*Conjugate(C0i(cc1,ME2,S,ME2,0,ME2,ME2)) + 
//              Conjugate(C0i(cc2,ME2,S,ME2,0,ME2,ME2)) + 
//              (-1. + 2.*Conjugate(B0i(bb0,ME2,0,ME2)) - 
//                 2.*Conjugate(B0i(bb1,ME2,0,ME2)))*Den(ME2,ME2)) + 
//           Conjugate(B0i(bb0,MP2,0,MP2))*(2. - 8.*MP2*Den(MP2,MP2)) + 
//           2.*Conjugate(B0i(bb1,MP2,0,MP2))*(1. - 2.*MP2*Den(MP2,MP2)))*Den(S,0) 
// + 2.*(-2.*Conjugate(A0(ME2)) - 2.*Conjugate(A0(ML2)) - 2.*Conjugate(A0(MM2)) + 
//           Conjugate(A0(MP2)) + 
//           4.*(Conjugate(B0i(bb00,S,ME2,ME2)) + 
//              Conjugate(B0i(bb00,S,ML2,ML2)) + Conjugate(B0i(bb00,S,MM2,MM2))
// ) - 2.*Conjugate(B0i(bb00,S,MP2,MP2)))*pow(Den(S,0),2)));

// }


Complex PionPionVirtual::BornBox(){
return 
-16*M_PI*(ME2*(-2.*MP2 + m_s + 2.*m_t)*
   Conjugate(C0i(cc0, m_s, ME2, ME2, 0, 0, ME2)) + 
  ME2*(2.*MP2 - m_s - 2.*m_t)*Conjugate(C0i(cc2, m_s, ME2, ME2, 0, 0, ME2)) + 
  ME2*ME2*ME2*Conjugate(D0i(dd0, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, ME2, 
     MP2)) - 3.*ME2*ME2*MP2*Conjugate(D0i(dd0, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) + 
  ME2*MP2*MP2*Conjugate(D0i(dd0, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) + MP2*MP2*MP2*Conjugate(D0i(dd0, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) + 
  2.*ME2*ME2*m_s*Conjugate(D0i(dd0, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) - 3.*ME2*MP2*m_s*Conjugate(D0i(dd0, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, 
     ME2, MP2, 0, 0, ME2, MP2)) - 
  MP2*MP2*m_s*Conjugate(D0i(dd0, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, ME2, 
     MP2)) + ME2*m_s*m_s*Conjugate(D0i(dd0, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)) + ME2*ME2*m_t*Conjugate(D0i(dd0, m_s, ME2, 2.*MP2 - m_s - m_t, 
     MP2, ME2, MP2, 0, 0, ME2, MP2)) - 
  3.*MP2*MP2*m_t*Conjugate(D0i(dd0, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) + 3.*MP2*m_s*m_t*Conjugate(D0i(dd0, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, 
     ME2, MP2, 0, 0, ME2, MP2)) - 
  m_s*m_s*m_t*Conjugate(D0i(dd0, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, ME2, 
     MP2)) - ME2*m_t*m_t*Conjugate(D0i(dd0, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)) + 3.*MP2*m_t*m_t*Conjugate(D0i(dd0, m_s, ME2, 2.*MP2 - m_s - m_t, 
     MP2, ME2, MP2, 0, 0, ME2, MP2)) - 
  2.*m_s*m_t*m_t*Conjugate(D0i(dd0, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, ME2, 
     MP2)) - m_t*m_t*m_t*Conjugate(D0i(dd0, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 
     0, ME2, MP2)) - ME2*ME2*ME2*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) + ME2*ME2*MP2*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 0, 
     0, ME2, MP2)) + 3.*ME2*MP2*MP2*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)) + MP2*MP2*MP2*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 0, 
     0, ME2, MP2)) - ME2*ME2*m_s*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 0, 
     0, ME2, MP2)) - ME2*MP2*m_s*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)) + ME2*ME2*m_t*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)) - 4.*ME2*MP2*m_t*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) - 
  3.*MP2*MP2*m_t*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) + 
  2.*ME2*m_s*m_t*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) + 
  MP2*m_s*m_t*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) + 
  ME2*m_t*m_t*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) + 
  3.*MP2*m_t*m_t*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) - 
  m_s*m_t*m_t*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) - 
  m_t*m_t*m_t*Conjugate(D0i(dd0, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) - 
  ME2*ME2*ME2*Conjugate(D0i(dd2, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, ME2, 
     MP2)) + 3.*ME2*ME2*MP2*Conjugate(D0i(dd2, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) - 
  ME2*MP2*MP2*Conjugate(D0i(dd2, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) + MP2*MP2*MP2*Conjugate(D0i(dd2, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) - 
  ME2*ME2*m_s*Conjugate(D0i(dd2, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, ME2, 
     MP2)) + ME2*MP2*m_s*Conjugate(D0i(dd2, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) - 
  MP2*MP2*m_s*Conjugate(D0i(dd2, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, ME2, 
     MP2)) - ME2*ME2*m_t*Conjugate(D0i(dd2, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)) + 2.*ME2*MP2*m_t*Conjugate(D0i(dd2, m_s, ME2, 2.*MP2 - m_s - m_t, 
     MP2, ME2, MP2, 0, 0, ME2, MP2)) - 
  3.*MP2*MP2*m_t*Conjugate(D0i(dd2, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) - ME2*m_s*m_t*Conjugate(D0i(dd2, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) + 
  3.*MP2*m_s*m_t*Conjugate(D0i(dd2, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) - m_s*m_s*m_t*Conjugate(D0i(dd2, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) - 
  ME2*m_t*m_t*Conjugate(D0i(dd2, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, ME2, 
     MP2)) + 3.*MP2*m_t*m_t*Conjugate(D0i(dd2, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) - 
  2.*m_s*m_t*m_t*Conjugate(D0i(dd2, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, ME2, 
     MP2)) - m_t*m_t*m_t*Conjugate(D0i(dd2, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 
     0, ME2, MP2)) + ME2*ME2*ME2*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) - ME2*ME2*MP2*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 0, 
     0, ME2, MP2)) + ME2*MP2*MP2*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)) + MP2*MP2*MP2*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 0, 
     0, ME2, MP2)) - ME2*MP2*m_s*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)) - ME2*ME2*m_t*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)) - 2.*ME2*MP2*m_t*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) - 
  3.*MP2*MP2*m_t*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) + 
  ME2*m_s*m_t*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) + 
  MP2*m_s*m_t*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) + 
  ME2*m_t*m_t*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) + 
  3.*MP2*m_t*m_t*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) - 
  m_s*m_t*m_t*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) - 
  m_t*m_t*m_t*Conjugate(D0i(dd2, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, ME2, MP2)) + 
  ME2*ME2*MP2*Conjugate(D0i(dd3, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) - 2.*ME2*MP2*MP2*Conjugate(D0i(dd3, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, 
     ME2, MP2, 0, 0, ME2, MP2)) - 
  MP2*MP2*MP2*Conjugate(D0i(dd3, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, ME2, 
     MP2)) + ME2*MP2*m_s*Conjugate(D0i(dd3, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) + 
  2.*MP2*MP2*m_t*Conjugate(D0i(dd3, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) - MP2*m_s*m_t*Conjugate(D0i(dd3, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, 
     MP2, 0, 0, ME2, MP2)) - 
  MP2*m_t*m_t*Conjugate(D0i(dd3, m_s, ME2, 2.*MP2 - m_s - m_t, MP2, ME2, MP2, 0, 0, ME2, 
     MP2)) - ME2*ME2*MP2*Conjugate(D0i(dd3, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) + 2.*ME2*MP2*MP2*Conjugate(D0i(dd3, m_s, ME2, m_t, MP2, ME2, MP2, 0, 
     0, ME2, MP2)) + MP2*MP2*MP2*Conjugate(D0i(dd3, m_s, ME2, m_t, MP2, ME2, MP2, 0, 0, 
     ME2, MP2)) - ME2*MP2*m_s*Conjugate(D0i(dd3, m_s, ME2, m_t, MP2, ME2, MP2, 0, 
     0, ME2, MP2)) - 2.*MP2*MP2*m_t*Conjugate(D0i(dd3, m_s, ME2, m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)) + MP2*m_s*m_t*Conjugate(D0i(dd3, m_s, ME2, m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)) + MP2*m_t*m_t*Conjugate(D0i(dd3, m_s, ME2, m_t, MP2, ME2, MP2, 
     0, 0, ME2, MP2)))/m_s;
}

Complex PionPionVirtual::BornSelf(){
   return 4.*M_PI*(ME2*ME2 - MP2*MP2 + T*U - ME2*(T + U))*pow(Den(S,0),2)*
  (8.*(B0i(bb1,S,ME2,ME2) + B0i(bb1,S,ML2,ML2) + B0i(bb1,S,MM2,MM2)) - 
    4.*(A0(MP2) + 4.*(B0i(bb00,S,ME2,ME2) + B0i(bb00,S,ML2,ML2) + 
          B0i(bb00,S,MM2,MM2)) - 
       2.*(A0(ME2) + A0(ML2) + A0(MM2) + B0i(bb00,S,MP2,MP2)))*Den(S,0) - 
    8.*(Re(B0i(bb1,0,ME2,ME2)) + Re(B0i(bb1,0,ML2,ML2)) + 
       Re(B0i(bb1,0,MM2,MM2)) - 
       2.*(Re(B0i(dbb00,0,ME2,ME2)) + Re(B0i(dbb00,0,ML2,ML2)) + 
          Re(B0i(dbb00,0,MM2,MM2))) + Re(B0i(dbb00,0,MP2,MP2))));
}


Complex PionPionVirtual::BornTriangle(){
return 8.*M_PI/m_s/m_s*
 (-4.*(ME2*ME2 - MP2*MP2 + ME2*(-2.*MP2 + m_s) + 2.*MP2*m_t - m_t*(m_s + m_t))*
   Conjugate(B0i(bb0, MP2, 0, MP2)) + 
  2.*(ME2*ME2 - MP2*MP2 + ME2*(-2.*MP2 + m_s) + 2.*MP2*m_t - m_t*(m_s + m_t))*
   Conjugate(B0i(bb0, m_s, ME2, ME2)) - 
  2.*ME2*ME2*Conjugate(B0i(bb1, MP2, 0, MP2)) + 
  4.*ME2*MP2*Conjugate(B0i(bb1, MP2, 0, MP2)) + 
  2.*MP2*MP2*Conjugate(B0i(bb1, MP2, 0, MP2)) - 
  2.*ME2*m_s*Conjugate(B0i(bb1, MP2, 0, MP2)) - 
  4.*MP2*m_t*Conjugate(B0i(bb1, MP2, 0, MP2)) + 
  2.*m_s*m_t*Conjugate(B0i(bb1, MP2, 0, MP2)) + 
  2.*m_t*m_t*Conjugate(B0i(bb1, MP2, 0, MP2)) + 
  4.*ME2*ME2*ME2*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) - 
  8.*ME2*ME2*MP2*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4.*ME2*MP2*MP2*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2.*ME2*ME2*m_s*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4.*ME2*MP2*m_s*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2.*MP2*MP2*m_s*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) - 
  2.*ME2*m_s*m_s*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) + 
  8.*ME2*MP2*m_t*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4.*ME2*m_s*m_t*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4.*MP2*m_s*m_t*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2.*m_s*m_s*m_t*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4.*ME2*m_t*m_t*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2.*m_s*m_t*m_t*Conjugate(C0i(cc0, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4.*ME2*ME2*MP2*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8.*ME2*MP2*MP2*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) - 
  4.*MP2*MP2*MP2*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) - 
  2.*ME2*ME2*m_s*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) + 
  8.*ME2*MP2*m_s*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) + 
  2.*MP2*MP2*m_s*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) - 
  2.*ME2*m_s*m_s*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) + 
  8.*MP2*MP2*m_t*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8.*MP2*m_s*m_t*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) + 
  2.*m_s*m_s*m_t*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) - 
  4.*MP2*m_t*m_t*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) + 
  2.*m_s*m_t*m_t*Conjugate(C0i(cc0, MP2, m_s, MP2, 0, MP2, MP2)) - 
  4.*ME2*ME2*Conjugate(C0i(cc00, ME2, m_s, ME2, 0, ME2, ME2)) + 
  8.*ME2*MP2*Conjugate(C0i(cc00, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4.*MP2*MP2*Conjugate(C0i(cc00, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4.*ME2*m_s*Conjugate(C0i(cc00, ME2, m_s, ME2, 0, ME2, ME2)) - 
  8.*MP2*m_t*Conjugate(C0i(cc00, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4.*m_s*m_t*Conjugate(C0i(cc00, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4.*m_t*m_t*Conjugate(C0i(cc00, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4.*ME2*ME2*Conjugate(C0i(cc00, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8.*ME2*MP2*Conjugate(C0i(cc00, MP2, m_s, MP2, 0, MP2, MP2)) - 
  4.*MP2*MP2*Conjugate(C0i(cc00, MP2, m_s, MP2, 0, MP2, MP2)) + 
  4.*ME2*m_s*Conjugate(C0i(cc00, MP2, m_s, MP2, 0, MP2, MP2)) + 
  8.*MP2*m_t*Conjugate(C0i(cc00, MP2, m_s, MP2, 0, MP2, MP2)) - 
  4.*m_s*m_t*Conjugate(C0i(cc00, MP2, m_s, MP2, 0, MP2, MP2)) - 
  4.*m_t*m_t*Conjugate(C0i(cc00, MP2, m_s, MP2, 0, MP2, MP2)) + 
  8.*ME2*ME2*ME2*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) - 
  16*ME2*ME2*MP2*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4.*ME2*MP2*MP2*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) + 
  6*ME2*ME2*m_s*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2.*MP2*MP2*m_s*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) - 
  ME2*m_s*m_s*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) + 
  8.*ME2*MP2*m_t*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4.*ME2*m_s*m_t*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4.*MP2*m_s*m_t*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2.*m_s*m_s*m_t*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4.*ME2*m_t*m_t*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2.*m_s*m_t*m_t*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) + 
  8.*ME2*ME2*MP2*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) - 
  16*ME2*MP2*MP2*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8.*MP2*MP2*MP2*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) - 
  3.*ME2*ME2*m_s*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) + 
  14.*ME2*MP2*m_s*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) + 
  3.*MP2*MP2*m_s*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) - 
  3.*ME2*m_s*m_s*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) + 
  16*MP2*MP2*m_t*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) - 
  14.*MP2*m_s*m_t*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) + 
  3.*m_s*m_s*m_t*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8.*MP2*m_t*m_t*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) + 
  3.*m_s*m_t*m_t*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) + 
  4.*ME2*MP2*MP2*Conjugate(C0i(cc11, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4.*ME2*MP2*m_s*Conjugate(C0i(cc11, ME2, m_s, ME2, 0, ME2, ME2)) + 
  ME2*m_s*m_s*Conjugate(C0i(cc11, ME2, m_s, ME2, 0, ME2, ME2)) - 
  8.*ME2*MP2*m_t*Conjugate(C0i(cc11, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4.*ME2*m_s*m_t*Conjugate(C0i(cc11, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4.*ME2*m_t*m_t*Conjugate(C0i(cc11, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4.*ME2*ME2*MP2*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8.*ME2*MP2*MP2*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) - 
  4.*MP2*MP2*MP2*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) - 
  ME2*ME2*m_s*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) + 
  6*ME2*MP2*m_s*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) + 
  MP2*MP2*m_s*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) - 
  ME2*m_s*m_s*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) + 
  8.*MP2*MP2*m_t*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) - 
  6*MP2*m_s*m_t*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) + 
  m_s*m_s*m_t*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) - 
  4.*MP2*m_t*m_t*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) + 
  m_s*m_t*m_t*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) + 
  8.*ME2*MP2*MP2*Conjugate(C0i(cc12, ME2, m_s, ME2, 0, ME2, ME2)) - 
  8.*ME2*MP2*m_s*Conjugate(C0i(cc12, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2.*ME2*m_s*m_s*Conjugate(C0i(cc12, ME2, m_s, ME2, 0, ME2, ME2)) - 
  16*ME2*MP2*m_t*Conjugate(C0i(cc12, ME2, m_s, ME2, 0, ME2, ME2)) + 
  8.*ME2*m_s*m_t*Conjugate(C0i(cc12, ME2, m_s, ME2, 0, ME2, ME2)) + 
  8.*ME2*m_t*m_t*Conjugate(C0i(cc12, ME2, m_s, ME2, 0, ME2, ME2)) + 
  8.*ME2*ME2*MP2*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) - 
  16*ME2*MP2*MP2*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8.*MP2*MP2*MP2*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) - 
  2.*ME2*ME2*m_s*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) + 
  12.*ME2*MP2*m_s*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) + 
  2.*MP2*MP2*m_s*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) - 
  2.*ME2*m_s*m_s*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) + 
  16*MP2*MP2*m_t*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) - 
  12.*MP2*m_s*m_t*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) + 
  2.*m_s*m_s*m_t*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8.*MP2*m_t*m_t*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) + 
  2.*m_s*m_t*m_t*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) + 
  4.*ME2*ME2*ME2*Conjugate(C0i(cc2, ME2, m_s, ME2, 0, ME2, ME2)) - 
  8.*ME2*ME2*MP2*Conjugate(C0i(cc2, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2.*ME2*ME2*m_s*Conjugate(C0i(cc2, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2.*MP2*MP2*m_s*Conjugate(C0i(cc2, ME2, m_s, ME2, 0, ME2, ME2)) - 
  ME2*m_s*m_s*Conjugate(C0i(cc2, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4.*MP2*m_s*m_t*Conjugate(C0i(cc2, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2.*m_s*m_s*m_t*Conjugate(C0i(cc2, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2.*m_s*m_t*m_t*Conjugate(C0i(cc2, ME2, m_s, ME2, 0, ME2, ME2)) + 
  8.*ME2*ME2*MP2*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) - 
  16*ME2*MP2*MP2*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8.*MP2*MP2*MP2*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) - 
  3.*ME2*ME2*m_s*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) + 
  14.*ME2*MP2*m_s*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) + 
  3.*MP2*MP2*m_s*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) - 
  3.*ME2*m_s*m_s*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) + 
  16*MP2*MP2*m_t*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) - 
  14.*MP2*m_s*m_t*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) + 
  3.*m_s*m_s*m_t*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8.*MP2*m_t*m_t*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) + 
  3.*m_s*m_t*m_t*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) + 
  4.*ME2*MP2*MP2*Conjugate(C0i(cc22, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4.*ME2*MP2*m_s*Conjugate(C0i(cc22, ME2, m_s, ME2, 0, ME2, ME2)) + 
  ME2*m_s*m_s*Conjugate(C0i(cc22, ME2, m_s, ME2, 0, ME2, ME2)) - 
  8.*ME2*MP2*m_t*Conjugate(C0i(cc22, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4.*ME2*m_s*m_t*Conjugate(C0i(cc22, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4.*ME2*m_t*m_t*Conjugate(C0i(cc22, ME2, m_s, ME2, 0, ME2, ME2)) + 
  4.*ME2*ME2*MP2*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8.*ME2*MP2*MP2*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) - 
  4.*MP2*MP2*MP2*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) - 
  ME2*ME2*m_s*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) + 
  6*ME2*MP2*m_s*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) + 
  MP2*MP2*m_s*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) - 
  ME2*m_s*m_s*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) + 
  8.*MP2*MP2*m_t*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) - 
  6*MP2*m_s*m_t*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) + 
  m_s*m_s*m_t*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) - 
  4.*MP2*m_t*m_t*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) + 
  m_s*m_t*m_t*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) - 
  2.*ME2*ME2*Re(B0i(bb0, ME2, 0, ME2)) + 4.*ME2*MP2*Re(B0i(bb0, ME2, 0, ME2)) + 
  2.*MP2*MP2*Re(B0i(bb0, ME2, 0, ME2)) - 2.*ME2*m_s*Re(B0i(bb0, ME2, 0, ME2)) - 
  4.*MP2*m_t*Re(B0i(bb0, ME2, 0, ME2)) + 2.*m_s*m_t*Re(B0i(bb0, ME2, 0, ME2)) + 
  2.*m_t*m_t*Re(B0i(bb0, ME2, 0, ME2)) + 3.*ME2*ME2*Re(B0i(bb0, MP2, 0, MP2)) - 
  6*ME2*MP2*Re(B0i(bb0, MP2, 0, MP2)) - 3.*MP2*MP2*Re(B0i(bb0, MP2, 0, MP2)) + 
  3.*ME2*m_s*Re(B0i(bb0, MP2, 0, MP2)) + 6*MP2*m_t*Re(B0i(bb0, MP2, 0, MP2)) - 
  3.*m_s*m_t*Re(B0i(bb0, MP2, 0, MP2)) - 3.*m_t*m_t*Re(B0i(bb0, MP2, 0, MP2)) - 
  2.*ME2*ME2*Re(B0i(bb1, ME2, 0, ME2)) + 4.*ME2*MP2*Re(B0i(bb1, ME2, 0, ME2)) + 
  2.*MP2*MP2*Re(B0i(bb1, ME2, 0, ME2)) - 2.*ME2*m_s*Re(B0i(bb1, ME2, 0, ME2)) - 
  4.*MP2*m_t*Re(B0i(bb1, ME2, 0, ME2)) + 2.*m_s*m_t*Re(B0i(bb1, ME2, 0, ME2)) + 
  2.*m_t*m_t*Re(B0i(bb1, ME2, 0, ME2)) + 2.*ME2*ME2*Re(B0i(bb1, MP2, 0, MP2)) - 
  4.*ME2*MP2*Re(B0i(bb1, MP2, 0, MP2)) - 2.*MP2*MP2*Re(B0i(bb1, MP2, 0, MP2)) + 
  2.*ME2*m_s*Re(B0i(bb1, MP2, 0, MP2)) + 4.*MP2*m_t*Re(B0i(bb1, MP2, 0, MP2)) - 
  2.*m_s*m_t*Re(B0i(bb1, MP2, 0, MP2)) - 2.*m_t*m_t*Re(B0i(bb1, MP2, 0, MP2)) + 
  4.*ME2*ME2*ME2*Re(B0i(dbb0, ME2, 0, ME2)) - 
  8.*ME2*ME2*MP2*Re(B0i(dbb0, ME2, 0, ME2)) - 
  4.*ME2*MP2*MP2*Re(B0i(dbb0, ME2, 0, ME2)) + 
  4.*ME2*ME2*m_s*Re(B0i(dbb0, ME2, 0, ME2)) + 
  8.*ME2*MP2*m_t*Re(B0i(dbb0, ME2, 0, ME2)) - 
  4.*ME2*m_s*m_t*Re(B0i(dbb0, ME2, 0, ME2)) - 
  4.*ME2*m_t*m_t*Re(B0i(dbb0, ME2, 0, ME2)) + 
  4.*ME2*ME2*MP2*Re(B0i(dbb0, MP2, 0, MP2)) - 
  8.*ME2*MP2*MP2*Re(B0i(dbb0, MP2, 0, MP2)) - 
  4.*MP2*MP2*MP2*Re(B0i(dbb0, MP2, 0, MP2)) + 
  4.*ME2*MP2*m_s*Re(B0i(dbb0, MP2, 0, MP2)) + 
  8.*MP2*MP2*m_t*Re(B0i(dbb0, MP2, 0, MP2)) - 
  4.*MP2*m_s*m_t*Re(B0i(dbb0, MP2, 0, MP2)) - 
  4.*MP2*m_t*m_t*Re(B0i(dbb0, MP2, 0, MP2)) - 4.*ME2*ME2*ME2*Re(B0i(dbb1, ME2, 0, ME2)) + 
  8.*ME2*ME2*MP2*Re(B0i(dbb1, ME2, 0, ME2)) + 
  4.*ME2*MP2*MP2*Re(B0i(dbb1, ME2, 0, ME2)) - 
  4.*ME2*ME2*m_s*Re(B0i(dbb1, ME2, 0, ME2)) - 
  8.*ME2*MP2*m_t*Re(B0i(dbb1, ME2, 0, ME2)) + 
  4.*ME2*m_s*m_t*Re(B0i(dbb1, ME2, 0, ME2)) + 
  4.*ME2*m_t*m_t*Re(B0i(dbb1, ME2, 0, ME2)) + 
  2.*ME2*ME2*MP2*Re(B0i(dbb1, MP2, 0, MP2)) - 
  4.*ME2*MP2*MP2*Re(B0i(dbb1, MP2, 0, MP2)) - 
  2.*MP2*MP2*MP2*Re(B0i(dbb1, MP2, 0, MP2)) + 
  2.*ME2*MP2*m_s*Re(B0i(dbb1, MP2, 0, MP2)) + 
  4.*MP2*MP2*m_t*Re(B0i(dbb1, MP2, 0, MP2)) - 
  2.*MP2*m_s*m_t*Re(B0i(dbb1, MP2, 0, MP2)) - 2.*MP2*m_t*m_t*Re(B0i(dbb1, MP2, 0, MP2)));
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
  if ( (fl[0]==Flavour(kf_e) || fl[1]==Flavour(kf_e))  && fl[1]==fl[0].Bar() &&
      ((fl[2].Kfcode()==kf_pi_plus || fl[2].Kfcode()==-kf_pi_plus)) && fl[3]==fl[2].Bar())
  {
    return new PionPionVirtual(pi,fl,0.0,0);
  }
  return NULL;
}