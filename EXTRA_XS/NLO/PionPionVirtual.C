#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "METOOLS/Loops/Master_Integrals.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "clooptools.h"
#include <iostream>
#include "ATOOLS/Phys/FormFactor.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace METOOLS;
using namespace MODEL;
using namespace EXTRAXS;
using namespace std;

namespace EXTRAXS {
  class PionPionVirtual : public PHASIC::Virtual_ME2_Base {
    double m_eps2, m_eps, m_fin, m_s, m_t, m_u;
    double ME2, MP2, m_alpha, ML2, MM2, ME, MP, ML, MM, MZ, MZ2, MW, MW2;
    double S, T, U,CW, SW, SW2, CW2;
    double m_photonmass;
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
  if      (m_stype&sbt::qcd) factor=2.*M_PI/AlphaQCD();
  else if (m_stype&sbt::qed) factor=2.*M_PI/AlphaQED();
  else THROW(fatal_error,"Unknown coupling.");
  // for (int i = 0; i < 10; ++i)
  // {
     // m_photonmass = pow(10,i);
     // Setlambda(m_photonmass*m_photonmass);
     Setlambda(0);
     m_s = (momenta[2]+momenta[3]).Abs2();
     m_t = (momenta[0]-momenta[2]).Abs2();
     m_u = (momenta[0]-momenta[3]).Abs2();
     S=m_s;
     T = m_t;
     U = m_u;
     if(IsBad(S)||IsBad(T)||IsBad(U)){
      msg_Error()<<"Invalid phasespace point in "<<METHOD<<endl;
     }
     m_alpha = AlphaQED();
     m_res.Finite() =  (ISR()+FSR()+IFI()).real();
     // m_res.Finite() = (ISR()).real();  
     clearcache();
  // Complex ffull=BornTriangle()+BornBox();
  // // m_res.Finite()= 2.*(Full()).real()/m_alpha;
  Setlambda(-1.);
  PRINT_VAR(m_s);
  PRINT_VAR(ISR());
  // PRINT_VAR(ISR()+FSR());
  // PRINT_VAR(ISR()+IFI());
  // PRINT_VAR(FSR()+IFI());
  // PRINT_VAR(ISR()+FSR()+IFI());
  m_res.IR() =  -(ISR()).real()*m_alpha;
  // }
  // exit(1);
}
  
Complex PionPionVirtual::C0e(const double x, const double y, const double mass){
   return  Cget(mass,mass,m_s,x,mass,y);
}


Complex PionPionVirtual::C0e(const double x, const double mass){
   return  Cget(mass,mass,m_s,mass,x,mass);
}


Complex PionPionVirtual::C0e(const double x, const double y, const double mass1, const double mass2){
   return  Cget(mass1,mass2,x,mass1,y,mass2);
}

Complex PionPionVirtual::D0e(const double z, const double x, const double y, const double mass1, const double mass2){
   return  Dget(mass1,mass1,mass2,mass2,m_s,z,x,mass1,y,mass2);
}

double PionPionVirtual::Kappa(){
   double num = 4.*ME2*(12.*MP2*m_s-3.*m_s*m_s+2.*pow(m_t-m_u,2.));
   double den = (4.*ME2-m_s)*(4.*MP2*m_s-m_s*m_s+pow(m_t-m_u,2.));
   return 3.*m_s/(4.*ME2-m_s)-num/den; 
} 

double PionPionVirtual::Kappat(const double t, const double u){
   return (2.*MP2-t+u)*ME2+MP2*MP2+pow(MP2-t,2)-t*u;
}

Complex PionPionVirtual::ISR(){
//    return (32.*pow(m_alpha,3)*M_PI*(((ME2 + MP2 - T)*(T - U)*
//          (Conjugate(A0(ME2)) + 
//            ME2*(1. - 2.*Conjugate(B0i(bb0,ME2,0,ME2)) + 
//               Conjugate(B0i(bb0,S,ME2,ME2)))))/(4.*ME2 - S) - 
//       ((ME2 + MP2 - U)*(T - U)*
//          (Conjugate(A0(ME2)) + 
//            ME2*(1. - 2.*Conjugate(B0i(bb0,ME2,0,ME2)) + 
//               Conjugate(B0i(bb0,S,ME2,ME2)))))/(4.*ME2 - S) + 
//       2.*ME2*MP2*(2. - 4.*Conjugate(B0i(bb0,ME2,0,ME2)) + 
//          3.*Conjugate(B0i(bb0,S,ME2,ME2)) + 
//          2.*(-2.*ME2 + S)*Conjugate(C0i(cc0,ME2,S,ME2,0,ME2,ME2))) - 
//       (ME2*ME2 - MP2*MP2 + ME2*(2.*MP2 - T - U) + T*U)*
//        (2. - 4.*Conjugate(B0i(bb0,ME2,0,ME2)) + 
//          3.*Conjugate(B0i(bb0,S,ME2,ME2)) + 
//          2.*(-2.*ME2 + S)*Conjugate(C0i(cc0,ME2,S,ME2,0,ME2,ME2)))))/(S*S)
// ;
   // 2409.03469 Eq. 2.18
   return m_alpha/2/M_PI*(Kappa()*(Bget(m_s,ME2,ME2)-Bget(ME2,m_photonmass*m_photonmass,ME2))
                          +2.*(2.*ME2-m_s)*C0e(m_photonmass*m_photonmass,ME2)
                          +4.*ME2*DB0(ME2,m_photonmass*m_photonmass,ME2)
                         );

}

Complex PionPionVirtual::FSR(){
  //  return (16.*pow(m_alpha,3)*M_PI*(ME2*ME2 - MP2*MP2 + T*U - ME2*(T + U))*
  //   (3.*(-4.*MP2 + S)*Conjugate(A0(MP2)) + 
  //     4.*MP2*(4.*(3.*MP2 - S)*Conjugate(B0i(bb0,MP2,0,MP2)) + 
  //        (2.*MP2 - S)*(-2.*Conjugate(B0i(bb0,S,MP2,MP2)) + 
  //           (4.*MP2 - S)*Conjugate(C0i(cc0,MP2,S,MP2,0,MP2,MP2))))))/
  // (MP2*(4.*MP2 - S)*(S*S));
   // 2409.03469 Eq. 2.2
   return m_alpha/M_PI/(4.*MP2-m_s)*((2.*MP2-m_s)*( (4.*MP2-m_s)*C0e(m_photonmass*m_photonmass,MP2)
                                                    -2.*( B0i(bb0,m_s,MP2,MP2) - B0i(bb0,MP2,m_photonmass*m_photonmass,MP2) ) )
                                    +2.*MP2*(4.*MP2-m_s)*DB0(MP2,m_photonmass*m_photonmass,MP2)
                                    );
}

Complex PionPionVirtual::IFI(){
   // return (16.*pow(m_alpha,3)*M_PI*(-(ME2*(ME2 + MP2 - T)*
   //       (2.*Conjugate(C0i(cc0,S,ME2,ME2,0,0,ME2)) + 
   //         (-4.*Conjugate(B0i(bb0,ME2,0,ME2)) + 
   //            4.*Conjugate(B0i(bb0,S,0,0)) + 
   //            2.*S*Conjugate(C0i(cc0,S,ME2,ME2,0,0,ME2)))/(-4.*ME2 + S) + 
   //         (2.*((ME2 - T)*(4.*MP2 - S + T - U)*
   //               Conjugate(C0i(cc0,ME2,T,MP2,0,ME2,MP2)) + 
   //              (ME2 - U)*(4.*MP2 - S - T + U)*
   //               Conjugate(C0i(cc0,ME2,U,MP2,0,ME2,MP2)) + 
   //              pow(T - U,2)*Conjugate(C0i(cc0,S,ME2,ME2,0,0,ME2)) + 
   //              (8*(MP2*MP2) - 6.*MP2*S + S*S)*
   //               Conjugate(C0i(cc0,S,MP2,MP2,0,0,MP2)) + 
   //              S*(ME2 - T)*(ME2 - MP2 - T)*
   //               Conjugate(D0i(dd0,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2)) + 
   //              S*(ME2 - U)*(ME2 - MP2 - U)*
   //               Conjugate(D0i(dd0,S,ME2,U,MP2,ME2,MP2,0,0,ME2,MP2))))/
   //          (pow(ME2 - MP2,2) + (-2.*(ME2 + MP2) + S)*U + U*U))) + 
   //    ME2*(ME2 + MP2 - U)*(2.*Conjugate(C0i(cc0,S,ME2,ME2,0,0,ME2)) + 
   //       (-4.*Conjugate(B0i(bb0,ME2,0,ME2)) + 
   //          4.*Conjugate(B0i(bb0,S,0,0)) + 
   //          2.*S*Conjugate(C0i(cc0,S,ME2,ME2,0,0,ME2)))/(-4.*ME2 + S) + 
   //       (2.*((ME2 - T)*(4.*MP2 - S + T - U)*
   //             Conjugate(C0i(cc0,ME2,T,MP2,0,ME2,MP2)) + 
   //            (ME2 - U)*(4.*MP2 - S - T + U)*
   //             Conjugate(C0i(cc0,ME2,U,MP2,0,ME2,MP2)) + 
   //            pow(T - U,2)*Conjugate(C0i(cc0,S,ME2,ME2,0,0,ME2)) + 
   //            (8*(MP2*MP2) - 6.*MP2*S + S*S)*
   //             Conjugate(C0i(cc0,S,MP2,MP2,0,0,MP2)) + 
   //            S*(ME2 - T)*(ME2 - MP2 - T)*
   //             Conjugate(D0i(dd0,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2)) + 
   //            S*(ME2 - U)*(ME2 - MP2 - U)*
   //             Conjugate(D0i(dd0,S,ME2,U,MP2,ME2,MP2,0,0,ME2,MP2))))/
   //        (pow(ME2 - MP2,2) + (-2.*(ME2 + MP2) + S)*U + U*U)) - 
   //    2.*ME2*MP2*(4.*(ME2 + MP2 - T)*
   //        Conjugate(D0i(dd0,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2)) - 
   //       4.*(ME2 + MP2 - U)*Conjugate(D0i(dd0,S,ME2,U,MP2,ME2,MP2,0,0,ME2,
   //          MP2)) + (2.*(-2.*(ME2*ME2 + T*T + ME2*(MP2 - S - T - U) + 
   //               MP2*(-2.*MP2 + S + U))*
   //             Conjugate(C0i(cc0,ME2,T,MP2,0,ME2,MP2)) + 
   //            (-4.*(MP2*MP2) + (S + T)*(-3.*ME2 + S + T) + ME2*U + 
   //               2.*MP2*U + U*U)*Conjugate(C0i(cc0,ME2,U,MP2,0,ME2,MP2)) - 
   //            (S*T + 2.*ME2*(-2.*T + U) + U*(-2.*MP2 + T + U))*
   //             Conjugate(C0i(cc0,S,ME2,ME2,0,0,ME2)) + 
   //            (-2.*ME2*MP2 - 2.*(MP2*MP2) + MP2*(S + 3.*T - U) + S*(-T + U))*
   //             Conjugate(C0i(cc0,S,MP2,MP2,0,0,MP2)) - 
   //            S*(ME2 - T)*(ME2 - MP2 + T)*
   //             Conjugate(D0i(dd0,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2)) + 
   //            S*(ME2*ME2 + (-MP2 + S + T)*U - ME2*(MP2 + 2.*U))*
   //             Conjugate(D0i(dd0,S,ME2,U,MP2,ME2,MP2,0,0,ME2,MP2))))/
   //        (pow(ME2 - MP2,2) + (-2.*(ME2 + MP2) + S)*U + U*U)) + 
   //    (ME2*ME2 - MP2*MP2 + ME2*(2.*MP2 - T - U) + T*U)*
   //     (4.*(ME2 + MP2 - T)*Conjugate(D0i(dd0,S,ME2,T,MP2,ME2,MP2,0,0,ME2,
   //          MP2)) - 4.*(ME2 + MP2 - U)*
   //        Conjugate(D0i(dd0,S,ME2,U,MP2,ME2,MP2,0,0,ME2,MP2)) + 
   //       (2.*(-2.*(ME2*ME2 + T*T + ME2*(MP2 - S - T - U) + 
   //               MP2*(-2.*MP2 + S + U))*
   //             Conjugate(C0i(cc0,ME2,T,MP2,0,ME2,MP2)) + 
   //            (-4.*(MP2*MP2) + (S + T)*(-3.*ME2 + S + T) + ME2*U + 
   //               2.*MP2*U + U*U)*Conjugate(C0i(cc0,ME2,U,MP2,0,ME2,MP2)) - 
   //            (S*T + 2.*ME2*(-2.*T + U) + U*(-2.*MP2 + T + U))*
   //             Conjugate(C0i(cc0,S,ME2,ME2,0,0,ME2)) + 
   //            (-2.*ME2*MP2 - 2.*(MP2*MP2) + MP2*(S + 3.*T - U) + S*(-T + U))*
   //             Conjugate(C0i(cc0,S,MP2,MP2,0,0,MP2)) - 
   //            S*(ME2 - T)*(ME2 - MP2 + T)*
   //             Conjugate(D0i(dd0,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2)) + 
   //            S*(ME2*ME2 + (-MP2 + S + T)*U - ME2*(MP2 + 2.*U))*
   //             Conjugate(D0i(dd0,S,ME2,U,MP2,ME2,MP2,0,0,ME2,MP2))))/
   //        (pow(ME2 - MP2,2) + (-2.*(ME2 + MP2) + S)*U + U*U))))/S;

   // 2409.03469 Eq. 2.21
   Complex t1 = 4.*ME2*(m_t-m_u)/(4.*ME2-m_s)*(B0i(bb0,ME2,m_photonmass*m_photonmass,ME2)-B0i(bb1,m_s,m_photonmass*m_photonmass,m_photonmass*m_photonmass));
   
   t1 += (8.*ME2*(ME2-m_s)+m_s*m_s)/(4.*ME2-m_s)*(m_t-m_u)*C0e(m_photonmass*m_photonmass,m_photonmass*m_photonmass,ME2);
   
   t1 += (2.*MP2-m_s)*(m_t-m_u)*C0e(m_photonmass*m_photonmass,m_photonmass*m_photonmass,MP2)-2.*(ME2-m_t)*(ME2+MP2-m_t)*C0e(m_t, m_photonmass*m_photonmass, ME2, MP2);

   t1 += 2.*(ME2-m_u)*(ME2+MP2-m_u)*C0e(m_u, m_photonmass*m_photonmass, ME2, MP2);

   t1 += Kappat(m_t, m_u)*(ME2+MP2-m_t)*D0e(m_t,m_photonmass*m_photonmass,m_photonmass*m_photonmass, ME2, MP2);

   t1 += -Kappat(m_u,m_t)*(ME2+MP2-m_u)*D0e(m_u,m_photonmass*m_photonmass,m_photonmass*m_photonmass, ME2, MP2);

   return m_alpha/2/M_PI*(4.*m_s/(4.*MP2*m_s-m_s*m_s+pow(m_t-m_u,2)));
}


Complex PionPionVirtual::Full(){
   // return BornTriangle();
   return (-4.*m_alpha*M_PI*(-2.*m_alpha*m_alpha*ME2*S*pow(T - U,2)*
       (C0i(cc1,ME2,S,ME2,0,ME2,ME2) + C0i(cc11,ME2,S,ME2,0,ME2,ME2) + 
         2.*C0i(cc12,ME2,S,ME2,0,ME2,ME2) + C0i(cc2,ME2,S,ME2,0,ME2,ME2) + 
         C0i(cc22,ME2,S,ME2,0,ME2,ME2)) + 
      4.*m_alpha*m_alpha*ME2*(S*S)*(-T + U)*
       ((ME2 - T)*D0i(dd2,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2) + 
         (ME2 - U)*D0i(dd2,S,ME2,U,MP2,ME2,MP2,0,0,ME2,MP2)) + 
      2.*m_alpha*m_alpha*S*(MP2*MP2 - (ME2 - T)*(ME2 - U))*
       (-1. + 2.*B0i(bb0,S,ME2,ME2) - 
         2.*(-2.*ME2 + S)*C0i(cc0,ME2,S,ME2,0,ME2,ME2) - 
         2.*(-2.*MP2 + S)*C0i(cc0,MP2,S,MP2,0,MP2,MP2) - 
         4.*C0i(cc00,ME2,S,ME2,0,ME2,ME2) + 
         4.*C0i(cc00,MP2,S,MP2,0,MP2,MP2) - 
         2.*(-4.*ME2 + S)*C0i(cc1,ME2,S,ME2,0,ME2,ME2) + 
         (8.*MP2 - 3.*S)*C0i(cc1,MP2,S,MP2,0,MP2,MP2) - 
         2.*(-4.*MP2 + S)*C0i(cc12,MP2,S,MP2,0,MP2,MP2) + 
         2.*(2.*ME2 - 2.*MP2 + T + U)*C0i(cc2,ME2,S,ME2,0,ME2,ME2) + 
         (8.*MP2 - 3.*S)*C0i(cc2,MP2,S,MP2,0,MP2,MP2) + 
         (4.*MP2 - S)*(C0i(cc11,MP2,S,MP2,0,MP2,MP2) + 
            C0i(cc22,MP2,S,MP2,0,MP2,MP2)) - 
         4.*(Re(B0i(bb1,0,ME2,ME2)) + Re(B0i(bb1,0,ML2,ML2)) + 
            Re(B0i(bb1,0,MM2,MM2)) - 
            2.*(Re(B0i(dbb00,0,ME2,ME2)) + Re(B0i(dbb00,0,ML2,ML2)) + 
               Re(B0i(dbb00,0,MM2,MM2))) + Re(B0i(dbb00,0,MP2,MP2)))) + 
      S*(MP2*MP2 - (ME2 - T)*(ME2 - U))*
       (m_alpha*m_alpha*S*(4.*((ME2 + MP2 - T)*
                D0i(dd0,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2) + 
               (ME2 - MP2 + U)*D0i(dd2,S,ME2,U,MP2,ME2,MP2,0,0,ME2,MP2) + 
               MP2*D0i(dd3,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2)) - 
            4.*((ME2 + MP2 - U)*D0i(dd0,S,ME2,U,MP2,ME2,MP2,0,0,ME2,MP2) + 
               (ME2 - MP2 + T)*D0i(dd2,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2) + 
               MP2*D0i(dd3,S,ME2,U,MP2,ME2,MP2,0,0,ME2,MP2))) - 
         8.*m_alpha*m_alpha*(Re(B0i(bb1,0,ME2,ME2)) + Re(B0i(bb1,0,ML2,ML2)) + 
            Re(B0i(bb1,0,MM2,MM2)) - 
            2.*(Re(B0i(dbb00,0,ME2,ME2)) + Re(B0i(dbb00,0,ML2,ML2)) + 
               Re(B0i(dbb00,0,MM2,MM2))) + Re(B0i(dbb00,0,MP2,MP2)))) - 
      (MP2*MP2 - (ME2 - T)*(ME2 - U))*
       (16.*m_alpha*m_alpha*(B0i(bb00,S,ME2,ME2) + B0i(bb00,S,ML2,ML2) + 
            B0i(bb00,S,MM2,MM2)) - 
         8.*m_alpha*m_alpha*(A0(ME2) + A0(ML2) + A0(MM2) + B0i(bb00,S,MP2,MP2)) - 
         8.*m_alpha*m_alpha*S*(B0i(bb1,S,ME2,ME2) + B0i(bb1,S,ML2,ML2) + 
            B0i(bb1,S,MM2,MM2)) + 
         8.*m_alpha*m_alpha*S*(Re(B0i(bb1,0,ME2,ME2)) + Re(B0i(bb1,0,ML2,ML2)) + 
            Re(B0i(bb1,0,MM2,MM2)) - 
            2.*(Re(B0i(dbb00,0,ME2,ME2)) + Re(B0i(dbb00,0,ML2,ML2)) + 
               Re(B0i(dbb00,0,MM2,MM2))) + Re(B0i(dbb00,0,MP2,MP2))))))/
  pow(S,3);

   // return 4.*M_PI*Den(S,0)*(ME2*(-T + U)*
   //   (4.*C0i(cc0,S,ME2,ME2,0,0,ME2) - 4.*C0i(cc2,S,ME2,ME2,0,0,ME2) + 
   //     4.*(-ME2 + T)*D0i(dd2,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2) + 
   //     4.*(-ME2 + U)*D0i(dd2,S,ME2,U,MP2,ME2,MP2,0,0,ME2,MP2) + 
   //     2.*(-T + U)*(C0i(cc1,ME2,S,ME2,0,ME2,ME2) + 
   //        C0i(cc11,ME2,S,ME2,0,ME2,ME2) + 2.*C0i(cc12,ME2,S,ME2,0,ME2,ME2) + 
   //        C0i(cc2,ME2,S,ME2,0,ME2,ME2) + C0i(cc22,ME2,S,ME2,0,ME2,ME2))*
   //      Den(S,0)) + (ME2*ME2 - MP2*MP2 + T*U - ME2*(T + U))*
   //   (4.*((ME2 + MP2 - T)*D0i(dd0,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2) + 
   //        (ME2 - MP2 + U)*D0i(dd2,S,ME2,U,MP2,ME2,MP2,0,0,ME2,MP2) + 
   //        MP2*D0i(dd3,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2)) - 
   //     4.*((ME2 + MP2 - U)*D0i(dd0,S,ME2,U,MP2,ME2,MP2,0,0,ME2,MP2) + 
   //        (ME2 - MP2 + T)*D0i(dd2,S,ME2,T,MP2,ME2,MP2,0,0,ME2,MP2) + 
   //        MP2*D0i(dd3,S,ME2,U,MP2,ME2,MP2,0,0,ME2,MP2)) + 
   //     2.*(-1. - 4.*C0i(cc00,ME2,S,ME2,0,ME2,ME2) + 
   //        4.*(B0i(bb1,S,ME2,ME2) + B0i(bb1,S,ML2,ML2) + 
   //           B0i(bb1,S,MM2,MM2) + C0i(cc00,MP2,S,MP2,0,MP2,MP2)) + 
   //        (8.*MP2 - 3.*S)*C0i(cc1,MP2,S,MP2,0,MP2,MP2) + 
   //        2.*(B0i(bb0,S,ME2,ME2) + 
   //           (2.*ME2 - S)*C0i(cc0,ME2,S,ME2,0,ME2,ME2) + 
   //           (2.*MP2 - S)*C0i(cc0,MP2,S,MP2,0,MP2,MP2) + 
   //           (4.*ME2 - S)*C0i(cc1,ME2,S,ME2,0,ME2,ME2) + 
   //           (4.*MP2 - S)*C0i(cc12,MP2,S,MP2,0,MP2,MP2)) - 
   //        2.*(B0i(bb0,MP2,0,MP2) + B0i(bb1,MP2,0,MP2) - 
   //           (2.*ME2 - 2.*MP2 + T + U)*C0i(cc2,ME2,S,ME2,0,ME2,ME2)) + 
   //        (8.*MP2 - 3.*S)*C0i(cc2,MP2,S,MP2,0,MP2,MP2) + 
   //        (4.*MP2 - S)*(C0i(cc11,MP2,S,MP2,0,MP2,MP2) + 
   //           C0i(cc22,MP2,S,MP2,0,MP2,MP2)) + 
   //        4.*ME2*(-1. + 2.*B0i(bb0,ME2,0,ME2) - 2.*B0i(bb1,ME2,0,ME2))*
   //         Den(ME2,ME2) + 4.*MP2*(2.*B0i(bb0,MP2,0,MP2) + B0i(bb1,MP2,0,MP2))*
   //         Den(MP2,MP2))*Den(S,0) - 
   //     4.*(A0(MP2) + 4.*(B0i(bb00,S,ME2,ME2) + B0i(bb00,S,ML2,ML2) + 
   //           B0i(bb00,S,MM2,MM2)) - 
   //        2.*(A0(ME2) + A0(ML2) + A0(MM2) + B0i(bb00,S,MP2,MP2)))*
   //      pow(Den(S,0),2)));
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
//           Conjugate(B0i(bb0,MP2,0,MP2))*(2. - 8.*MP2*Den(MP2,MP2)) + 
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
//           Conjugate(B0i(bb0,MP2,0,MP2))*(2. - 8.*MP2*Den(MP2,MP2)) + 
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
-16.*M_PI*(ME2*(-2.*MP2 + m_s + 2.*m_t)*
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
return 8.*M_PI*
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
  16.*ME2*ME2*MP2*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4.*ME2*MP2*MP2*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) + 
  6.*ME2*ME2*m_s*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2.*MP2*MP2*m_s*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) - 
  ME2*m_s*m_s*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) + 
  8.*ME2*MP2*m_t*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4.*ME2*m_s*m_t*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4.*MP2*m_s*m_t*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2.*m_s*m_s*m_t*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) - 
  4.*ME2*m_t*m_t*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2.*m_s*m_t*m_t*Conjugate(C0i(cc1, ME2, m_s, ME2, 0, ME2, ME2)) + 
  8.*ME2*ME2*MP2*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) - 
  16.*ME2*MP2*MP2*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8.*MP2*MP2*MP2*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) - 
  3.*ME2*ME2*m_s*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) + 
  14.*ME2*MP2*m_s*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) + 
  3.*MP2*MP2*m_s*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) - 
  3.*ME2*m_s*m_s*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) + 
  16.*MP2*MP2*m_t*Conjugate(C0i(cc1, MP2, m_s, MP2, 0, MP2, MP2)) - 
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
  6.*ME2*MP2*m_s*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) + 
  MP2*MP2*m_s*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) - 
  ME2*m_s*m_s*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) + 
  8.*MP2*MP2*m_t*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) - 
  6.*MP2*m_s*m_t*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) + 
  m_s*m_s*m_t*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) - 
  4.*MP2*m_t*m_t*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) + 
  m_s*m_t*m_t*Conjugate(C0i(cc11, MP2, m_s, MP2, 0, MP2, MP2)) + 
  8.*ME2*MP2*MP2*Conjugate(C0i(cc12, ME2, m_s, ME2, 0, ME2, ME2)) - 
  8.*ME2*MP2*m_s*Conjugate(C0i(cc12, ME2, m_s, ME2, 0, ME2, ME2)) + 
  2.*ME2*m_s*m_s*Conjugate(C0i(cc12, ME2, m_s, ME2, 0, ME2, ME2)) - 
  16.*ME2*MP2*m_t*Conjugate(C0i(cc12, ME2, m_s, ME2, 0, ME2, ME2)) + 
  8.*ME2*m_s*m_t*Conjugate(C0i(cc12, ME2, m_s, ME2, 0, ME2, ME2)) + 
  8.*ME2*m_t*m_t*Conjugate(C0i(cc12, ME2, m_s, ME2, 0, ME2, ME2)) + 
  8.*ME2*ME2*MP2*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) - 
  16.*ME2*MP2*MP2*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8.*MP2*MP2*MP2*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) - 
  2.*ME2*ME2*m_s*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) + 
  12.*ME2*MP2*m_s*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) + 
  2.*MP2*MP2*m_s*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) - 
  2.*ME2*m_s*m_s*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) + 
  16.*MP2*MP2*m_t*Conjugate(C0i(cc12, MP2, m_s, MP2, 0, MP2, MP2)) - 
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
  16.*ME2*MP2*MP2*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) - 
  8.*MP2*MP2*MP2*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) - 
  3.*ME2*ME2*m_s*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) + 
  14.*ME2*MP2*m_s*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) + 
  3.*MP2*MP2*m_s*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) - 
  3.*ME2*m_s*m_s*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) + 
  16.*MP2*MP2*m_t*Conjugate(C0i(cc2, MP2, m_s, MP2, 0, MP2, MP2)) - 
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
  6.*ME2*MP2*m_s*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) + 
  MP2*MP2*m_s*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) - 
  ME2*m_s*m_s*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) + 
  8.*MP2*MP2*m_t*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) - 
  6.*MP2*m_s*m_t*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) + 
  m_s*m_s*m_t*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) - 
  4.*MP2*m_t*m_t*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) + 
  m_s*m_t*m_t*Conjugate(C0i(cc22, MP2, m_s, MP2, 0, MP2, MP2)) - 
  2.*ME2*ME2*Re(B0i(bb0, ME2, 0, ME2)) + 4.*ME2*MP2*Re(B0i(bb0, ME2, 0, ME2)) + 
  2.*MP2*MP2*Re(B0i(bb0, ME2, 0, ME2)) - 2.*ME2*m_s*Re(B0i(bb0, ME2, 0, ME2)) - 
  4.*MP2*m_t*Re(B0i(bb0, ME2, 0, ME2)) + 2.*m_s*m_t*Re(B0i(bb0, ME2, 0, ME2)) + 
  2.*m_t*m_t*Re(B0i(bb0, ME2, 0, ME2)) + 3.*ME2*ME2*Re(B0i(bb0, MP2, 0, MP2)) - 
  6.*ME2*MP2*Re(B0i(bb0, MP2, 0, MP2)) - 3.*MP2*MP2*Re(B0i(bb0, MP2, 0, MP2)) + 
  3.*ME2*m_s*Re(B0i(bb0, MP2, 0, MP2)) + 6.*MP2*m_t*Re(B0i(bb0, MP2, 0, MP2)) - 
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
   PRINT_VAR(pi.m_loopgenerator);
  if (pi.m_loopgenerator.find("Internal")!=0) return NULL;
  // if (pi.m_fi.m_nlotype==nlo_type::loop) {
  Flavour_Vector fl(pi.ExtractFlavours());
   PRINT_INFO(fl);
  if (fl.size()!=4) return NULL;
  if ( (fl[0]==Flavour(kf_e) || fl[1]==Flavour(kf_e))  && fl[1]==fl[0].Bar() &&
      ((fl[2].Kfcode()==kf_pi_plus || fl[3].Kfcode()==kf_pi_plus)) && fl[3]==fl[2].Bar())
  {
    return new PionPionVirtual(pi,fl,0.0,0);
  }
  return NULL;
}