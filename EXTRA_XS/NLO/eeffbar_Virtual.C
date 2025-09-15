#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "METOOLS/Loops/Master_Integrals.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "clooptools.h"
#include <iostream>

using namespace PHASIC;
using namespace ATOOLS;
using namespace METOOLS;
using namespace MODEL;
using namespace std;

namespace EXTRAXS {
  class eeffbar_Virtual : public PHASIC::Virtual_ME2_Base {
    double m_eps2, m_eps, m_fin;
    double ME2, MP2, m_alpha, ML2, MM2, MM, ME, ML, MZ, MZ2;
    double m_s, m_t, m_u, CW, SW, SW2, CW2;
    double dZAA1, dZZA1;
  public:
    eeffbar_Virtual(const Process_Info& pi, const Flavour_Vector& flavs,
                  const double& ep2, const double& ep);

    ~eeffbar_Virtual() {
    }

    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    Complex BornTriangle();
    Complex BornBox();
    Complex BornSelf();
    Complex Full();
    inline double Den(const double &a, const double &b) {return (a==b?0.:1./(a-b));}
    inline double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom) {return 4*M_PI;}
    // Complex MBub(double a, double b, double c);
    // Complex Mm_tri(double a, double b, double c,
    //             double aa, double bb, double cc);
    // Complex MBox(double a, double b, double c);
    // inline Complex Conjugate(Complex a) {return conj(a);}
    // inline double Re(Complex a) {return a.real();}

  };
}

using namespace EXTRAXS;

eeffbar_Virtual::eeffbar_Virtual(const Process_Info& pi, const Flavour_Vector& flavs,
                  const double& ep2, const double& ep) :
      Virtual_ME2_Base(pi, flavs),
      m_eps2(ep2), m_eps(ep)
   {
      // if (!s_loader->LoadLibrary("ooptools")) m_tHROW(fatal_error, "Failed to load libooptools.");
      ltini();
      Setlambda(0.);
      ME = Flavour(kf_e).Mass();
      ME2 = ME*ME;
      MM = Flavour(kf_mu).Mass();
      MM2 = MM*MM;
      ML = Flavour(kf_tau).Mass();
      ML2 = ML*ML;
      MZ = Flavour(kf_Z).Mass();
      MZ2 = MZ*MZ;
      double  MW= Flavour(kf_Wplus).Mass();
      double MW2 = MW*MW;

      MP2 = (Flavour(kf_pi).Mass()*Flavour(kf_pi).Mass());
      SW2 = std::abs(MODEL::s_model->ComplexConstant(string("csin2_thetaW")));
      CW2 = 1.-SW*SW;
      CW = sqrt(CW2);
      SW = sqrt(SW2);
      // double MU = 0.1;
      double MU =Flavour(kf_u).Mass(); double MU2 = MU*MU;
      double MD =Flavour(kf_d).Mass(); double MD2 = MD*MD;
      double MS =Flavour(kf_s).Mass(); double MS2 = MS*MS;
      double MC =Flavour(kf_c).Mass(); double MC2 = MC*MC;
      double MB =Flavour(kf_b).Mass(); double MB2 = MB*MB;
      double MT =Flavour(kf_t).Mass(); double MT2 = MT*MT;
      m_mode=2;
      m_IRscale=1;
      m_UVscale=1.;
      // PRINT_VAR(getdelta());
      // setdelta(-2);
      // setmudim(1);
      // setuvdiv(0);
     //  dZAA1=(m_alpha*
     //         (1./6. + (5*Re(B0i(bb0,0,MW2,MW2)))/4. + 
     //           Re(B0i(bb1,0,ME2,ME2)) + Re(B0i(bb1,0,ML2,ML2)) + 
     //           Re(B0i(bb1,0,MM2,MM2)) + 
     //           (Re(B0i(bb1,0,MB2,MB2)) + Re(B0i(bb1,0,MD2,MD2)) + 
     //              Re(B0i(bb1,0,MS2,MS2)))/3. + 
     //           (4*(Re(B0i(bb1,0,MC2,MC2)) + Re(B0i(bb1,0,MT2,MT2)) + 
     //                Re(B0i(bb1,0,MU2,MU2))))/3. + Re(B0i(bb1,0,MW2,MW2))/2. - 
     //           2*(Re(B0i(dbb00,0,ME2,ME2)) + Re(B0i(dbb00,0,ML2,ML2)) + 
     //              Re(B0i(dbb00,0,MM2,MM2))) - 
     //           (2*(Re(B0i(dbb00,0,MB2,MB2)) + Re(B0i(dbb00,0,MD2,MD2)) + 
     //                Re(B0i(dbb00,0,MS2,MS2))))/3. - 
     //           (8*(Re(B0i(dbb00,0,MC2,MC2)) + Re(B0i(dbb00,0,MT2,MT2)) + 
     //                Re(B0i(dbb00,0,MU2,MU2))))/3. + 3*Re(B0i(dbb00,0,MW2,MW2))))/M_PI;

     //  dZZA1=(2*m_alpha*(-1/8.*(((-4*CW*MW2)/SW - (4*MW2*SW)/CW)*
     //        Re(B0i(bb0,0,MW2,MW2))) - 
     //     ((2*SW - (1 - 2*SW2)/SW)*
     //         ((Re(A0(ME2)) + Re(A0(ML2)) + Re(A0(MM2)))/4. + 
     //           (-Re(B0i(bb00,0,ME2,ME2)) - Re(B0i(bb00,0,ML2,ML2)) - 
     //              Re(B0i(bb00,0,MM2,MM2)))/2.) + 
     //        (2*SW - (3 - 2*SW2)/SW)*
     //         ((Re(A0(MB2)) + Re(A0(MD2)) + Re(A0(MS2)))/12. + 
     //           (-Re(B0i(bb00,0,MB2,MB2)) - Re(B0i(bb00,0,MD2,MD2)) - 
     //              Re(B0i(bb00,0,MS2,MS2)))/6.) + 
     //        (4*SW - (3 - 4*SW2)/SW)*
     //         ((Re(A0(MC2)) + Re(A0(MT2)) + Re(A0(MU2)))/6. + 
     //           (-Re(B0i(bb00,0,MC2,MC2)) - Re(B0i(bb00,0,MT2,MT2)) - 
     //              Re(B0i(bb00,0,MU2,MU2)))/3.))/CW - 
     //     (-1/4.*((1/CW - 6*CW)*Re(A0(MW2))) - 
     //        ((8*CW - (3*CW2 + SW2)/CW)*Re(B0i(bb00,0,MW2,MW2)))/2.)/SW))/
     // (MZ2*M_PI);
      // setmudim(m_IRscale);
      // setdelta(4.*M_PI);
      // m_setminmass(0.);
   }


void eeffbar_Virtual::Calc(const Vec4D_Vector& momenta) {
  double factor(1.);
  if      (m_stype&sbt::qcd) factor=2.*M_PI/AlphaQCD();
  else if (m_stype&sbt::qed) factor=2.*M_PI/AlphaQED();
  else THROW(fatal_error,"m_unknown coupling.");
  Setlambda(0.);
  m_s = (momenta[0]+momenta[1]).Abs2();
  m_t = (momenta[0]-momenta[2]).Abs2();
  m_u = (momenta[0]-momenta[3]).Abs2();
  m_alpha = AlphaQED();
  m_res.Finite()=  2*(Full()).real()/m_alpha;
  Setlambda(-1.);
  m_res.IR() = 4*M_PI*(Full()).real()/m_alpha;
  clearcache();// m_since in general s will be different we gain nothing from cache
}
  
Complex eeffbar_Virtual::Full(){
   // return BornTriangle();
   return BornBox()+BornTriangle()+BornSelf();
}

Complex eeffbar_Virtual::BornBox(){
return 
4.*(2.*m_alpha*m_alpha*(-1/2.*(m_alpha*ME2*MM2*M_PI*(ME2 + MM2 - m_t))/m_s + 
       (m_alpha*ME2*MM2*M_PI*(-ME2 - MM2 + m_s + m_t))/m_s - 
       (m_alpha*M_PI*(pow(ME2,3) + pow(MM2,3) - 
            (ME2*ME2 - ME2*MM2 + MM2*MM2 + 
               (ME2 + MM2)*(2.*ME2 + 2.*MM2 - m_s - m_t))*m_t + 
            (2.*ME2 + 2.*MM2 - m_s - m_t)*(m_t*m_t)))/(2.*m_s))*
     (8.*(D0i(dd1,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd1,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd11,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd11,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd12,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd12,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd13,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd13,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
       4.*D0i(dd23,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
    2.*m_alpha*m_alpha*ME*((m_alpha*ME*MM2*M_PI*(ME2 + MM2 - m_t))/m_s - 
       (m_alpha*ME*MM2*M_PI*(-ME2 - MM2 + m_s + m_t))/m_s)*
     (2.*(D0i(dd2,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd2,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
       4.*(D0i(dd12,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd12,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd22,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd22,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd23,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2))) + 
    2.*m_alpha*m_alpha*ME*((3.*m_alpha*ME*MM2*M_PI*(ME2 + MM2 - m_t))/m_s + 
       (3.*m_alpha*ME*MM2*M_PI*(-ME2 - MM2 + m_s + m_t))/m_s)*
     (2.*(D0i(dd2,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd2,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
       4.*(D0i(dd12,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd12,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd22,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd22,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd23,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2))) + 
    8.*m_alpha*m_alpha*ME*((m_alpha*ME*ME2*MM2*M_PI)/m_s + 
       (m_alpha*ME*MM2*M_PI*(2.*ME2 - m_s))/(2.*m_s) + 
       (m_alpha*ME*M_PI*(ME2*ME2 + MM2*MM2 - m_t*m_u))/(2.*m_s) + 
       (m_alpha*ME*M_PI*(ME2*ME2 - MM2*MM2 - m_t*m_u + MM2*(-2.*ME2 + m_t + m_u)))/(2.*m_s))*
     (D0i(dd12,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
       D0i(dd12,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
       D0i(dd2,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
       D0i(dd22,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
       D0i(dd23,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
    2.*m_alpha*m_alpha*ME*MM*((12.*m_alpha*ME*MM*M_PI*(ME2 + MM2 - m_t))/m_s + 
       (12.*m_alpha*ME*MM*M_PI*(-ME2 - MM2 + m_s + m_t))/m_s)*
     (-D0i(dd23,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
       D0i(dd23,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2)) - 
    2.*m_alpha*m_alpha*((3.*m_alpha*ME2*MM2*M_PI*(ME2 + MM2 - m_t))/m_s + 
       (m_alpha*MM2*M_PI*(-3.*(ME2*ME2 + ME2*MM2) + (2.*ME2 + MM2)*m_s + 
            (3.*ME2 - m_s)*m_t))/m_s + 
       (m_alpha*ME2*M_PI*(-3.*(ME2*MM2 + MM2*MM2) + (ME2 + 2.*MM2)*m_s + 
            (3.*MM2 - m_s)*m_t))/m_s - 
       (m_alpha*M_PI*(-pow(ME2,3) - pow(MM2,3) - 
            2.*(ME2*ME2*MM2 + ME2*(MM2*MM2)) - 
            (ME2 + MM2 - m_t)*pow(2.*ME2 + 2.*MM2 - m_s - m_t,2) + 
            (ME2*ME2 + ME2*MM2 + MM2*MM2)*m_t + 
            (2.*ME2 + 2.*MM2 - m_s - m_t)*
             (2.*(ME2*ME2) + 4.*ME2*MM2 + 2.*(MM2*MM2) - 2.*(ME2 + MM2)*m_t)))/m_s)*
     (-D0i(dd23,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
       D0i(dd23,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
    2.*m_alpha*m_alpha*(-((m_alpha*ME2*MM2*M_PI*(ME2 + MM2 - m_t))/m_s) + 
       (m_alpha*ME2*MM2*M_PI*(-ME2 - MM2 + m_s + m_t))/m_s)*
     (-2.*D0i(dd23,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
       2.*D0i(dd23,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2)) - 
    2.*m_alpha*m_alpha*(-1/2.*(m_alpha*ME2*MM2*M_PI*(ME2 + MM2 - m_t))/m_s + 
       (m_alpha*ME2*MM2*M_PI*(-ME2 - MM2 + m_s + m_t))/m_s - 
       (m_alpha*M_PI*(-pow(ME2,3) - pow(MM2,3) + 
            2.*(ME2*ME2*MM2 + ME2*(MM2*MM2)) - 
            (ME2*MM2 + pow(2.*ME2 + 2.*MM2 - m_s - m_t,2))*m_t + 
            (2.*ME2 + 2.*MM2 - m_s - m_t)*
             (ME2*ME2 - 2.*ME2*MM2 + MM2*MM2 + (ME2 + MM2)*m_t)))/(2.*m_s))*
     (-8.*(D0i(dd1,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd1,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd11,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd11,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd12,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd12,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd13,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd13,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd23,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
       4.*D0i(dd23,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
    2.*m_alpha*m_alpha*ME*((m_alpha*ME*ME2*MM2*M_PI)/m_s + 
       (m_alpha*ME*MM2*M_PI*(2.*ME2 - m_s))/(2.*m_s) + 
       (m_alpha*ME*M_PI*(ME2*ME2 + MM2*MM2 - m_t*m_u))/(2.*m_s) + 
       (m_alpha*ME*M_PI*(ME2*ME2 - MM2*MM2 - m_t*m_u + MM2*(-2.*ME2 + m_t + m_u)))/(2.*m_s))*
     (8.*D0i(dd22,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
       4.*(D0i(dd12,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd12,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd2,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd22,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd23,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2))) + 
    2.*m_alpha*m_alpha*ME*((m_alpha*ME*MM2*M_PI*(ME2 + MM2 - m_t))/m_s - 
       (m_alpha*ME*MM2*M_PI*(-ME2 - MM2 + m_s + m_t))/m_s)*
     (2.*(D0i(dd2,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd2,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
       4.*(D0i(dd12,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd12,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd22,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd22,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd23,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2))) + 
    2.*m_alpha*m_alpha*ME*((3.*m_alpha*ME*MM2*M_PI*(ME2 + MM2 - m_t))/(2.*m_s) + 
       (3.*m_alpha*ME*MM2*M_PI*(-ME2 - MM2 + m_s + m_t))/(2.*m_s) + 
       (m_alpha*ME*M_PI*(-3.*(ME2*MM2 + MM2*MM2) + (2.*ME2 + MM2)*m_s + 
            (3.*MM2 - 2.*m_s)*m_t))/(2.*m_s) + 
       (m_alpha*ME*M_PI*(-3.*(ME2*MM2 + MM2*MM2) + (2.*ME2 + MM2)*m_s + 
            (3.*MM2 - 2.*m_s)*m_u))/(2.*m_s))*
     (2.*(D0i(dd2,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd2,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
       4.*(D0i(dd12,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd12,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd22,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd22,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd23,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2))) - 
    2.*m_alpha*m_alpha*((-4.*m_alpha*ME2*MM2*M_PI)/m_s + (m_alpha*MM2*M_PI*(2.*ME2 - m_s))/m_s + 
       (m_alpha*ME2*M_PI*(2.*MM2 - m_s))/m_s - 
       (m_alpha*M_PI*(ME2*ME2 + 2.*ME2*MM2 + MM2*MM2 - 
            2.*(ME2 + MM2)*(2.*ME2 + 2.*MM2 - m_s - m_t) + 
            pow(2.*ME2 + 2.*MM2 - m_s - m_t,2)))/m_s)*
     (-16.*D0i(dd00,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) - 
       12.*D0i(dd00,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
       4.*(C0i(cc0,m_s,MM2,MM2,0,0,MM2) + 
          (ME2 + 2.*MM2 - m_s - m_t)*
           D0i(dd13,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2)) - 
       4.*((-ME2 - MM2 + m_s + m_t)*
           (D0i(dd11,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,
              MM2) + D0i(dd11,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
          (-MM2 + m_s + m_t)*(D0i(dd12,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,
              MM2,0,0,ME2,MM2) + D0i(dd12,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2)) \
+ (-ME2 + m_s + m_t)*D0i(dd13,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2)) - 
       8.*ME2*(D0i(dd22,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,
           MM2) + D0i(dd22,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2)) - 
       (ME2 - 3.*MM2 + 3.*m_t - 4.*m_u)*
        D0i(dd23,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) - 
       (7.*(ME2 + MM2) - 3.*m_t)*D0i(dd23,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) - 
       2.*((-ME2 - MM2 + m_s + m_t)*
           D0i(dd0,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          (-ME2 - MM2 + m_t)*D0i(dd0,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          (-2.*ME2 - 2.*MM2 + 2.*m_s + 2.*m_t)*
           (D0i(dd1,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
             D0i(dd1,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
          (-MM2 + m_s + m_t)*D0i(dd2,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,
            0,ME2,MM2) + (6.*ME2 + MM2 - m_s - m_u)*
           D0i(dd2,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          (-ME2 - 2.*MM2 + m_s + m_t)*
           D0i(dd3,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          (ME2 + 2.*MM2 - m_t)*D0i(dd3,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
       MM2*(-4.*D0i(dd33,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,
            MM2) - 8.*D0i(dd33,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2))) - 
    2.*m_alpha*m_alpha*MM*(-((m_alpha*ME2*MM*MM2*M_PI)/m_s) - 
       (m_alpha*ME2*MM*M_PI*(2.*MM2 - m_s))/(2.*m_s) - 
       (m_alpha*MM*M_PI*(ME2*ME2 + MM2*MM2 - m_t*m_u))/(2.*m_s) + 
       (m_alpha*MM*M_PI*(ME2*ME2 - MM2*MM2 + ME2*(2.*MM2 - m_t - m_u) + m_t*m_u))/(2.*m_s))*
     (-4.*(D0i(dd13,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd23,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd3,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
       4.*(D0i(dd13,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd23,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd3,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
       8.*D0i(dd33,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
    2.*m_alpha*m_alpha*MM*((-3.*m_alpha*ME2*MM*M_PI*(ME2 + MM2 - m_t))/m_s - 
       (3.*m_alpha*ME2*MM*M_PI*(-ME2 - MM2 + m_s + m_t))/m_s)*
     (4.*(D0i(dd13,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd23,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
       2.*D0i(dd3,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) - 
       2.*D0i(dd3,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) - 
       4.*(D0i(dd13,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd23,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd33,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2))) + 
    4.*m_alpha*m_alpha*MM*(-((m_alpha*ME2*MM*M_PI*(ME2 + MM2 - m_t))/m_s) + 
       (m_alpha*ME2*MM*M_PI*(-ME2 - MM2 + m_s + m_t))/m_s)*
     (4.*(D0i(dd13,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd23,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
       2.*D0i(dd3,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) - 
       2.*D0i(dd3,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) - 
       4.*(D0i(dd13,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd23,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd33,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2))) + 
    2.*m_alpha*m_alpha*MM*((-3.*m_alpha*ME2*MM*M_PI*(ME2 + MM2 - m_t))/(2.*m_s) - 
       (3.*m_alpha*ME2*MM*M_PI*(-ME2 - MM2 + m_s + m_t))/(2.*m_s) - 
       (m_alpha*MM*M_PI*(-3.*(ME2*ME2 + ME2*MM2) + (ME2 + 2.*MM2)*m_s + 
            (3.*ME2 - 2.*m_s)*m_t))/(2.*m_s) - 
       (m_alpha*MM*M_PI*(-3.*(ME2*ME2 + ME2*MM2) + (ME2 + 2.*MM2)*m_s + 
            (3.*ME2 - 2.*m_s)*m_u))/(2.*m_s))*
     (4.*(D0i(dd13,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd23,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
       2.*D0i(dd3,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) - 
       2.*D0i(dd3,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) - 
       4.*(D0i(dd13,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd23,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd33,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2))) - 
    2.*m_alpha*m_alpha*MM*(-((m_alpha*ME2*MM*MM2*M_PI)/m_s) - 
       (m_alpha*ME2*MM*M_PI*(2.*MM2 - m_s))/(2.*m_s) - 
       (m_alpha*MM*M_PI*(ME2*ME2 + MM2*MM2 - m_t*m_u))/(2.*m_s) + 
       (m_alpha*MM*M_PI*(ME2*ME2 - MM2*MM2 + ME2*(2.*MM2 - m_t - m_u) + m_t*m_u))/(2.*m_s))*
     (-4.*(D0i(dd13,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd23,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd3,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd33,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
       4.*(D0i(dd13,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd23,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd3,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd33,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2))) + 
    4.*m_alpha*m_alpha*ME*MM*((m_alpha*ME*MM*M_PI*(ME2 + MM2 - m_t))/m_s - 
       (m_alpha*ME*MM*M_PI*(-ME2 - MM2 + m_s + m_t))/m_s)*
     (4.*(D0i(dd13,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd23,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2)) - 
       4.*(D0i(dd12,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd12,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd13,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd2,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd2,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd22,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd22,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd23,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd3,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          D0i(dd33,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2))) - 
    2.*m_alpha*m_alpha*((-4.*m_alpha*ME2*MM2*M_PI)/m_s + (m_alpha*MM2*M_PI*(2.*ME2 - m_s))/m_s + 
       (m_alpha*ME2*M_PI*(2.*MM2 - m_s))/m_s - 
       (m_alpha*M_PI*(ME2*ME2 + 2.*ME2*MM2 + MM2*MM2 - 2.*(ME2 + MM2)*m_t + m_t*m_t))/m_s)*
     (4.*(C0i(cc0,m_s,MM2,MM2,0,0,MM2) + 
          (ME2 + 2.*MM2 - m_s - m_t)*
           D0i(dd13,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2)) - 
       2.*((-ME2 - MM2 + m_s + m_t)*
           D0i(dd0,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          (-ME2 - MM2 + m_t)*D0i(dd0,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          (-2.*ME2 - 2.*MM2 + 2.*m_s + 2.*m_t)*
           (D0i(dd1,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
             D0i(dd1,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
          (-MM2 + m_s + m_t)*D0i(dd2,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,
            0,ME2,MM2) + (6.*ME2 + MM2 - m_s - m_u)*
           D0i(dd2,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          (-ME2 - 2.*MM2 + m_s + m_t)*
           D0i(dd3,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          (ME2 + 2.*MM2 - m_t)*D0i(dd3,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2)) - 
       4.*(D0i(dd00,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          (-ME2 - MM2 + m_s + m_t)*
           (D0i(dd11,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
             D0i(dd11,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
          (-MM2 + m_s + m_t)*(D0i(dd12,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,
              0,0,ME2,MM2) + D0i(dd12,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
          (-ME2 + m_s + m_t)*D0i(dd13,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          ME2*(D0i(dd22,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,
              MM2) + D0i(dd22,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2)) + 
          (-ME2 - 2.*MM2 + m_s + m_t)*
           D0i(dd23,m_s,ME2,2.*ME2 + 2.*MM2 - m_s - m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          (ME2 + MM2)*D0i(dd23,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2) + 
          MM2*D0i(dd33,m_s,ME2,m_t,MM2,ME2,MM2,0,0,ME2,MM2))));
}

Complex eeffbar_Virtual::BornSelf(){
   return 4.*(2.*((-4.*m_alpha*ME2*MM2*M_PI)/m_s + (m_alpha*MM2*M_PI*(2.*ME2 - m_s))/m_s + 
       (m_alpha*ME2*M_PI*(2.*MM2 - m_s))/m_s - 
       (m_alpha*M_PI*(ME2*ME2 + 2.*ME2*MM2 + MM2*MM2 - 
            2.*(ME2 + MM2)*(2.*ME2 + 2.*MM2 - m_s - m_t) + 
            pow(2.*ME2 + 2.*MM2 - m_s - m_t,2)))/m_s)*
     (-(m_alpha*m_alpha*((4.*A0(ME2) - 8.*B0i(bb00,m_s,ME2,ME2))/(m_s*m_s) + 
            (4.*B0i(bb1,m_s,ME2,ME2))/m_s)) - 
       m_alpha*m_alpha*((4.*A0(ML2) - 8.*B0i(bb00,m_s,ML2,ML2))/(m_s*m_s) + 
          (4.*B0i(bb1,m_s,ML2,ML2))/m_s) - 
       m_alpha*m_alpha*((4.*A0(MM2) - 8.*B0i(bb00,m_s,MM2,MM2))/(m_s*m_s) + 
          (4.*B0i(bb1,m_s,MM2,MM2))/m_s) - 
       (4.*m_alpha*m_alpha*(-Re(B0i(bb1,0,ME2,ME2)) + 2.*Re(B0i(dbb00,0,ME2,ME2))))/m_s - 
       (4.*m_alpha*m_alpha*(-Re(B0i(bb1,0,ML2,ML2)) + 2.*Re(B0i(dbb00,0,ML2,ML2))))/m_s - 
       (4.*m_alpha*m_alpha*(-Re(B0i(bb1,0,MM2,MM2)) + 2.*Re(B0i(dbb00,0,MM2,MM2))))/m_s) + 
    2.*((-4.*m_alpha*ME2*MM2*M_PI)/m_s + (m_alpha*MM2*M_PI*(2.*ME2 - m_s))/m_s + 
       (m_alpha*ME2*M_PI*(2.*MM2 - m_s))/m_s - 
       (m_alpha*M_PI*(ME2*ME2 + 2.*ME2*MM2 + MM2*MM2 - 2.*(ME2 + MM2)*m_t + m_t*m_t))/m_s)*
     (-(m_alpha*m_alpha*((4.*A0(ME2) - 8.*B0i(bb00,m_s,ME2,ME2))/(m_s*m_s) + 
            (4.*B0i(bb1,m_s,ME2,ME2))/m_s)) - 
       m_alpha*m_alpha*((4.*A0(ML2) - 8.*B0i(bb00,m_s,ML2,ML2))/(m_s*m_s) + 
          (4.*B0i(bb1,m_s,ML2,ML2))/m_s) - 
       m_alpha*m_alpha*((4.*A0(MM2) - 8.*B0i(bb00,m_s,MM2,MM2))/(m_s*m_s) + 
          (4.*B0i(bb1,m_s,MM2,MM2))/m_s) - 
       (4.*m_alpha*m_alpha*(-Re(B0i(bb1,0,ME2,ME2)) + 2.*Re(B0i(dbb00,0,ME2,ME2))))/m_s - 
       (4.*m_alpha*m_alpha*(-Re(B0i(bb1,0,ML2,ML2)) + 2.*Re(B0i(dbb00,0,ML2,ML2))))/m_s - 
       (4.*m_alpha*m_alpha*(-Re(B0i(bb1,0,MM2,MM2)) + 2.*Re(B0i(dbb00,0,MM2,MM2))))/m_s));
}


Complex eeffbar_Virtual::BornTriangle(){
return 4.*((4.*((m_alpha*ME*ME2*MM2*M_PI)/m_s + (m_alpha*ME*MM2*M_PI*(2*ME2 - m_s))/(2.*m_s) + 
 (m_alpha*ME*M_PI*(ME2*ME2 + MM2*MM2 - m_t*m_u))/(2.*m_s) + 
 (m_alpha*ME*M_PI*(ME2*ME2 - MM2*MM2 - m_t*m_u + MM2*(-2*ME2 + m_t + m_u)))/
  (2.*m_s))*(-8*m_alpha*m_alpha*ME*C0i(cc12,ME2,m_s,ME2,0,ME2,ME2) - 
 4*m_alpha*m_alpha*ME*(C0i(cc1,ME2,m_s,ME2,0,ME2,ME2) + 
    C0i(cc11,ME2,m_s,ME2,0,ME2,ME2) + C0i(cc2,ME2,m_s,ME2,0,ME2,ME2) + 
    C0i(cc22,ME2,m_s,ME2,0,ME2,ME2))))/m_s + 
    (4*(-((m_alpha*ME2*MM*MM2*M_PI)/m_s) - (m_alpha*ME2*MM*M_PI*(2*MM2 - m_s))/(2.*m_s) - 
 (m_alpha*MM*M_PI*(ME2*ME2 + MM2*MM2 - m_t*m_u))/(2.*m_s) + 
 (m_alpha*MM*M_PI*(ME2*ME2 - MM2*MM2 + ME2*(2*MM2 - m_t - m_u) + m_t*m_u))/
  (2.*m_s))*(8*m_alpha*m_alpha*MM*C0i(cc12,MM2,m_s,MM2,0,MM2,MM2) + 
 4*m_alpha*m_alpha*MM*(C0i(cc1,MM2,m_s,MM2,0,MM2,MM2) + 
    C0i(cc11,MM2,m_s,MM2,0,MM2,MM2) + C0i(cc2,MM2,m_s,MM2,0,MM2,MM2) + 
    C0i(cc22,MM2,m_s,MM2,0,MM2,MM2))))/m_s + 
    (2*((-4*m_alpha*ME2*MM2*M_PI)/m_s + (m_alpha*MM2*M_PI*(2*ME2 - m_s))/m_s + 
 (m_alpha*ME2*M_PI*(2*MM2 - m_s))/m_s - 
 (m_alpha*M_PI*(ME2*ME2 + 2*ME2*MM2 + MM2*MM2 - 2*(ME2 + MM2)*m_t + m_t*m_t))/
  m_s)*(2*m_alpha*m_alpha - 2*m_alpha*m_alpha*(B0i(bb0,m_s,ME2,ME2) + B0i(bb0,m_s,MM2,MM2)) + 
 2*m_alpha*m_alpha*(-2*ME2 + m_s)*C0i(cc0,ME2,m_s,ME2,0,ME2,ME2) + 
 2*m_alpha*m_alpha*(-2*MM2 + m_s)*C0i(cc0,MM2,m_s,MM2,0,MM2,MM2) + 
 4*m_alpha*m_alpha*(C0i(cc00,ME2,m_s,ME2,0,ME2,ME2) + 
    C0i(cc00,MM2,m_s,MM2,0,MM2,MM2)) + 
 2*m_alpha*m_alpha*(-4*ME2 + m_s)*(C0i(cc1,ME2,m_s,ME2,0,ME2,ME2) + 
    C0i(cc2,ME2,m_s,ME2,0,ME2,ME2)) + 
 2*m_alpha*m_alpha*(-4*MM2 + m_s)*(C0i(cc1,MM2,m_s,MM2,0,MM2,MM2) + 
    C0i(cc2,MM2,m_s,MM2,0,MM2,MM2)) - 
 m_alpha*M_PI*(2*((m_alpha*(1/4. + 
    (-Re(B0i(bb0,ME2,0,ME2)) - Re(B0i(bb1,ME2,0,ME2)))/2.\
))/M_PI + (m_alpha*(1/4. + (-Re(B0i(bb0,MM2,0,MM2)) - Re(B0i(bb1,MM2,0,MM2)))/
     2.))/M_PI + (m_alpha*ME2*Re(B0i(dbb0,ME2,0,ME2)))/M_PI + 
       (m_alpha*MM2*Re(B0i(dbb0,MM2,0,MM2)))/M_PI - 
       (m_alpha*ME2*Re(B0i(dbb1,ME2,0,ME2)))/M_PI - 
       (m_alpha*MM2*Re(B0i(dbb1,MM2,0,MM2)))/M_PI) - 
    ((2*m_alpha*(2*SW - (1 - 2*SW2)/SW)*
  (-1/4.*Re(A0(ME2)) + Re(B0i(bb00,0,ME2,ME2))/2.))/
        (CW*MZ2*M_PI) - (8*m_alpha*(2*SW - (1 - 2*SW2)/SW)*SW2*
  (-1/4.*Re(A0(ME2)) + Re(B0i(bb00,0,ME2,ME2))/2.))/
        (CW*MZ2*M_PI) + (2*m_alpha*(2*SW - (1 - 2*SW2)/SW)*
  (-1/4.*Re(A0(ML2)) + Re(B0i(bb00,0,ML2,ML2))/2.))/
        (CW*MZ2*M_PI) - (8*m_alpha*(2*SW - (1 - 2*SW2)/SW)*SW2*
  (-1/4.*Re(A0(ML2)) + Re(B0i(bb00,0,ML2,ML2))/2.))/
        (CW*MZ2*M_PI) + (2*m_alpha*(2*SW - (1 - 2*SW2)/SW)*
  (-1/4.*Re(A0(MM2)) + Re(B0i(bb00,0,MM2,MM2))/2.))/
        (CW*MZ2*M_PI) - (8*m_alpha*(2*SW - (1 - 2*SW2)/SW)*SW2*
  (-1/4.*Re(A0(MM2)) + Re(B0i(bb00,0,MM2,MM2))/2.))/
        (CW*MZ2*M_PI) - 2*CW*SW*
        ((m_alpha*(1/4. + 
       (-Re(B0i(bb0,ME2,0,ME2)) - 
        Re(B0i(bb1,ME2,0,ME2)))/2.))/M_PI + 
  (m_alpha*(1/4. + 
       (-Re(B0i(bb0,MM2,0,MM2)) - 
        Re(B0i(bb1,MM2,0,MM2)))/2.))/M_PI + 
  (m_alpha*ME2*Re(B0i(dbb0,ME2,0,ME2)))/M_PI + 
  (m_alpha*MM2*Re(B0i(dbb0,MM2,0,MM2)))/M_PI - 
  (2*m_alpha*(-Re(B0i(bb1,0,ME2,ME2)) + 
       2*Re(B0i(dbb00,0,ME2,ME2))))/M_PI - 
  (2*m_alpha*(-Re(B0i(bb1,0,ML2,ML2)) + 
       2*Re(B0i(dbb00,0,ML2,ML2))))/M_PI - 
  (2*m_alpha*(-Re(B0i(bb1,0,MM2,MM2)) + 
       2*Re(B0i(dbb00,0,MM2,MM2))))/M_PI + 
  2*((-2*m_alpha*SW*(2*SW - (1 - 2*SW2)/SW)*
        (-1/4.*Re(A0(ME2)) + Re(B0i(bb00,0,ME2,ME2))/2.)\
)/(CW2*MZ2*M_PI) - (2*m_alpha*SW*(2*SW - (1 - 2*SW2)/SW)*
        (-1/4.*Re(A0(ML2)) + Re(B0i(bb00,0,ML2,ML2))/2.)\
)/(CW2*MZ2*M_PI) - (2*m_alpha*SW*(2*SW - (1 - 2*SW2)/SW)*
        (-1/4.*Re(A0(MM2)) + Re(B0i(bb00,0,MM2,MM2))/2.)\
)/(CW2*MZ2*M_PI) + (m_alpha*(-Re(B0i(bb1,0,ME2,ME2)) + 
         2*Re(B0i(dbb00,0,ME2,ME2))))/M_PI + 
     (m_alpha*(-Re(B0i(bb1,0,ML2,ML2)) + 
         2*Re(B0i(dbb00,0,ML2,ML2))))/M_PI + 
     (m_alpha*(-Re(B0i(bb1,0,MM2,MM2)) + 
         2*Re(B0i(dbb00,0,MM2,MM2))))/M_PI) - 
  (m_alpha*ME2*Re(B0i(dbb1,ME2,0,ME2)))/M_PI - 
  (m_alpha*MM2*Re(B0i(dbb1,MM2,0,MM2)))/M_PI))/(CW*SW))))/m_s + 
    (((-4*m_alpha*ME2*MM2*M_PI)/m_s + (m_alpha*MM2*M_PI*(2*ME2 - m_s))/m_s + 
 (m_alpha*ME2*M_PI*(2*MM2 - m_s))/m_s - 
 (m_alpha*M_PI*(ME2*ME2 + 2*ME2*MM2 + MM2*MM2 - 
      2*(ME2 + MM2)*(2*ME2 + 2*MM2 - m_s - m_t) + 
      pow(2*ME2 + 2*MM2 - m_s - m_t,2)))/m_s)*
       (-2*m_alpha*m_alpha*(B0i(bb0,m_s,ME2,ME2) + B0i(bb0,m_s,MM2,MM2)) + 
 2*m_alpha*m_alpha*(-2*ME2 + m_s)*C0i(cc0,ME2,m_s,ME2,0,ME2,ME2) + 
 2*m_alpha*m_alpha*(-2*MM2 + m_s)*C0i(cc0,MM2,m_s,MM2,0,MM2,MM2) + 
 4*m_alpha*m_alpha*(C0i(cc00,ME2,m_s,ME2,0,ME2,ME2) + 
    C0i(cc00,MM2,m_s,MM2,0,MM2,MM2)) + 
 2*m_alpha*m_alpha*(-4*ME2 + m_s)*(C0i(cc1,ME2,m_s,ME2,0,ME2,ME2) + 
    C0i(cc2,ME2,m_s,ME2,0,ME2,ME2)) + 
 2*m_alpha*m_alpha*(-4*MM2 + m_s)*(C0i(cc1,MM2,m_s,MM2,0,MM2,MM2) + 
    C0i(cc2,MM2,m_s,MM2,0,MM2,MM2)) + 
 2*(m_alpha*m_alpha - m_alpha*M_PI*((m_alpha*
  (1/4. + (-Re(B0i(bb0,ME2,0,ME2)) - 
       Re(B0i(bb1,ME2,0,ME2)))/2.))/M_PI + 
       (m_alpha*(1/4. + (-Re(B0i(bb0,MM2,0,MM2)) - 
       Re(B0i(bb1,MM2,0,MM2)))/2.))/M_PI + 
       (m_alpha*ME2*Re(B0i(dbb0,ME2,0,ME2)))/M_PI + 
       (m_alpha*MM2*Re(B0i(dbb0,MM2,0,MM2)))/M_PI - 
       (m_alpha*ME2*Re(B0i(dbb1,ME2,0,ME2)))/M_PI - 
       (m_alpha*MM2*Re(B0i(dbb1,MM2,0,MM2)))/M_PI + 
       ((4*m_alpha*SW*(2*SW - (1 - 2*SW2)/SW)*
     (-1/4.*Re(A0(ME2)) + Re(B0i(bb00,0,ME2,ME2))/2.))/
   (CW*MZ2*M_PI) + 
  (4*m_alpha*SW*(2*SW - (1 - 2*SW2)/SW)*
     (-1/4.*Re(A0(ML2)) + Re(B0i(bb00,0,ML2,ML2))/2.))/
   (CW*MZ2*M_PI) + 
  (4*m_alpha*SW*(2*SW - (1 - 2*SW2)/SW)*
     (-1/4.*Re(A0(MM2)) + Re(B0i(bb00,0,MM2,MM2))/2.))/
   (CW*MZ2*M_PI) + 
  CW*((m_alpha*(1/4. + 
         (-Re(B0i(bb0,ME2,0,ME2)) - 
         Re(B0i(bb1,ME2,0,ME2)))/2.))/M_PI + 
     (m_alpha*(1/4. + 
         (-Re(B0i(bb0,MM2,0,MM2)) - 
         Re(B0i(bb1,MM2,0,MM2)))/2.))/M_PI + 
     (m_alpha*ME2*Re(B0i(dbb0,ME2,0,ME2)))/M_PI + 
     (m_alpha*MM2*Re(B0i(dbb0,MM2,0,MM2)))/M_PI - 
     (2*m_alpha*(-Re(B0i(bb1,0,ME2,ME2)) + 
         2*Re(B0i(dbb00,0,ME2,ME2))))/M_PI - 
     (2*m_alpha*(-Re(B0i(bb1,0,ML2,ML2)) + 
         2*Re(B0i(dbb00,0,ML2,ML2))))/M_PI - 
     (2*m_alpha*(-Re(B0i(bb1,0,MM2,MM2)) + 
         2*Re(B0i(dbb00,0,MM2,MM2))))/M_PI + 
     2*((-2*m_alpha*SW*(2*SW - (1 - 2*SW2)/SW)*
         (-1/4.*Re(A0(ME2)) + 
        Re(B0i(bb00,0,ME2,ME2))/2.))/(CW2*MZ2*M_PI) - 
        (2*m_alpha*SW*(2*SW - (1 - 2*SW2)/SW)*
         (-1/4.*Re(A0(ML2)) + 
        Re(B0i(bb00,0,ML2,ML2))/2.))/(CW2*MZ2*M_PI) - 
        (2*m_alpha*SW*(2*SW - (1 - 2*SW2)/SW)*
         (-1/4.*Re(A0(MM2)) + 
        Re(B0i(bb00,0,MM2,MM2))/2.))/(CW2*MZ2*M_PI) + 
        (m_alpha*
         (-Re(B0i(bb1,0,ME2,ME2)) + 
         2*Re(B0i(dbb00,0,ME2,ME2))))/M_PI + 
        (m_alpha*
         (-Re(B0i(bb1,0,ML2,ML2)) + 
         2*Re(B0i(dbb00,0,ML2,ML2))))/M_PI + 
        (m_alpha*
         (-Re(B0i(bb1,0,MM2,MM2)) + 
         2*Re(B0i(dbb00,0,MM2,MM2))))/M_PI) - 
     (m_alpha*ME2*Re(B0i(dbb1,ME2,0,ME2)))/M_PI - 
     (m_alpha*MM2*Re(B0i(dbb1,MM2,0,MM2)))/M_PI))/CW))))/m_s + 
    (((-4*m_alpha*ME2*MM2*M_PI)/m_s + (m_alpha*MM2*M_PI*(2*ME2 - m_s))/m_s + 
 (m_alpha*ME2*M_PI*(2*MM2 - m_s))/m_s - 
 (m_alpha*M_PI*(ME2*ME2 + 2*ME2*MM2 + MM2*MM2 - 
      2*(ME2 + MM2)*(2*ME2 + 2*MM2 - m_s - m_t) + 
      pow(2*ME2 + 2*MM2 - m_s - m_t,2)))/m_s)*
       (-2*m_alpha*m_alpha*(B0i(bb0,m_s,ME2,ME2) + B0i(bb0,m_s,MM2,MM2)) + 
 2*m_alpha*m_alpha*(-2*ME2 + m_s)*C0i(cc0,ME2,m_s,ME2,0,ME2,ME2) + 
 2*m_alpha*m_alpha*(-2*MM2 + m_s)*C0i(cc0,MM2,m_s,MM2,0,MM2,MM2) + 
 4*m_alpha*m_alpha*(C0i(cc00,ME2,m_s,ME2,0,ME2,ME2) + 
    C0i(cc00,MM2,m_s,MM2,0,MM2,MM2)) + 
 2*m_alpha*m_alpha*(-4*ME2 + m_s)*(C0i(cc1,ME2,m_s,ME2,0,ME2,ME2) + 
    C0i(cc2,ME2,m_s,ME2,0,ME2,ME2)) + 
 2*m_alpha*m_alpha*(-4*MM2 + m_s)*(C0i(cc1,MM2,m_s,MM2,0,MM2,MM2) + 
    C0i(cc2,MM2,m_s,MM2,0,MM2,MM2)) + 
 2*(m_alpha*m_alpha - m_alpha*M_PI*((m_alpha*
  (1/4. + (-Re(B0i(bb0,ME2,0,ME2)) - 
       Re(B0i(bb1,ME2,0,ME2)))/2.))/M_PI + 
       (m_alpha*(1/4. + (-Re(B0i(bb0,MM2,0,MM2)) - 
       Re(B0i(bb1,MM2,0,MM2)))/2.))/M_PI + 
       (m_alpha*ME2*Re(B0i(dbb0,ME2,0,ME2)))/M_PI + 
       (m_alpha*MM2*Re(B0i(dbb0,MM2,0,MM2)))/M_PI - 
       (m_alpha*ME2*Re(B0i(dbb1,ME2,0,ME2)))/M_PI - 
       (m_alpha*MM2*Re(B0i(dbb1,MM2,0,MM2)))/M_PI - 
       ((2*m_alpha*(2*SW - (1 - 2*SW2)/SW)*
     (-1/4.*Re(A0(ME2)) + Re(B0i(bb00,0,ME2,ME2))/2.))/
   (CW*MZ2*M_PI) - 
  (4*m_alpha*(2*SW - (1 - 2*SW2)/SW)*SW2*
     (-1/4.*Re(A0(ME2)) + Re(B0i(bb00,0,ME2,ME2))/2.))/
   (CW*MZ2*M_PI) + 
  (2*m_alpha*(2*SW - (1 - 2*SW2)/SW)*
     (-1/4.*Re(A0(ML2)) + Re(B0i(bb00,0,ML2,ML2))/2.))/
   (CW*MZ2*M_PI) - 
  (4*m_alpha*(2*SW - (1 - 2*SW2)/SW)*SW2*
     (-1/4.*Re(A0(ML2)) + Re(B0i(bb00,0,ML2,ML2))/2.))/
   (CW*MZ2*M_PI) + 
  (2*m_alpha*(2*SW - (1 - 2*SW2)/SW)*
     (-1/4.*Re(A0(MM2)) + Re(B0i(bb00,0,MM2,MM2))/2.))/
   (CW*MZ2*M_PI) - 
  (4*m_alpha*(2*SW - (1 - 2*SW2)/SW)*SW2*
     (-1/4.*Re(A0(MM2)) + Re(B0i(bb00,0,MM2,MM2))/2.))/
   (CW*MZ2*M_PI) - 
  CW*SW*((m_alpha*
        (1/4. + 
         (-Re(B0i(bb0,ME2,0,ME2)) - 
         Re(B0i(bb1,ME2,0,ME2)))/2.))/M_PI + 
     (m_alpha*(1/4. + 
         (-Re(B0i(bb0,MM2,0,MM2)) - 
         Re(B0i(bb1,MM2,0,MM2)))/2.))/M_PI + 
     (m_alpha*ME2*Re(B0i(dbb0,ME2,0,ME2)))/M_PI + 
     (m_alpha*MM2*Re(B0i(dbb0,MM2,0,MM2)))/M_PI - 
     (2*m_alpha*(-Re(B0i(bb1,0,ME2,ME2)) + 
         2*Re(B0i(dbb00,0,ME2,ME2))))/M_PI - 
     (2*m_alpha*(-Re(B0i(bb1,0,ML2,ML2)) + 
         2*Re(B0i(dbb00,0,ML2,ML2))))/M_PI - 
     (2*m_alpha*(-Re(B0i(bb1,0,MM2,MM2)) + 
         2*Re(B0i(dbb00,0,MM2,MM2))))/M_PI + 
     2*((-2*m_alpha*SW*(2*SW - (1 - 2*SW2)/SW)*
         (-1/4.*Re(A0(ME2)) + Re(B0i(bb00,0,ME2,ME2))/2.)\
)/(CW2*MZ2*M_PI) - (2*m_alpha*SW*(2*SW - (1 - 2*SW2)/SW)*
         (-1/4.*Re(A0(ML2)) + Re(B0i(bb00,0,ML2,ML2))/2.)\
)/(CW2*MZ2*M_PI) - (2*m_alpha*SW*(2*SW - (1 - 2*SW2)/SW)*
         (-1/4.*Re(A0(MM2)) + Re(B0i(bb00,0,MM2,MM2))/2.)\
)/(CW2*MZ2*M_PI) + (m_alpha*(-Re(B0i(bb1,0,ME2,ME2)) + 
         2*Re(B0i(dbb00,0,ME2,ME2))))/M_PI + 
        (m_alpha*
         (-Re(B0i(bb1,0,ML2,ML2)) + 
         2*Re(B0i(dbb00,0,ML2,ML2))))/M_PI + 
        (m_alpha*
         (-Re(B0i(bb1,0,MM2,MM2)) + 
         2*Re(B0i(dbb00,0,MM2,MM2))))/M_PI) - 
     (m_alpha*ME2*Re(B0i(dbb1,ME2,0,ME2)))/M_PI - 
     (m_alpha*MM2*Re(B0i(dbb1,MM2,0,MM2)))/M_PI))/(CW*SW)))))/m_s);
}


DECLARE_VIRTUALME2_GETTER(EXTRAXS::eeffbar_Virtual,"eeffbar_Virtual")
Virtual_ME2_Base *ATOOLS::Getter
<PHASIC::Virtual_ME2_Base,PHASIC::Process_Info,EXTRAXS::eeffbar_Virtual>::
operator()(const Process_Info &pi) const
{
 

  PRINT_VAR(pi);
  if (pi.m_loopgenerator.find("Internal")!=0) return NULL;
  // if (pi.m_fi.m_nlotype==nlo_type::loop) {
  Flavour_Vector fl(pi.ExtractFlavours());
   PRINT_INFO(fl);
  if (fl.size()!=4) return NULL;
  if ( (fl[0]==Flavour(kf_e) || fl[1]==Flavour(kf_e))  && fl[1]==fl[0].Bar() &&
      ((fl[2].Kfcode()==kf_mu || fl[3].Kfcode()==kf_mu)) && fl[3]==fl[2].Bar())
  {
    return new eeffbar_Virtual(pi,fl,0.0,0);
  }
  return NULL;
}