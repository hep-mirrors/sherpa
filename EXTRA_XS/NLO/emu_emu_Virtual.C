#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "METOOLS/Loops/Master_Integrals.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "clooptools.h"

using namespace PHASIC;
using namespace ATOOLS;
using namespace METOOLS;
using namespace MODEL;

namespace EXTRAXS {
  class emu_emu_Virtual : public PHASIC::Virtual_ME2_Base {
    double m_eps2, m_eps, m_fin;
    double ME, MM, ME2, MM2, m_alpha;
    double m_u, m_t, m_s, m_photon_mass;
  public:
    emu_emu_Virtual(const Process_Info& pi, const Flavour_Vector& flavs,
                  const double& ep2, const double& ep);

    ~emu_emu_Virtual() {
    }

    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    Complex BornTriangle();
    Complex BornBox();
    Complex Virtual();
    // Complex MBub(double a, double b, double c);
    // Complex MTri(double a, double b, double c,
    //             double aa, double bb, double cc);
    // Complex MBox(double a, double b, double c);
    // inline Complex Conjugate(Complex a) {return conj(a);}
    // inline double Re(Complex a) {return a.real();}

  };
}

using namespace EXTRAXS;

emu_emu_Virtual::emu_emu_Virtual(const Process_Info& pi, const Flavour_Vector& flavs,
                  const double& ep2, const double& ep) :
      Virtual_ME2_Base(pi, flavs),
      m_eps2(ep2), m_eps(ep){
         Scoped_Settings s{ Settings::GetMainSettings()["YFS"] };
         m_photon_mass = s["PHOTON_MASS"].Get<double>();
         FORTRAN(ltini)();
         Setlambda(0);
         ME = Flavour(kf_e).Mass();
         MM = Flavour(kf_mu).Mass();
         ME2 = ME*ME;
         MM2 = MM*MM;
         m_mode = 0;
      }


void emu_emu_Virtual::Calc(const Vec4D_Vector& momenta) {
  double factor(1.);
  if      (m_stype&sbt::qcd) factor=2.*M_PI/AlphaQCD();
  else if (m_stype&sbt::qed) factor=2.*M_PI/AlphaQED();
  else THROW(fatal_error,"Unknown coupling.");
  // 1/epsIR
  // m_res.IR()=m_eps*factor;
  // // 1/epsIR2
  // m_res.IR2()=m_eps2.*factor;
  // // finite
  // if(m_count==1000){
    // m_count=0;
  m_s = (momenta[0]+momenta[1]).Abs2();
  m_t = (momenta[0]-momenta[2]).Abs2();
  m_u = (momenta[0]-momenta[3]).Abs2();
  m_alpha = (*aqed)(m_s);
  m_res.Finite()= 32.*pow(M_PI,3)*m_alpha*(BornBox()+BornTriangle()).real();
  Setlambda(-1);
  m_res[4]= 32.*pow(M_PI,3)*(Virtual()).real();
  // PRINT_VAR(m_res);
  // PRINT_VAR(pow(m_alpha,3)*ME*M_PI);
  // PRINT_VAR(pow(m_alpha,3));
  clearcache();// Since in general s will be different we gain nothing from cache
  // PRINT_VAR(Virtual());
  // FORTRAN(ltexi)();
}

Complex emu_emu_Virtual::Virtual(){
return (pow(m_alpha,3)*ME*M_PI*(-(MM2*(m_s + m_t - m_u)) + m_t*(-2*ME2 + m_s + m_u))*
   ((-4*ME*(C0i(cc1, ME2, m_t, ME2,0.0, ME2, ME2) + C0i(cc11, ME2, m_t, ME2,0.0, 
        ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2,0.0, ME2, ME2) + 
       C0i(cc2, ME2, m_t, ME2,0.0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2,0.0, ME2, 
        ME2)))/m_t + 4.*(MM*(D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + 
        D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
        D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
        D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
        D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
        D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
      ME*(D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
        D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
      (ME + MM)*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
        D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
        D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) - 
    2.*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))/(2.*m_t);
}
  
Complex emu_emu_Virtual::BornBox(){
return (-32.*pow(m_alpha,3)*pow(ME2,2)*M_PI*(D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
     MM2) - D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
  64.*pow(m_alpha,3)*ME*ME2*MM*M_PI*(D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
     MM2) - D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
  96.*pow(m_alpha,3)*ME2*MM2*M_PI*(D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
     MM2) - D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
  64.*pow(m_alpha,3)*ME*MM*MM2*M_PI*(D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
     MM2) - D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
  64.*pow(m_alpha,3)*pow(MM2,2)*M_PI*(D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
     MM2) - D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
  64.*pow(m_alpha,3)*ME2*M_PI*m_s*(D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
     MM2) - D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
  64.*pow(m_alpha,3)*ME*MM*M_PI*m_s*(D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
     MM2) - D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
  32.*pow(m_alpha,3)*M_PI*pow(m_s,2)*(D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
  (128.*pow(m_alpha,3)*ME*ME2*(ME - MM)*MM2*M_PI*
    (D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t + 
  (128.*pow(m_alpha,3)*ME*(ME - MM)*pow(MM2,2)*M_PI*(D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 
      0, 0, ME2, MM2) - D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))/
   m_t - (256.*pow(m_alpha,3)*ME*(ME - MM)*MM2*M_PI*m_s*
    (D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t - 
  16.*pow(m_alpha,3)*ME2*M_PI*m_t*(D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
     MM2) - D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
  16.*pow(m_alpha,3)*MM2*M_PI*m_t*(D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
     MM2) - D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
  32.*pow(m_alpha,3)*MM2*M_PI*m_u*(D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
     MM2) - D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
  (128.*pow(m_alpha,3)*ME*(ME - MM)*MM2*M_PI*m_u*(D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 
      0, 0, ME2, MM2) - D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))/
   m_t + 192.*pow(m_alpha,3)*pow(ME2,2)*M_PI*
   (-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
  64.*pow(m_alpha,3)*ME*ME2*MM*M_PI*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
      MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
  32.*pow(m_alpha,3)*ME2*MM2*M_PI*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
      MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
  64.*pow(m_alpha,3)*ME*MM*MM2*M_PI*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
      MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
  32.*pow(m_alpha,3)*pow(MM2,2)*M_PI*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
      MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
  128.*pow(m_alpha,3)*ME2*M_PI*m_s*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
      MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
  64.*pow(m_alpha,3)*ME*MM*M_PI*m_s*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
      MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
  64.*pow(m_alpha,3)*MM2*M_PI*m_s*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
      MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
  64.*pow(m_alpha,3)*M_PI*pow(m_s,2)*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
      MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
  (64.*pow(m_alpha,3)*pow(ME2,2)*MM2*M_PI*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, 
       ME2, MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t - 
  (640*pow(m_alpha,3)*ME2*pow(MM2,2)*M_PI*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, 
       ME2, MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t + 
  (64.*pow(m_alpha,3)*pow(MM2,3)*M_PI*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t - 
  (128.*pow(m_alpha,3)*ME2*MM2*M_PI*m_s*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, 
       ME2, MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t - 
  (128.*pow(m_alpha,3)*pow(MM2,2)*M_PI*m_s*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t + 
  (64.*pow(m_alpha,3)*MM2*M_PI*pow(m_s,2)*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t + 
  16.*pow(m_alpha,3)*ME2*M_PI*m_t*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
      MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
  16.*pow(m_alpha,3)*MM2*M_PI*m_t*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
      MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
  (256.*pow(m_alpha,3)*ME2*M_PI*(4.*pow(MM2,2) - 4.*MM2*m_t + pow(m_t,2))*
    (-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t + 
  (16.*pow(m_alpha,3)*ME2*M_PI*(12.*pow(MM2,2) - 4.*MM2*m_t + pow(m_t,2))*
    (-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t - 
  256.*pow(m_alpha,3)*ME2*M_PI*m_u*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
      MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
  32.*pow(m_alpha,3)*MM2*M_PI*m_u*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
      MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
  128.*pow(m_alpha,3)*M_PI*pow(m_u,2)*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
      MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
  (32.*pow(m_alpha,3)*ME*MM*M_PI*(-6.*(ME2*MM2 + pow(MM2,2)) + (5*ME2 + MM2)*m_t + 
     (6.*MM2 - 5*m_t)*m_u)*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t + 
  (256.*pow(m_alpha,3)*MM2*M_PI*(pow(ME2,2) + 2.*ME2*MM2 + pow(MM2,2) - 2.*(ME2 + MM2)*m_u + pow(m_u,2))*
    (-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t + 
  (16.*pow(m_alpha,3)*MM2*M_PI*(pow(m_t,2) - 2.*(ME2 + MM2)*(m_t + 4.*m_u) + 
     4.*(pow(ME2,2) + ME2*MM2 + pow(MM2,2) + m_t*m_u + pow(m_u,2)))*
    (-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t + 
  (96.*pow(m_alpha,3)*ME*ME2*MM2*M_PI*
    (-2.*(MM*(D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       ME*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) - 
     2.*(ME + MM)*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t + 
  (96.*pow(m_alpha,3)*ME*pow(MM2,2)*M_PI*
    (-2.*(MM*(D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       ME*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) - 
     2.*(ME + MM)*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t - 
  (32.*pow(m_alpha,3)*ME*MM2*M_PI*m_s*
    (-2.*(MM*(D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       ME*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) - 
     2.*(ME + MM)*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t + 
  (32.*pow(m_alpha,3)*ME*M_PI*(-2.*(ME2*MM2 + pow(MM2,2)) + 2.*MM2*m_s + (ME2 + MM2 - m_s)*m_t)*
    (-2.*(MM*(D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       ME*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) - 
     2.*(ME + MM)*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t + 
  4.*pow(m_alpha,2)*ME*((8.*m_alpha*MM2*M_PI*(ME2 + MM2 - m_s))/m_t + 
    (8.*m_alpha*M_PI*(-2.*(ME2*MM2 + pow(MM2,2)) + 2.*MM2*m_s + (ME2 + MM2 - m_s)*m_t))/m_t - 
    (8.*m_alpha*M_PI*(pow(MM2,2) + (MM2 - m_t)*(ME2 - m_u)))/m_t + 
    (16.*m_alpha*MM2*M_PI*(ME2 + MM2 - m_u))/m_t)*
   (-2.*(MM*(D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
        D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
        D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
        D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
        D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
        D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
        D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
        D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
      ME*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
        D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
        D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
        D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) - 
    2.*(ME + MM)*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
  (64.*pow(m_alpha,3)*ME*MM2*M_PI*m_u*
    (-2.*(MM*(D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       ME*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) - 
     2.*(ME + MM)*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t + 
  16.*pow(m_alpha,3)*MM2*M_PI*(8.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 
    4.*(ME2 + MM2 - m_s)*D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    8.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    4.*(2.*ME2 + ME*MM - 2.*(ME2 - MM2) - 4.*MM2)*D0i(dd13, m_t, ME2, m_s, MM2, ME2, 
      MM2, 0, 0, ME2, MM2) + 4.*ME*MM*(D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
       0, ME2, MM2) + D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
        D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
      D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
    4.*(ME*MM + MM2 - m_s - m_t)*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
      MM2) + (4.*ME*MM - 6.*MM2 + 3.*m_t)*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
      0, ME2, MM2) - (4.*ME2 - 2.*ME*MM + 3.*m_t - 4.*m_u)*D0i(dd23, m_t, ME2, m_u, MM2, 
      ME2, MM2, 0, 0, ME2, MM2) + 
    4.*((ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
        MM2) - (MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
        MM2) - MM2*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*ME2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      ME2*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      ME2*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      ME2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      (ME2 + MM2 - m_s - m_t)*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
        MM2))) + (16.*pow(m_alpha,3)*pow(ME2,2)*M_PI*
    (8.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 4.*(ME2 + MM2 - m_s)*
      D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     8.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*(2.*ME2 + ME*MM - 2.*(ME2 - MM2) - 4.*MM2)*D0i(dd13, m_t, ME2, m_s, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 4.*ME*MM*(D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
        0, ME2, MM2) + D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     4.*(ME*MM + MM2 - m_s - m_t)*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) + (4.*ME*MM - 6.*MM2 + 3.*m_t)*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
       0, ME2, MM2) - (4.*ME2 - 2.*ME*MM + 3.*m_t - 4.*m_u)*D0i(dd23, m_t, ME2, m_u, MM2, 
       ME2, MM2, 0, 0, ME2, MM2) + 
     4.*((ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - (MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - MM2*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME2*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME2*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (ME2 + MM2 - m_s - m_t)*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2))))/m_t + (64.*pow(m_alpha,3)*ME2*MM2*M_PI*
    (8.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 4.*(ME2 + MM2 - m_s)*
      D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     8.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*(2.*ME2 + ME*MM - 2.*(ME2 - MM2) - 4.*MM2)*D0i(dd13, m_t, ME2, m_s, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 4.*ME*MM*(D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
        0, ME2, MM2) + D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     4.*(ME*MM + MM2 - m_s - m_t)*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) + (4.*ME*MM - 6.*MM2 + 3.*m_t)*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
       0, ME2, MM2) - (4.*ME2 - 2.*ME*MM + 3.*m_t - 4.*m_u)*D0i(dd23, m_t, ME2, m_u, MM2, 
       ME2, MM2, 0, 0, ME2, MM2) + 
     4.*((ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - (MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - MM2*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME2*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME2*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (ME2 + MM2 - m_s - m_t)*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2))))/m_t + (16.*pow(m_alpha,3)*pow(MM2,2)*M_PI*
    (8.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 4.*(ME2 + MM2 - m_s)*
      D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     8.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*(2.*ME2 + ME*MM - 2.*(ME2 - MM2) - 4.*MM2)*D0i(dd13, m_t, ME2, m_s, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 4.*ME*MM*(D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
        0, ME2, MM2) + D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     4.*(ME*MM + MM2 - m_s - m_t)*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) + (4.*ME*MM - 6.*MM2 + 3.*m_t)*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
       0, ME2, MM2) - (4.*ME2 - 2.*ME*MM + 3.*m_t - 4.*m_u)*D0i(dd23, m_t, ME2, m_u, MM2, 
       ME2, MM2, 0, 0, ME2, MM2) + 
     4.*((ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - (MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - MM2*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME2*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME2*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (ME2 + MM2 - m_s - m_t)*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2))))/m_t - (32.*pow(m_alpha,3)*ME2*M_PI*m_s*
    (8.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 4.*(ME2 + MM2 - m_s)*
      D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     8.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*(2.*ME2 + ME*MM - 2.*(ME2 - MM2) - 4.*MM2)*D0i(dd13, m_t, ME2, m_s, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 4.*ME*MM*(D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
        0, ME2, MM2) + D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     4.*(ME*MM + MM2 - m_s - m_t)*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) + (4.*ME*MM - 6.*MM2 + 3.*m_t)*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
       0, ME2, MM2) - (4.*ME2 - 2.*ME*MM + 3.*m_t - 4.*m_u)*D0i(dd23, m_t, ME2, m_u, MM2, 
       ME2, MM2, 0, 0, ME2, MM2) + 
     4.*((ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - (MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - MM2*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME2*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME2*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (ME2 + MM2 - m_s - m_t)*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2))))/m_t - (32.*pow(m_alpha,3)*MM2*M_PI*m_s*
    (8.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 4.*(ME2 + MM2 - m_s)*
      D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     8.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*(2.*ME2 + ME*MM - 2.*(ME2 - MM2) - 4.*MM2)*D0i(dd13, m_t, ME2, m_s, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 4.*ME*MM*(D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
        0, ME2, MM2) + D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     4.*(ME*MM + MM2 - m_s - m_t)*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) + (4.*ME*MM - 6.*MM2 + 3.*m_t)*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
       0, ME2, MM2) - (4.*ME2 - 2.*ME*MM + 3.*m_t - 4.*m_u)*D0i(dd23, m_t, ME2, m_u, MM2, 
       ME2, MM2, 0, 0, ME2, MM2) + 
     4.*((ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - (MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - MM2*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME2*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME2*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (ME2 + MM2 - m_s - m_t)*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2))))/m_t + (16.*pow(m_alpha,3)*M_PI*pow(m_s,2)*
    (8.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 4.*(ME2 + MM2 - m_s)*
      D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     8.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*(2.*ME2 + ME*MM - 2.*(ME2 - MM2) - 4.*MM2)*D0i(dd13, m_t, ME2, m_s, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 4.*ME*MM*(D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
        0, ME2, MM2) + D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     4.*(ME*MM + MM2 - m_s - m_t)*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) + (4.*ME*MM - 6.*MM2 + 3.*m_t)*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
       0, ME2, MM2) - (4.*ME2 - 2.*ME*MM + 3.*m_t - 4.*m_u)*D0i(dd23, m_t, ME2, m_u, MM2, 
       ME2, MM2, 0, 0, ME2, MM2) + 
     4.*((ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - (MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - MM2*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME2*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME2*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (ME2 + MM2 - m_s - m_t)*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2))))/m_t - 2.*pow(m_alpha,2)*((-32.*m_alpha*ME2*MM2*M_PI)/m_t - 
    (8.*m_alpha*M_PI*pow(ME2 + MM2 - m_s,2))/m_t - (8.*m_alpha*MM2*M_PI*(-2.*ME2 + m_t))/m_t - 
    (8.*m_alpha*ME2*M_PI*(-2.*MM2 + m_t))/m_t)*(8.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 
    4.*(ME2 + MM2 - m_s)*D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    8.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    4.*(2.*ME2 + ME*MM - 2.*(ME2 - MM2) - 4.*MM2)*D0i(dd13, m_t, ME2, m_s, MM2, ME2, 
      MM2, 0, 0, ME2, MM2) + 4.*ME*MM*(D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
       0, ME2, MM2) + D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
        D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
      D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
    4.*(ME*MM + MM2 - m_s - m_t)*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
      MM2) + (4.*ME*MM - 6.*MM2 + 3.*m_t)*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
      0, ME2, MM2) - (4.*ME2 - 2.*ME*MM + 3.*m_t - 4.*m_u)*D0i(dd23, m_t, ME2, m_u, MM2, 
      ME2, MM2, 0, 0, ME2, MM2) + 
    4.*((ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
        MM2) - (MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
        MM2) - MM2*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*ME2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      ME2*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      ME2*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      ME2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      (ME2 + MM2 - m_s - m_t)*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
        MM2))) - 128.*pow(m_alpha,3)*ME2*MM*M_PI*
   (-2.*ME*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
    MM*(2.*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      6.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      2.*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      4.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) + 
  (128.*pow(m_alpha,3)*ME2*MM*MM2*M_PI*
    (-2.*ME*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     MM*(2.*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       6.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       4.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t + 
  32.*pow(m_alpha,3)*MM*M_PI*m_u*(-2.*ME*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, 
       ME2, MM2) + D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
    MM*(2.*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      6.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      2.*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      4.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) - 
  (64.*pow(m_alpha,3)*MM*M_PI*(pow(ME2,2) + 2.*ME2*MM2 + pow(MM2,2) + m_s*m_u - 
     (ME2 + MM2)*(m_s + m_u))*(-2.*ME*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, 
        ME2, MM2) + D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     MM*(2.*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       6.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       4.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t - 
  (160*pow(m_alpha,3)*ME*ME2*MM2*M_PI*
    (-2.*ME*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     MM*(2.*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t - 
  (160*pow(m_alpha,3)*ME*pow(MM2,2)*M_PI*
    (-2.*ME*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     MM*(2.*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t + 
  (32.*pow(m_alpha,3)*ME*MM2*M_PI*m_s*
    (-2.*ME*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     MM*(2.*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t - 
  (64.*pow(m_alpha,3)*ME*M_PI*(-2.*(ME2*MM2 + pow(MM2,2)) + 2.*MM2*m_s + (ME2 + MM2 - m_s)*m_t)*
    (-2.*ME*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     MM*(2.*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t + 
  (32.*pow(m_alpha,3)*ME*M_PI*(pow(MM2,2) + (MM2 - m_t)*(ME2 - m_u))*
    (-2.*ME*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     MM*(2.*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t + 
  (128.*pow(m_alpha,3)*ME*MM2*M_PI*m_u*
    (-2.*ME*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     MM*(2.*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t - 
  64.*pow(m_alpha,3)*ME2*MM*M_PI*(2.*(ME - MM)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
      0, ME2, MM2) + 2.*(ME - MM)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
      ME2, MM2) + 2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    2.*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    2.*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    4.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    2.*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    4.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    2.*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    4.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*ME*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    4.*ME*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    6.*ME*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    4.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    3.*ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    5*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    2.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    4.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
  16.*pow(m_alpha,3)*MM*MM2*M_PI*(2.*(ME - MM)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
      0, ME2, MM2) + 2.*(ME - MM)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
      ME2, MM2) + 2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    2.*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    2.*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    4.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    2.*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    4.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    2.*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    4.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*ME*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    4.*ME*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    6.*ME*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    4.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    3.*ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    5*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    2.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    4.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
  (64.*pow(m_alpha,3)*ME*ME2*MM2*M_PI*(2.*(ME - MM)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 2.*(ME - MM)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 
       0, 0, ME2, MM2) + 2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) - 2.*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     4.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     6.*ME*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     3.*ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     5*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t + 
  (32.*pow(m_alpha,3)*ME2*MM*MM2*M_PI*(2.*(ME - MM)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 2.*(ME - MM)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 
       0, 0, ME2, MM2) + 2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) - 2.*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     4.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     6.*ME*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     3.*ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     5*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t - 
  (64.*pow(m_alpha,3)*ME*pow(MM2,2)*M_PI*(2.*(ME - MM)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 
       0, 0, ME2, MM2) + 2.*(ME - MM)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
       ME2, MM2) + 2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     4.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     6.*ME*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     3.*ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     5*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t + 
  (128.*pow(m_alpha,3)*ME*MM2*M_PI*m_s*(2.*(ME - MM)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 
       0, 0, ME2, MM2) + 2.*(ME - MM)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
       ME2, MM2) + 2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     4.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     6.*ME*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     3.*ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     5*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t - 
  (32.*pow(m_alpha,3)*MM*M_PI*(pow(ME2,2) + 2.*ME2*MM2 + pow(MM2,2) - 2.*(ME2 + MM2)*m_s + pow(m_s,2))*
    (2.*(ME - MM)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*(ME - MM)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     4.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     6.*ME*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     3.*ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     5*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t - 
  (64.*pow(m_alpha,3)*ME*MM2*M_PI*m_u*(2.*(ME - MM)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 
       0, 0, ME2, MM2) + 2.*(ME - MM)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
       ME2, MM2) + 2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     4.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     6.*ME*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     3.*ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     5*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t + 
  (16.*pow(m_alpha,3)*MM*M_PI*(-2.*(pow(ME2,2) + ME2*(2.*MM2 - m_s)) + MM2*m_s + 
     (2.*ME2 + MM2 - 2.*m_s)*m_u)*(2.*(ME - MM)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 
       0, 0, ME2, MM2) + 2.*(ME - MM)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
       ME2, MM2) + 2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     4.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     6.*ME*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     3.*ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     5*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t + 
  (32.*pow(m_alpha,3)*MM*M_PI*(pow(ME2,2) + pow(MM2,2) + 2.*(ME2 + MM2)*m_s - pow(m_s,2) - 
     2.*(ME2*MM2 + m_s*m_u))*(2.*(ME - MM)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, 
       ME2, MM2) + 2.*(ME - MM)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) + 2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     4.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     6.*ME*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     3.*ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     5*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t - 
  64.*pow(m_alpha,3)*ME*ME2*M_PI*
   (2.*ME*(D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
    MM*(2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      2.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      2.*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) + 
  64.*pow(m_alpha,3)*ME*MM2*M_PI*
   (2.*ME*(D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
    MM*(2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      2.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      2.*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) + 
  64.*pow(m_alpha,3)*ME*M_PI*m_s*(2.*ME*(D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) + D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
    MM*(2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      2.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      2.*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) - 
  (96.*pow(m_alpha,3)*ME*ME2*MM2*M_PI*
    (2.*ME*(D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     MM*(2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t - 
  (96.*pow(m_alpha,3)*ME*pow(MM2,2)*M_PI*
    (2.*ME*(D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     MM*(2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t + 
  (96.*pow(m_alpha,3)*ME*MM2*M_PI*m_s*
    (2.*ME*(D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     MM*(2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t - 
  (32.*pow(m_alpha,3)*ME*M_PI*(-3.*(ME2*MM2 + pow(MM2,2)) + (ME2 + 2.*MM2)*m_t + 
     (3.*MM2 - m_t)*m_u)*(2.*ME*(D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
        MM2) + D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     MM*(2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t + 
  128.*pow(m_alpha,3)*ME2*MM*M_PI*
   (ME*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
    MM*(2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      2.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) + 
  32.*pow(m_alpha,3)*MM*MM2*M_PI*
   (ME*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
    MM*(2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      2.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) - 
  (64.*pow(m_alpha,3)*ME2*MM*MM2*M_PI*
    (ME*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     MM*(2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t + 
  (64.*pow(m_alpha,3)*MM*M_PI*(pow(ME2,2) + 2.*ME2*MM2 + pow(MM2,2) - 2.*(ME2 + MM2)*m_s + pow(m_s,2))*
    (ME*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     MM*(2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t - 
  (32.*pow(m_alpha,3)*MM*M_PI*(-2.*(pow(ME2,2) + ME2*(2.*MM2 - m_s)) + MM2*m_s + 
     (2.*ME2 + MM2 - 2.*m_s)*m_u)*
    (ME*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     MM*(2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t - 
  (64.*pow(m_alpha,3)*MM*M_PI*(pow(ME2,2) + pow(MM2,2) + 2.*(ME2 + MM2)*m_s - pow(m_s,2) - 
     2.*(ME2*MM2 + m_s*m_u))*
    (ME*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     MM*(2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t - 
  128.*pow(m_alpha,3)*ME2*MM*M_PI*(2.*ME*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, 
      ME2, MM2) + 2.*ME*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    6.*MM*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    6.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    6.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    6.*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    4.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    4.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    6.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    2.*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
  (128.*pow(m_alpha,3)*ME2*MM*MM2*M_PI*(2.*ME*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
       0, ME2, MM2) + 2.*ME*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) + 2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     6.*MM*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     6.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     6.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     6.*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     4.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     6.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t + 
  32.*pow(m_alpha,3)*MM*M_PI*m_u*(2.*ME*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
      MM2) + 2.*ME*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    6.*MM*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    6.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    6.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    6.*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    4.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    4.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    6.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    2.*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
  (64.*pow(m_alpha,3)*MM*M_PI*(pow(ME2,2) + 2.*ME2*MM2 + pow(MM2,2) + m_s*m_u - 
     (ME2 + MM2)*(m_s + m_u))*(2.*ME*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) + 2.*ME*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     6.*MM*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     6.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     6.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     6.*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     4.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     6.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))/m_t - 
  (160*pow(m_alpha,3)*ME*ME2*MM*M_PI*(2.*(ME2 + 5*MM2 - m_s - 2.*m_t)*
      (D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     2.*(2.*ME2 - ME*MM + 5*MM2 - m_s - 2.*m_t)*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 
        0, 0, ME2, MM2) + D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
        MM2)) + 2.*(ME2 - ME*MM + 4.*MM2 - m_s - 2.*m_t)*D0i(dd13, m_t, ME2, m_s, MM2, 
       ME2, MM2, 0, 0, ME2, MM2) + 2.*(ME2 + ME*MM + 6.*MM2 - m_s - 2.*m_t)*
      D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*MM*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     (2.*ME2 - 4.*ME*MM + 4.*MM2 - 2.*m_s - 5*m_t)*D0i(dd23, m_t, ME2, m_s, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + (3.*(ME2 + MM2) + 2.*m_t - m_u)*
      D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*((ME2 - 2.*ME*MM)*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME2*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME2*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       4.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t - 
  (160*pow(m_alpha,3)*ME*MM*MM2*M_PI*(2.*(ME2 + 5*MM2 - m_s - 2.*m_t)*
      (D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     2.*(2.*ME2 - ME*MM + 5*MM2 - m_s - 2.*m_t)*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 
        0, 0, ME2, MM2) + D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
        MM2)) + 2.*(ME2 - ME*MM + 4.*MM2 - m_s - 2.*m_t)*D0i(dd13, m_t, ME2, m_s, MM2, 
       ME2, MM2, 0, 0, ME2, MM2) + 2.*(ME2 + ME*MM + 6.*MM2 - m_s - 2.*m_t)*
      D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*MM*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     (2.*ME2 - 4.*ME*MM + 4.*MM2 - 2.*m_s - 5*m_t)*D0i(dd23, m_t, ME2, m_s, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + (3.*(ME2 + MM2) + 2.*m_t - m_u)*
      D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*((ME2 - 2.*ME*MM)*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME2*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME2*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       4.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t + 
  (128.*pow(m_alpha,3)*ME*MM*M_PI*m_s*(2.*(ME2 + 5*MM2 - m_s - 2.*m_t)*
      (D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     2.*(2.*ME2 - ME*MM + 5*MM2 - m_s - 2.*m_t)*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 
        0, 0, ME2, MM2) + D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
        MM2)) + 2.*(ME2 - ME*MM + 4.*MM2 - m_s - 2.*m_t)*D0i(dd13, m_t, ME2, m_s, MM2, 
       ME2, MM2, 0, 0, ME2, MM2) + 2.*(ME2 + ME*MM + 6.*MM2 - m_s - 2.*m_t)*
      D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*MM*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     (2.*ME2 - 4.*ME*MM + 4.*MM2 - 2.*m_s - 5*m_t)*D0i(dd23, m_t, ME2, m_s, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + (3.*(ME2 + MM2) + 2.*m_t - m_u)*
      D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*((ME2 - 2.*ME*MM)*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME2*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME2*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       4.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t + 
  (32.*pow(m_alpha,3)*ME*MM*M_PI*m_u*(2.*(ME2 + 5*MM2 - m_s - 2.*m_t)*
      (D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     2.*(2.*ME2 - ME*MM + 5*MM2 - m_s - 2.*m_t)*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 
        0, 0, ME2, MM2) + D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
        MM2)) + 2.*(ME2 - ME*MM + 4.*MM2 - m_s - 2.*m_t)*D0i(dd13, m_t, ME2, m_s, MM2, 
       ME2, MM2, 0, 0, ME2, MM2) + 2.*(ME2 + ME*MM + 6.*MM2 - m_s - 2.*m_t)*
      D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     4.*ME*MM*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     (2.*ME2 - 4.*ME*MM + 4.*MM2 - 2.*m_s - 5*m_t)*D0i(dd23, m_t, ME2, m_s, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + (3.*(ME2 + MM2) + 2.*m_t - m_u)*
      D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     2.*MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*((ME2 - 2.*ME*MM)*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME2*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME2*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       4.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t + 
  (96.*pow(m_alpha,3)*ME*ME2*MM*M_PI*(2.*(ME2 + MM2 - m_s)*D0i(dd1, m_t, ME2, m_s, MM2, 
       ME2, MM2, 0, 0, ME2, MM2) + 2.*(ME2 + MM2 - m_s)*D0i(dd1, m_t, ME2, m_u, MM2, 
       ME2, MM2, 0, 0, ME2, MM2) + 2.*(2.*ME2 - ME*MM + MM2 - m_s)*
      (D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     2.*(ME2 + ME*MM + 2.*MM2 - m_s)*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
       ME2, MM2) + (2.*ME2 - 4.*ME*MM + 4.*MM2 + 3.*(ME2 + MM2) - m_s - 2.*m_u)*
      D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*((ME2 + MM2 - m_s)*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (ME2 + MM2 - m_s)*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME2*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME*MM*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME2*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME2*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME2*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME2*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*MM2*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       m_u*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t + 
  (96.*pow(m_alpha,3)*ME*MM*MM2*M_PI*(2.*(ME2 + MM2 - m_s)*D0i(dd1, m_t, ME2, m_s, MM2, 
       ME2, MM2, 0, 0, ME2, MM2) + 2.*(ME2 + MM2 - m_s)*D0i(dd1, m_t, ME2, m_u, MM2, 
       ME2, MM2, 0, 0, ME2, MM2) + 2.*(2.*ME2 - ME*MM + MM2 - m_s)*
      (D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     2.*(ME2 + ME*MM + 2.*MM2 - m_s)*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
       ME2, MM2) + (2.*ME2 - 4.*ME*MM + 4.*MM2 + 3.*(ME2 + MM2) - m_s - 2.*m_u)*
      D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*((ME2 + MM2 - m_s)*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (ME2 + MM2 - m_s)*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME2*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME*MM*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME2*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME2*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME2*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME2*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*MM2*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       m_u*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t - 
  (96.*pow(m_alpha,3)*ME*MM*M_PI*m_u*(2.*(ME2 + MM2 - m_s)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 2.*(ME2 + MM2 - m_s)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 2.*(2.*ME2 - ME*MM + MM2 - m_s)*
      (D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     2.*(ME2 + ME*MM + 2.*MM2 - m_s)*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
       ME2, MM2) + (2.*ME2 - 4.*ME*MM + 4.*MM2 + 3.*(ME2 + MM2) - m_s - 2.*m_u)*
      D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*((ME2 + MM2 - m_s)*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (ME2 + MM2 - m_s)*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME2*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME*MM*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME2*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME2*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME2*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME2*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*MM2*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       m_u*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t - 
  (64.*pow(m_alpha,3)*ME*ME2*MM*M_PI*(2.*(ME*MM - MM2)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 2.*(ME*MM - MM2)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 2.*ME2*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
        0, ME2, MM2) + D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     (2.*ME2 + 5*ME*MM - 5*MM2)*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) + 2.*((ME*MM - MM2)*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + (ME*MM - MM2)*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + ME*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*MM2*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME*MM*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*MM2*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t - 
  (64.*pow(m_alpha,3)*ME*MM*MM2*M_PI*(2.*(ME*MM - MM2)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 2.*(ME*MM - MM2)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 2.*ME2*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
        0, ME2, MM2) + D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     (2.*ME2 + 5*ME*MM - 5*MM2)*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) + 2.*((ME*MM - MM2)*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + (ME*MM - MM2)*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + ME*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*MM2*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME*MM*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*MM2*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t - 
  (64.*pow(m_alpha,3)*ME*MM*M_PI*m_s*(2.*(ME*MM - MM2)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 2.*(ME*MM - MM2)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 2.*ME2*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
        0, ME2, MM2) + D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     (2.*ME2 + 5*ME*MM - 5*MM2)*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) + 2.*((ME*MM - MM2)*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + (ME*MM - MM2)*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + ME*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*MM2*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME*MM*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*MM2*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t + 
  (128.*pow(m_alpha,3)*ME*MM*M_PI*m_u*(2.*(ME*MM - MM2)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 2.*(ME*MM - MM2)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 2.*ME2*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
        0, ME2, MM2) + D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     (2.*ME2 + 5*ME*MM - 5*MM2)*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) + 2.*((ME*MM - MM2)*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + (ME*MM - MM2)*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + ME*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*MM2*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME*MM*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*MM2*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t + 
  64.*pow(m_alpha,3)*ME2*M_PI*(4.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 
    2.*(ME2 + MM2 - m_s)*D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*(ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
      MM2) - 12.*D0i(dd00, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    16.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*(ME*MM + 2.*(MM2 - m_t))*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
      MM2) + 2.*(ME*MM + 2.*MM2 - 2.*(MM2 + m_t))*D0i(dd13, m_t, ME2, m_u, MM2, ME2, 
      MM2, 0, 0, ME2, MM2) + 2.*(ME2 + ME*MM - MM2 + m_t - m_u)*
     D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    (4.*ME2 - 3.*ME*MM + MM2 - 4.*(2.*ME2 + m_t - m_u))*D0i(dd23, m_t, ME2, m_u, MM2, 
      ME2, MM2, 0, 0, ME2, MM2) + 
    2.*ME2*(2.*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
    2.*((ME*MM + MM2)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      (ME*MM + MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
        MM2) + ME*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      2.*m_t*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      ME*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      2.*m_t*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      2.*m_t*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      2.*m_t*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      ME*MM*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      MM2*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      m_s*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      m_t*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      m_s*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      m_t*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      3.*MM2*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))) + 
  64.*pow(m_alpha,3)*MM2*M_PI*(4.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 
    2.*(ME2 + MM2 - m_s)*D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*(ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
      MM2) - 12.*D0i(dd00, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
    16.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    2.*(ME*MM + 2.*(MM2 - m_t))*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
      MM2) + 2.*(ME*MM + 2.*MM2 - 2.*(MM2 + m_t))*D0i(dd13, m_t, ME2, m_u, MM2, ME2, 
      MM2, 0, 0, ME2, MM2) + 2.*(ME2 + ME*MM - MM2 + m_t - m_u)*
     D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
    (4.*ME2 - 3.*ME*MM + MM2 - 4.*(2.*ME2 + m_t - m_u))*D0i(dd23, m_t, ME2, m_u, MM2, 
      ME2, MM2, 0, 0, ME2, MM2) + 
    2.*ME2*(2.*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
    2.*((ME*MM + MM2)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      (ME*MM + MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
        MM2) + ME*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      2.*m_t*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      ME*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      2.*m_t*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      2.*m_t*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      2.*m_t*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      ME*MM*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      MM2*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      m_s*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      m_t*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
      2.*MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      m_s*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      m_t*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
      3.*MM2*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))) + 
  (32.*pow(m_alpha,3)*pow(ME2,2)*M_PI*(4.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 
     2.*(ME2 + MM2 - m_s)*D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*(ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) - 12.*D0i(dd00, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     16.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*(ME*MM + 2.*(MM2 - m_t))*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) + 2.*(ME*MM + 2.*MM2 - 2.*(MM2 + m_t))*D0i(dd13, m_t, ME2, m_u, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 2.*(ME2 + ME*MM - MM2 + m_t - m_u)*
      D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     (4.*ME2 - 3.*ME*MM + MM2 - 4.*(2.*ME2 + m_t - m_u))*D0i(dd23, m_t, ME2, m_u, MM2, 
       ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME2*(2.*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     2.*((ME*MM + MM2)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (ME*MM + MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + ME*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME*MM*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_t*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_t*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*MM2*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t - 
  (64.*pow(m_alpha,3)*ME2*MM2*M_PI*(4.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 
     2.*(ME2 + MM2 - m_s)*D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*(ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) - 12.*D0i(dd00, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     16.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*(ME*MM + 2.*(MM2 - m_t))*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) + 2.*(ME*MM + 2.*MM2 - 2.*(MM2 + m_t))*D0i(dd13, m_t, ME2, m_u, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 2.*(ME2 + ME*MM - MM2 + m_t - m_u)*
      D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     (4.*ME2 - 3.*ME*MM + MM2 - 4.*(2.*ME2 + m_t - m_u))*D0i(dd23, m_t, ME2, m_u, MM2, 
       ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME2*(2.*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     2.*((ME*MM + MM2)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (ME*MM + MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + ME*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME*MM*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_t*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_t*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*MM2*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t + 
  (32.*pow(m_alpha,3)*pow(MM2,2)*M_PI*(4.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 
     2.*(ME2 + MM2 - m_s)*D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*(ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) - 12.*D0i(dd00, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     16.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*(ME*MM + 2.*(MM2 - m_t))*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) + 2.*(ME*MM + 2.*MM2 - 2.*(MM2 + m_t))*D0i(dd13, m_t, ME2, m_u, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 2.*(ME2 + ME*MM - MM2 + m_t - m_u)*
      D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     (4.*ME2 - 3.*ME*MM + MM2 - 4.*(2.*ME2 + m_t - m_u))*D0i(dd23, m_t, ME2, m_u, MM2, 
       ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME2*(2.*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     2.*((ME*MM + MM2)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (ME*MM + MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + ME*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME*MM*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_t*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_t*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*MM2*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t - 
  (64.*pow(m_alpha,3)*ME2*M_PI*m_u*(4.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 
     2.*(ME2 + MM2 - m_s)*D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*(ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) - 12.*D0i(dd00, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     16.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*(ME*MM + 2.*(MM2 - m_t))*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) + 2.*(ME*MM + 2.*MM2 - 2.*(MM2 + m_t))*D0i(dd13, m_t, ME2, m_u, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 2.*(ME2 + ME*MM - MM2 + m_t - m_u)*
      D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     (4.*ME2 - 3.*ME*MM + MM2 - 4.*(2.*ME2 + m_t - m_u))*D0i(dd23, m_t, ME2, m_u, MM2, 
       ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME2*(2.*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     2.*((ME*MM + MM2)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (ME*MM + MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + ME*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME*MM*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_t*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_t*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*MM2*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t - 
  (64.*pow(m_alpha,3)*MM2*M_PI*m_u*(4.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 
     2.*(ME2 + MM2 - m_s)*D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*(ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) - 12.*D0i(dd00, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     16.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*(ME*MM + 2.*(MM2 - m_t))*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) + 2.*(ME*MM + 2.*MM2 - 2.*(MM2 + m_t))*D0i(dd13, m_t, ME2, m_u, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 2.*(ME2 + ME*MM - MM2 + m_t - m_u)*
      D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     (4.*ME2 - 3.*ME*MM + MM2 - 4.*(2.*ME2 + m_t - m_u))*D0i(dd23, m_t, ME2, m_u, MM2, 
       ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME2*(2.*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     2.*((ME*MM + MM2)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (ME*MM + MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + ME*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME*MM*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_t*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_t*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*MM2*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t + 
  (32.*pow(m_alpha,3)*M_PI*pow(m_u,2)*(4.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 
     2.*(ME2 + MM2 - m_s)*D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*(ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) - 12.*D0i(dd00, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
     16.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     2.*(ME*MM + 2.*(MM2 - m_t))*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
       MM2) + 2.*(ME*MM + 2.*MM2 - 2.*(MM2 + m_t))*D0i(dd13, m_t, ME2, m_u, MM2, ME2, 
       MM2, 0, 0, ME2, MM2) + 2.*(ME2 + ME*MM - MM2 + m_t - m_u)*
      D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
     (4.*ME2 - 3.*ME*MM + MM2 - 4.*(2.*ME2 + m_t - m_u))*D0i(dd23, m_t, ME2, m_u, MM2, 
       ME2, MM2, 0, 0, ME2, MM2) + 
     2.*ME2*(2.*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     2.*((ME*MM + MM2)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (ME*MM + MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + ME*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*m_t*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       ME*MM*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM2*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_t*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_s*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       m_t*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       3.*MM2*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))))/m_t + 
  (8.*m_alpha*M_PI*(-16.*pow(m_alpha,2)*ME*ME2*(ME - MM)*MM2*
      (D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
     16.*pow(m_alpha,2)*ME*(ME - MM)*pow(MM2,2)*(D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, 
        ME2, MM2) - D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     12.*pow(m_alpha,2)*ME2*MM2*m_t*(D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
        MM2) - D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     4.*pow(m_alpha,2)*pow(MM2,2)*m_t*(D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
     2.*pow(m_alpha,2)*ME2*pow(m_t,2)*(D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
     2.*pow(m_alpha,2)*MM2*pow(m_t,2)*(D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     16.*pow(m_alpha,2)*ME*(ME - MM)*MM2*m_u*(D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, 
        ME2, MM2) - D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
     4.*pow(m_alpha,2)*MM2*m_t*m_u*(D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
     64.*pow(m_alpha,2)*ME2*pow(MM2,2)*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     16.*pow(m_alpha,2)*pow(ME2,2)*m_t*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     4.*pow(m_alpha,2)*ME2*MM2*m_t*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
     4.*pow(m_alpha,2)*pow(MM2,2)*m_t*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     2.*pow(m_alpha,2)*ME2*pow(m_t,2)*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     2.*pow(m_alpha,2)*MM2*pow(m_t,2)*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     2.*pow(m_alpha,2)*ME2*(12.*pow(MM2,2) - 4.*MM2*m_t + pow(m_t,2))*
      (-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
     32.*pow(m_alpha,2)*ME2*m_t*m_u*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     4.*pow(m_alpha,2)*MM2*m_t*m_u*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     16.*pow(m_alpha,2)*m_t*pow(m_u,2)*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
     4.*pow(m_alpha,2)*ME*MM*(-6.*(ME2*MM2 + pow(MM2,2)) + (5*ME2 + MM2)*m_t + (6.*MM2 - 5*m_t)*m_u)*
      (-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     2.*pow(m_alpha,2)*MM2*(pow(m_t,2) - 2.*(ME2 + MM2)*(m_t + 4.*m_u) + 
       4.*(pow(ME2,2) + ME2*MM2 + pow(MM2,2) + m_t*m_u + pow(m_u,2)))*
      (-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
     4.*pow(m_alpha,2)*ME*(pow(MM2,2) + (MM2 - m_t)*(ME2 - m_u))*
      (-2.*(MM*(D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
           D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
           D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
           D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
           D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
           D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
           D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
           D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
         ME*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
           D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
           D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
           D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) - 
       2.*(ME + MM)*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
     4.*pow(m_alpha,2)*ME2*MM2*(8.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 
       4.*(ME2 + MM2 - m_s)*D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       8.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*(2.*ME2 + ME*MM - 2.*(ME2 - MM2) - 4.*MM2)*D0i(dd13, m_t, ME2, m_s, MM2, 
         ME2, MM2, 0, 0, ME2, MM2) + 4.*ME*MM*(D0i(dd1, m_t, ME2, m_s, MM2, ME2, 
          MM2, 0, 0, ME2, MM2) + D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
          ME2, MM2) + D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
           D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
         D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       4.*(ME*MM + MM2 - m_s - m_t)*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + (4.*ME*MM - 6.*MM2 + 3.*m_t)*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 
         0, 0, ME2, MM2) - (4.*ME2 - 2.*ME*MM + 3.*m_t - 4.*m_u)*
        D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*((ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) - (MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) - MM2*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*ME2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME2*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         ME2*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         (ME2 + MM2 - m_s - m_t)*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
           MM2))) + 2.*pow(m_alpha,2)*ME2*m_t*(8.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 
       4.*(ME2 + MM2 - m_s)*D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       8.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*(2.*ME2 + ME*MM - 2.*(ME2 - MM2) - 4.*MM2)*D0i(dd13, m_t, ME2, m_s, MM2, 
         ME2, MM2, 0, 0, ME2, MM2) + 4.*ME*MM*(D0i(dd1, m_t, ME2, m_s, MM2, ME2, 
          MM2, 0, 0, ME2, MM2) + D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
          ME2, MM2) + D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
           D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
         D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       4.*(ME*MM + MM2 - m_s - m_t)*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + (4.*ME*MM - 6.*MM2 + 3.*m_t)*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 
         0, 0, ME2, MM2) - (4.*ME2 - 2.*ME*MM + 3.*m_t - 4.*m_u)*
        D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*((ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) - (MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) - MM2*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*ME2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME2*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         ME2*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         (ME2 + MM2 - m_s - m_t)*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
           MM2))) - 16.*pow(m_alpha,2)*ME2*MM*MM2*
      (-2.*ME*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       MM*(2.*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         6.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         4.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) + 
     4.*pow(m_alpha,2)*MM*m_t*m_u*(-2.*ME*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
          MM2) + D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       MM*(2.*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         6.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         4.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) - 
     4.*pow(m_alpha,2)*ME*ME2*MM2*(-2.*ME*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, 
          ME2, MM2) + D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       MM*(2.*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         4.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) - 
     4.*pow(m_alpha,2)*ME*pow(MM2,2)*(-2.*ME*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
          MM2) + D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       MM*(2.*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         4.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) + 
     4.*pow(m_alpha,2)*ME*MM2*m_s*(-2.*ME*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
          MM2) + D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       MM*(2.*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         4.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) + 
     4.*pow(m_alpha,2)*ME*(pow(MM2,2) + (MM2 - m_t)*(ME2 - m_u))*
      (-2.*ME*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       MM*(2.*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         4.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) + 
     4.*pow(m_alpha,2)*ME*ME2*MM2*(2.*(ME - MM)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) + 2.*(ME - MM)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) + 2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - 2.*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       4.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       6.*ME*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       5*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     20*pow(m_alpha,2)*ME2*MM*MM2*(2.*(ME - MM)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
         0, ME2, MM2) + 2.*(ME - MM)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) + 2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - 2.*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       4.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       6.*ME*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       5*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     4.*pow(m_alpha,2)*ME*pow(MM2,2)*(2.*(ME - MM)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) + 2.*(ME - MM)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) + 2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - 2.*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       4.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       6.*ME*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       5*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
     8.*pow(m_alpha,2)*ME2*MM*m_t*(2.*(ME - MM)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) + 2.*(ME - MM)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) + 2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - 2.*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       4.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       6.*ME*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       5*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
     2.*pow(m_alpha,2)*MM*MM2*m_t*(2.*(ME - MM)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) + 2.*(ME - MM)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) + 2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - 2.*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       4.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       6.*ME*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       5*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     4.*pow(m_alpha,2)*ME*MM2*(ME2 + MM2 - m_u)*(2.*(ME - MM)*D0i(dd1, m_t, ME2, m_s, MM2, 
         ME2, MM2, 0, 0, ME2, MM2) + 2.*(ME - MM)*D0i(dd1, m_t, ME2, m_u, MM2, 
         ME2, MM2, 0, 0, ME2, MM2) + 2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 
         0, 0, ME2, MM2) - 2.*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) + 2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - 2.*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       4.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       6.*ME*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       5*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
     4.*pow(m_alpha,2)*ME*MM2*m_u*(2.*(ME - MM)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) + 2.*(ME - MM)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) + 2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - 2.*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       4.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       6.*ME*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       5*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     2.*pow(m_alpha,2)*MM*(-2.*(pow(ME2,2) + ME2*(2.*MM2 - m_s)) + MM2*m_s + 
       (2.*ME2 + MM2 - 2.*m_s)*m_u)*(2.*(ME - MM)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 
         0, 0, ME2, MM2) + 2.*(ME - MM)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 
         0, ME2, MM2) + 2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - 2.*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       4.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*ME*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       6.*ME*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       3.*ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       5*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
     12.*pow(m_alpha,2)*ME*ME2*MM2*(2.*ME*(D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, 
          ME2, MM2) + D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       MM*(2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) - 
     12.*pow(m_alpha,2)*ME*pow(MM2,2)*(2.*ME*(D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
          MM2) + D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       MM*(2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) + 
     12.*pow(m_alpha,2)*ME*MM2*m_s*(2.*ME*(D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
          MM2) + D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       MM*(2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) - 
     4.*pow(m_alpha,2)*ME*(-3.*(ME2*MM2 + pow(MM2,2)) + (ME2 + 2.*MM2)*m_t + (3.*MM2 - m_t)*m_u)*
      (2.*ME*(D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       MM*(2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) - 
     40*pow(m_alpha,2)*ME2*MM*MM2*(ME*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       MM*(2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) + 
     16.*pow(m_alpha,2)*ME2*MM*m_t*(ME*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       MM*(2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) + 
     4.*pow(m_alpha,2)*MM*MM2*m_t*(ME*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) + D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       MM*(2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) - 
     4.*pow(m_alpha,2)*MM*(-2.*(pow(ME2,2) + ME2*(2.*MM2 - m_s)) + MM2*m_s + 
       (2.*ME2 + MM2 - 2.*m_s)*m_u)*
      (ME*(-D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       MM*(2.*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2))) - 
     16.*pow(m_alpha,2)*ME2*MM*MM2*(2.*ME*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + 2.*ME*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       6.*MM*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       6.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       6.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       6.*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       4.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       6.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
     4.*pow(m_alpha,2)*MM*m_t*m_u*(2.*ME*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + 2.*ME*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       6.*MM*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       6.*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       6.*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       ME*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       6.*MM*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       4.*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       4.*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       6.*MM*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       2.*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) - 
     4.*pow(m_alpha,2)*ME*ME2*MM*(2.*(ME2 + 5*MM2 - m_s - 2.*m_t)*
        (D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       2.*(2.*ME2 - ME*MM + 5*MM2 - m_s - 2.*m_t)*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, 
          MM2, 0, 0, ME2, MM2) + D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
          ME2, MM2)) + 2.*(ME2 - ME*MM + 4.*MM2 - m_s - 2.*m_t)*
        D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*(ME2 + ME*MM + 6.*MM2 - m_s - 2.*m_t)*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 
         0, 0, ME2, MM2) + 4.*ME*MM*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) + (2.*ME2 - 4.*ME*MM + 4.*MM2 - 2.*m_s - 5*m_t)*
        D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (3.*(ME2 + MM2) + 2.*m_t - m_u)*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) - 2.*MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + 2.*((ME2 - 2.*ME*MM)*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
           ME2, MM2) + 2.*ME2*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) + 3.*ME2*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         4.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))) - 
     4.*pow(m_alpha,2)*ME*MM*MM2*(2.*(ME2 + 5*MM2 - m_s - 2.*m_t)*
        (D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       2.*(2.*ME2 - ME*MM + 5*MM2 - m_s - 2.*m_t)*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, 
          MM2, 0, 0, ME2, MM2) + D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
          ME2, MM2)) + 2.*(ME2 - ME*MM + 4.*MM2 - m_s - 2.*m_t)*
        D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*(ME2 + ME*MM + 6.*MM2 - m_s - 2.*m_t)*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 
         0, 0, ME2, MM2) + 4.*ME*MM*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) + (2.*ME2 - 4.*ME*MM + 4.*MM2 - 2.*m_s - 5*m_t)*
        D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (3.*(ME2 + MM2) + 2.*m_t - m_u)*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) - 2.*MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + 2.*((ME2 - 2.*ME*MM)*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
           ME2, MM2) + 2.*ME2*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) + 3.*ME2*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         4.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))) + 
     4.*pow(m_alpha,2)*ME*MM*m_u*(2.*(ME2 + 5*MM2 - m_s - 2.*m_t)*
        (D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       2.*(2.*ME2 - ME*MM + 5*MM2 - m_s - 2.*m_t)*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, 
          MM2, 0, 0, ME2, MM2) + D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
          ME2, MM2)) + 2.*(ME2 - ME*MM + 4.*MM2 - m_s - 2.*m_t)*
        D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*(ME2 + ME*MM + 6.*MM2 - m_s - 2.*m_t)*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 
         0, 0, ME2, MM2) + 4.*ME*MM*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) + (2.*ME2 - 4.*ME*MM + 4.*MM2 - 2.*m_s - 5*m_t)*
        D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (3.*(ME2 + MM2) + 2.*m_t - m_u)*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) - 2.*MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + 2.*((ME2 - 2.*ME*MM)*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
           ME2, MM2) + 2.*ME2*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) + 3.*ME2*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         4.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))) + 
     12.*pow(m_alpha,2)*ME*ME2*MM*(2.*(ME2 + MM2 - m_s)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 
         0, 0, ME2, MM2) + 2.*(ME2 + MM2 - m_s)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, 
         MM2, 0, 0, ME2, MM2) + 2.*(2.*ME2 - ME*MM + MM2 - m_s)*
        (D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       2.*(ME2 + ME*MM + 2.*MM2 - m_s)*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) + (2.*ME2 - 4.*ME*MM + 4.*MM2 + 3.*(ME2 + MM2) - m_s - 2.*m_u)*
        D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*((ME2 + MM2 - m_s)*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) + (ME2 + MM2 - m_s)*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
           ME2, MM2) + ME2*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) - ME*MM*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME2*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*ME2*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         3.*ME2*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         ME2*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*MM2*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         m_u*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))) + 
     12.*pow(m_alpha,2)*ME*MM*MM2*(2.*(ME2 + MM2 - m_s)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 
         0, 0, ME2, MM2) + 2.*(ME2 + MM2 - m_s)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, 
         MM2, 0, 0, ME2, MM2) + 2.*(2.*ME2 - ME*MM + MM2 - m_s)*
        (D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       2.*(ME2 + ME*MM + 2.*MM2 - m_s)*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) + (2.*ME2 - 4.*ME*MM + 4.*MM2 + 3.*(ME2 + MM2) - m_s - 2.*m_u)*
        D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*((ME2 + MM2 - m_s)*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) + (ME2 + MM2 - m_s)*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
           ME2, MM2) + ME2*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) - ME*MM*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME2*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*ME2*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         3.*ME2*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         ME2*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*MM2*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         m_u*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))) - 
     12.*pow(m_alpha,2)*ME*MM*m_u*(2.*(ME2 + MM2 - m_s)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 
         0, 0, ME2, MM2) + 2.*(ME2 + MM2 - m_s)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, 
         MM2, 0, 0, ME2, MM2) + 2.*(2.*ME2 - ME*MM + MM2 - m_s)*
        (D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       2.*(ME2 + ME*MM + 2.*MM2 - m_s)*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) + (2.*ME2 - 4.*ME*MM + 4.*MM2 + 3.*(ME2 + MM2) - m_s - 2.*m_u)*
        D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*((ME2 + MM2 - m_s)*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) + (ME2 + MM2 - m_s)*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
           ME2, MM2) + ME2*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) - ME*MM*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME2*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*ME2*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         3.*ME2*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         ME2*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*MM2*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         m_u*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME*MM*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))) + 
     8.*pow(m_alpha,2)*ME*ME2*MM*(2.*(ME*MM - MM2)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
         0, ME2, MM2) + 2.*(ME*MM - MM2)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 
         0, ME2, MM2) + 2.*ME2*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
          MM2) + D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       (2.*ME2 + 5*ME*MM - 5*MM2)*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) + 2.*((ME*MM - MM2)*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
           0, ME2, MM2) + (ME*MM - MM2)*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 
           0, 0, ME2, MM2) + ME*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, 
           ME2, MM2) - MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) + ME*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*ME*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*MM2*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         ME*MM*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         ME*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*MM2*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         ME*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         MM2*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))) + 
     8.*pow(m_alpha,2)*ME*MM*MM2*(2.*(ME*MM - MM2)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
         0, ME2, MM2) + 2.*(ME*MM - MM2)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 
         0, ME2, MM2) + 2.*ME2*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
          MM2) + D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       (2.*ME2 + 5*ME*MM - 5*MM2)*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) + 2.*((ME*MM - MM2)*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
           0, ME2, MM2) + (ME*MM - MM2)*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 
           0, 0, ME2, MM2) + ME*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, 
           ME2, MM2) - MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) + ME*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*ME*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*MM2*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         ME*MM*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         ME*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*MM2*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         ME*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         MM2*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))) - 
     8.*pow(m_alpha,2)*ME*MM*m_s*(2.*(ME*MM - MM2)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
         0, ME2, MM2) + 2.*(ME*MM - MM2)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 
         0, ME2, MM2) + 2.*ME2*(D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
          MM2) + D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       (2.*ME2 + 5*ME*MM - 5*MM2)*D0i(dd23, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, 
         ME2, MM2) + 2.*((ME*MM - MM2)*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 
           0, ME2, MM2) + (ME*MM - MM2)*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 
           0, 0, ME2, MM2) + ME*MM*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, 
           ME2, MM2) - MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) + ME*MM*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*ME*MM*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*MM2*D0i(dd13, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         ME*MM*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         ME*MM*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*MM2*D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         ME*MM*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         MM2*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME*MM*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))) + 
     4.*pow(m_alpha,2)*pow(ME2,2)*(4.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 
       2.*(ME2 + MM2 - m_s)*D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*(ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - 12.*D0i(dd00, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       16.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*(ME*MM + 2.*(MM2 - m_t))*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + 2.*(ME*MM + 2.*MM2 - 2.*(MM2 + m_t))*D0i(dd13, m_t, ME2, m_u, MM2, 
         ME2, MM2, 0, 0, ME2, MM2) + 2.*(ME2 + ME*MM - MM2 + m_t - m_u)*
        D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (4.*ME2 - 3.*ME*MM + MM2 - 4.*(2.*ME2 + m_t - m_u))*D0i(dd23, m_t, ME2, m_u, MM2, 
         ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME2*(2.*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       2.*((ME*MM + MM2)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         (ME*MM + MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) + ME*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         ME*MM*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_t*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_t*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*MM2*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))) + 
     24.*pow(m_alpha,2)*ME2*MM2*(4.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 
       2.*(ME2 + MM2 - m_s)*D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*(ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - 12.*D0i(dd00, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       16.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*(ME*MM + 2.*(MM2 - m_t))*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + 2.*(ME*MM + 2.*MM2 - 2.*(MM2 + m_t))*D0i(dd13, m_t, ME2, m_u, MM2, 
         ME2, MM2, 0, 0, ME2, MM2) + 2.*(ME2 + ME*MM - MM2 + m_t - m_u)*
        D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (4.*ME2 - 3.*ME*MM + MM2 - 4.*(2.*ME2 + m_t - m_u))*D0i(dd23, m_t, ME2, m_u, MM2, 
         ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME2*(2.*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       2.*((ME*MM + MM2)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         (ME*MM + MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) + ME*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         ME*MM*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_t*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_t*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*MM2*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))) + 
     4.*pow(m_alpha,2)*pow(MM2,2)*(4.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 
       2.*(ME2 + MM2 - m_s)*D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*(ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - 12.*D0i(dd00, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       16.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*(ME*MM + 2.*(MM2 - m_t))*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + 2.*(ME*MM + 2.*MM2 - 2.*(MM2 + m_t))*D0i(dd13, m_t, ME2, m_u, MM2, 
         ME2, MM2, 0, 0, ME2, MM2) + 2.*(ME2 + ME*MM - MM2 + m_t - m_u)*
        D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (4.*ME2 - 3.*ME*MM + MM2 - 4.*(2.*ME2 + m_t - m_u))*D0i(dd23, m_t, ME2, m_u, MM2, 
         ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME2*(2.*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       2.*((ME*MM + MM2)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         (ME*MM + MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) + ME*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         ME*MM*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_t*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_t*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*MM2*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))) - 
     8.*pow(m_alpha,2)*ME2*m_u*(4.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 
       2.*(ME2 + MM2 - m_s)*D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*(ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - 12.*D0i(dd00, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       16.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*(ME*MM + 2.*(MM2 - m_t))*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + 2.*(ME*MM + 2.*MM2 - 2.*(MM2 + m_t))*D0i(dd13, m_t, ME2, m_u, MM2, 
         ME2, MM2, 0, 0, ME2, MM2) + 2.*(ME2 + ME*MM - MM2 + m_t - m_u)*
        D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (4.*ME2 - 3.*ME*MM + MM2 - 4.*(2.*ME2 + m_t - m_u))*D0i(dd23, m_t, ME2, m_u, MM2, 
         ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME2*(2.*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       2.*((ME*MM + MM2)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         (ME*MM + MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) + ME*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         ME*MM*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_t*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_t*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*MM2*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))) - 
     8.*pow(m_alpha,2)*MM2*m_u*(4.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 
       2.*(ME2 + MM2 - m_s)*D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*(ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - 12.*D0i(dd00, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       16.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*(ME*MM + 2.*(MM2 - m_t))*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + 2.*(ME*MM + 2.*MM2 - 2.*(MM2 + m_t))*D0i(dd13, m_t, ME2, m_u, MM2, 
         ME2, MM2, 0, 0, ME2, MM2) + 2.*(ME2 + ME*MM - MM2 + m_t - m_u)*
        D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (4.*ME2 - 3.*ME*MM + MM2 - 4.*(2.*ME2 + m_t - m_u))*D0i(dd23, m_t, ME2, m_u, MM2, 
         ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME2*(2.*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       2.*((ME*MM + MM2)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         (ME*MM + MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) + ME*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         ME*MM*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_t*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_t*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*MM2*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2))) + 
     4.*pow(m_alpha,2)*pow(m_u,2)*(4.*C0i(cc0, ME2, m_s, MM2, 0, ME2, MM2) + 
       2.*(ME2 + MM2 - m_s)*D0i(dd0, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*(ME2 + MM2 - m_s - m_t)*D0i(dd0, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) - 12.*D0i(dd00, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
       16.*D0i(dd00, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       2.*(ME*MM + 2.*(MM2 - m_t))*D0i(dd13, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, 
         MM2) + 2.*(ME*MM + 2.*MM2 - 2.*(MM2 + m_t))*D0i(dd13, m_t, ME2, m_u, MM2, 
         ME2, MM2, 0, 0, ME2, MM2) + 2.*(ME2 + ME*MM - MM2 + m_t - m_u)*
        D0i(dd23, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
       (4.*ME2 - 3.*ME*MM + MM2 - 4.*(2.*ME2 + m_t - m_u))*D0i(dd23, m_t, ME2, m_u, MM2, 
         ME2, MM2, 0, 0, ME2, MM2) + 
       2.*ME2*(2.*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)) + 
       2.*((ME*MM + MM2)*D0i(dd1, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         (ME*MM + MM2 - 2.*m_t)*D0i(dd1, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, 
           MM2) + ME*MM*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd11, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         ME*MM*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd11, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd12, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*m_t*D0i(dd12, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd2, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         ME*MM*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         MM2*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_t*D0i(dd2, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         2.*ME*MM*D0i(dd22, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*ME*MM*D0i(dd22, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) + 
         2.*MM2*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd3, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_s*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         m_t*D0i(dd3, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         MM2*D0i(dd33, m_t, ME2, m_s, MM2, ME2, MM2, 0, 0, ME2, MM2) - 
         3.*MM2*D0i(dd33, m_t, ME2, m_u, MM2, ME2, MM2, 0, 0, ME2, MM2)))))/m_t)/32.;

}


Complex emu_emu_Virtual::BornTriangle(){
return 1/16.*(-64.*pow(m_alpha,3)*ME*ME2*MM*(ME2 + ME*MM)*M_PI*
    (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc11, ME2, m_t, ME2, 0, ME2, 
      ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
     C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, 
      ME2)))/pow(m_t,2) + (256*pow(m_alpha,3)*pow(ME2,2)*MM2*M_PI*
    (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc11, ME2, m_t, ME2, 0, ME2, 
      ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
     C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, 
      ME2)))/pow(m_t,2) - (128.*pow(m_alpha,3)*ME*ME2*MM*MM2*M_PI*
    (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc11, ME2, m_t, ME2, 0, ME2, 
      ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
     C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, 
      ME2)))/pow(m_t,2) - (64.*pow(m_alpha,3)*ME*MM*(ME2 + ME*MM)*MM2*M_PI*
    (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc11, ME2, m_t, ME2, 0, ME2, 
      ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
     C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, 
      ME2)))/pow(m_t,2) + (256*pow(m_alpha,3)*ME2*pow(MM2,2)*M_PI*
    (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc11, ME2, m_t, ME2, 0, ME2, 
      ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
     C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, 
      ME2)))/pow(m_t,2) - (64.*pow(m_alpha,3)*ME*MM*(ME2 + ME*MM)*M_PI*m_s*
    (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc11, ME2, m_t, ME2, 0, ME2, 
      ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
     C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, 
      ME2)))/pow(m_t,2) - (64.*pow(m_alpha,3)*ME2*MM2*M_PI*m_s*
    (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc11, ME2, m_t, ME2, 0, ME2, 
      ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
     C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, 
      ME2)))/pow(m_t,2) + (128.*pow(m_alpha,3)*ME*ME2*MM*M_PI*
    (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc11, ME2, m_t, ME2, 0, ME2, 
      ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
     C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, 
      ME2)))/m_t + (96*pow(m_alpha,3)*ME2*M_PI*(-2.*(ME2*MM2 + pow(MM2,2)) + 2.*MM2*m_s + 
     (ME2 + MM2 - m_s)*m_t)*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
     C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 
     2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, 
      ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)))/pow(m_t,2) - 
  (32.*pow(m_alpha,3)*ME2*M_PI*(pow(MM2,2) + (MM2 - m_t)*(ME2 - m_u))*
    (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc11, ME2, m_t, ME2, 0, ME2, 
      ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
     C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, 
      ME2)))/pow(m_t,2) + (4.*pow(m_alpha,2)*ME2*((8.*m_alpha*MM2*M_PI*(ME2 + MM2 - m_s))/m_t + 
     (8.*m_alpha*M_PI*(-2.*(ME2*MM2 + pow(MM2,2)) + 2.*MM2*m_s + (ME2 + MM2 - m_s)*m_t))/m_t - 
     (8.*m_alpha*M_PI*(pow(MM2,2) + (MM2 - m_t)*(ME2 - m_u)))/m_t + 
     (16.*m_alpha*MM2*M_PI*(ME2 + MM2 - m_u))/m_t)*
    (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc11, ME2, m_t, ME2, 0, ME2, 
      ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
     C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, 
      ME2)))/m_t + (128.*pow(m_alpha,3)*ME*MM*(ME2 + ME*MM)*M_PI*m_u*
    (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc11, ME2, m_t, ME2, 0, ME2, 
      ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
     C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, 
      ME2)))/pow(m_t,2) - (192.*pow(m_alpha,3)*ME2*MM2*M_PI*m_u*
    (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc11, ME2, m_t, ME2, 0, ME2, 
      ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
     C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, 
      ME2)))/pow(m_t,2) - (32.*pow(m_alpha,3)*ME*MM*M_PI*m_u*
    (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc11, ME2, m_t, ME2, 0, ME2, 
      ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
     C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, 
      ME2)))/m_t + (64.*pow(m_alpha,3)*ME*MM*M_PI*(pow(ME2,2) + 2.*ME2*MM2 + pow(MM2,2) + m_s*m_u - 
     (ME2 + MM2)*(m_s + m_u))*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
     C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 
     2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, 
      ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)))/pow(m_t,2) - 
  (128.*pow(m_alpha,3)*ME2*pow(MM2,2)*M_PI*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
     C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 
     2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, 
      MM2) + C0i(cc22, MM2, m_t, MM2, 0, MM2, MM2)))/pow(m_t,2) + 
  (128.*pow(m_alpha,3)*ME2*MM2*M_PI*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
     C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 
     2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, 
      MM2) + C0i(cc22, MM2, m_t, MM2, 0, MM2, MM2)))/m_t - 
  (32.*pow(m_alpha,3)*MM2*M_PI*m_u*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
     C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 
     2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, 
      MM2) + C0i(cc22, MM2, m_t, MM2, 0, MM2, MM2)))/m_t + 
  (64.*pow(m_alpha,3)*MM2*M_PI*(pow(ME2,2) + 2.*ME2*MM2 + pow(MM2,2) + m_s*m_u - 
     (ME2 + MM2)*(m_s + m_u))*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
     C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 
     2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, 
      MM2) + C0i(cc22, MM2, m_t, MM2, 0, MM2, MM2)))/pow(m_t,2) + 
  (32.*m_alpha*ME*ME2*MM2*M_PI*(-2.*pow(m_alpha,2)*ME*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
       C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, 
         ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 
        0, ME2, ME2)) + 2.*pow(m_alpha,2)*MM*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
       C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, 
         MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 
        0, MM2, MM2))))/pow(m_t,2) - 
  (16.*m_alpha*ME2*MM*MM2*M_PI*(-2.*pow(m_alpha,2)*ME*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
       C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, 
         ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 
        0, ME2, ME2)) + 2.*pow(m_alpha,2)*MM*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
       C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, 
         MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 
        0, MM2, MM2))))/pow(m_t,2) + 
  (32.*m_alpha*ME*pow(MM2,2)*M_PI*(-2.*pow(m_alpha,2)*ME*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
       C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, 
         ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 
        0, ME2, ME2)) + 2.*pow(m_alpha,2)*MM*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
       C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, 
         MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 
        0, MM2, MM2))))/pow(m_t,2) - 
  (64.*m_alpha*ME*MM2*M_PI*m_s*(-2.*pow(m_alpha,2)*ME*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
       C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, 
         ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 
        0, ME2, ME2)) + 2.*pow(m_alpha,2)*MM*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
       C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, 
         MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 
        0, MM2, MM2))))/pow(m_t,2) + 
  (16.*m_alpha*MM*M_PI*(pow(ME2,2) + 2.*ME2*MM2 + pow(MM2,2) - 2.*(ME2 + MM2)*m_s + pow(m_s,2))*
    (-2.*pow(m_alpha,2)*ME*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
       C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, 
         ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 
        0, ME2, ME2)) + 2.*pow(m_alpha,2)*MM*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
       C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, 
         MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 
        0, MM2, MM2))))/pow(m_t,2) + 
  (32.*m_alpha*ME2*MM*M_PI*(-2.*pow(m_alpha,2)*ME*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
       C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, 
         ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 
        0, ME2, ME2)) + 2.*pow(m_alpha,2)*MM*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
       C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, 
         MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 
        0, MM2, MM2))))/m_t + 
  (8.*m_alpha*MM*MM2*M_PI*(-2.*pow(m_alpha,2)*ME*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
       C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, 
         ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 
        0, ME2, ME2)) + 2.*pow(m_alpha,2)*MM*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
       C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, 
         MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 
        0, MM2, MM2))))/m_t + 
  (32.*m_alpha*ME*MM2*M_PI*m_u*(-2.*pow(m_alpha,2)*ME*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
       C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, 
         ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 
        0, ME2, ME2)) + 2.*pow(m_alpha,2)*MM*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
       C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, 
         MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 
        0, MM2, MM2))))/pow(m_t,2) - 
  (8.*m_alpha*MM*M_PI*(-2.*(pow(ME2,2) + ME2*(2.*MM2 - m_s)) + MM2*m_s + 
     (2.*ME2 + MM2 - 2.*m_s)*m_u)*
    (-2.*pow(m_alpha,2)*ME*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
       C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, 
         ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 
        0, ME2, ME2)) + 2.*pow(m_alpha,2)*MM*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
       C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, 
         MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 
        0, MM2, MM2))))/pow(m_t,2) - 
  (16.*m_alpha*MM*M_PI*(pow(ME2,2) + pow(MM2,2) + 2.*(ME2 + MM2)*m_s - pow(m_s,2) - 2.*(ME2*MM2 + m_s*m_u))*
    (-2.*pow(m_alpha,2)*ME*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
       C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, 
         ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 
        0, ME2, ME2)) + 2.*pow(m_alpha,2)*MM*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
       C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, 
         MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 
        0, MM2, MM2))))/pow(m_t,2) + (96*pow(m_alpha,3)*ME*ME2*MM*M_PI*
    (ME2*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc11, ME2, m_t, ME2, 0, ME2, 
        ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
       C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, 
        ME2)) + MM2*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
       C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, 
         MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 
        0, MM2, MM2))))/pow(m_t,2) + (96*pow(m_alpha,3)*ME*MM*MM2*M_PI*
    (ME2*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc11, ME2, m_t, ME2, 0, ME2, 
        ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
       C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, 
        ME2)) + MM2*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
       C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, 
         MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 
        0, MM2, MM2))))/pow(m_t,2) - (96*pow(m_alpha,3)*ME*MM*M_PI*m_u*
    (ME2*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc11, ME2, m_t, ME2, 0, ME2, 
        ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
       C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, 
        ME2)) + MM2*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
       C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, 
         MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 
        0, MM2, MM2))))/pow(m_t,2) - (160*pow(m_alpha,3)*ME*ME2*MM*M_PI*
    ((ME2 - 4.*ME*MM)*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
       C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, 
         ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 
        0, ME2, ME2)) + MM2*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
       C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, 
         MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 
        0, MM2, MM2))))/pow(m_t,2) - (160*pow(m_alpha,3)*ME*MM*MM2*M_PI*
    ((ME2 - 4.*ME*MM)*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
       C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, 
         ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 
        0, ME2, ME2)) + MM2*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
       C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, 
         MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 
        0, MM2, MM2))))/pow(m_t,2) + (128.*pow(m_alpha,3)*ME*MM*M_PI*m_s*
    ((ME2 - 4.*ME*MM)*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
       C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, 
         ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 
        0, ME2, ME2)) + MM2*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
       C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, 
         MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 
        0, MM2, MM2))))/pow(m_t,2) + (32.*pow(m_alpha,3)*ME*MM*M_PI*m_u*
    ((ME2 - 4.*ME*MM)*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
       C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, 
         ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 
        0, ME2, ME2)) + MM2*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
       C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, 
         MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 
        0, MM2, MM2))))/pow(m_t,2) - 
  (32.*m_alpha*pow(ME2,2)*M_PI*(pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - 
       B0i(bb0, m_t, MM2, MM2) + (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, 
         ME2) + (-2.*MM2 + m_t)*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
       2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
       2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 + ME*MM)*
        C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 2.*ME2*C0i(cc11, ME2, m_t, ME2, 0, 
         ME2, ME2) - ME*MM*C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 
       4.*ME2*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) - 
       2.*ME*MM*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) - 
       2.*ME2*C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) - 
       ME*MM*C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
       m_t*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc1, MM2, m_t, MM2, 0, MM2, 
          MM2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, MM2, m_t, MM2, 
          0, MM2, MM2)) + (2.*ME2 - ME*MM)*C0i(cc22, ME2, m_t, ME2, 0, ME2, 
         ME2) + MM2*(-3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, 
           MM2, MM2) - 3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc22, MM2, m_t, MM2, 0, MM2, MM2))) - 
     m_alpha*M_PI*((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - 
            Re(B0i(bb1, ME2, 0, ME2)))/2))/M_PI + 
       (2.*m_alpha*(1./4 + (-Re(B0i(bb0, MM2, 0, MM2)) - 
            Re(B0i(bb1, MM2, 0, MM2)))/2))/M_PI + 
       (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
       (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
       (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
       (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)))/pow(m_t,2) - 
  (128.*m_alpha*ME2*MM2*M_PI*(pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - 
       B0i(bb0, m_t, MM2, MM2) + (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, 
         ME2) + (-2.*MM2 + m_t)*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
       2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
       2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 + ME*MM)*
        C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 2.*ME2*C0i(cc11, ME2, m_t, ME2, 0, 
         ME2, ME2) - ME*MM*C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 
       4.*ME2*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) - 
       2.*ME*MM*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) - 
       2.*ME2*C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) - 
       ME*MM*C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
       m_t*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc1, MM2, m_t, MM2, 0, MM2, 
          MM2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, MM2, m_t, MM2, 
          0, MM2, MM2)) + (2.*ME2 - ME*MM)*C0i(cc22, ME2, m_t, ME2, 0, ME2, 
         ME2) + MM2*(-3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, 
           MM2, MM2) - 3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc22, MM2, m_t, MM2, 0, MM2, MM2))) - 
     m_alpha*M_PI*((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - 
            Re(B0i(bb1, ME2, 0, ME2)))/2))/M_PI + 
       (2.*m_alpha*(1./4 + (-Re(B0i(bb0, MM2, 0, MM2)) - 
            Re(B0i(bb1, MM2, 0, MM2)))/2))/M_PI + 
       (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
       (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
       (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
       (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)))/pow(m_t,2) - 
  (32.*m_alpha*pow(MM2,2)*M_PI*(pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - 
       B0i(bb0, m_t, MM2, MM2) + (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, 
         ME2) + (-2.*MM2 + m_t)*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
       2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
       2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 + ME*MM)*
        C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 2.*ME2*C0i(cc11, ME2, m_t, ME2, 0, 
         ME2, ME2) - ME*MM*C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 
       4.*ME2*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) - 
       2.*ME*MM*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) - 
       2.*ME2*C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) - 
       ME*MM*C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
       m_t*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc1, MM2, m_t, MM2, 0, MM2, 
          MM2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, MM2, m_t, MM2, 
          0, MM2, MM2)) + (2.*ME2 - ME*MM)*C0i(cc22, ME2, m_t, ME2, 0, ME2, 
         ME2) + MM2*(-3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, 
           MM2, MM2) - 3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc22, MM2, m_t, MM2, 0, MM2, MM2))) - 
     m_alpha*M_PI*((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - 
            Re(B0i(bb1, ME2, 0, ME2)))/2))/M_PI + 
       (2.*m_alpha*(1./4 + (-Re(B0i(bb0, MM2, 0, MM2)) - 
            Re(B0i(bb1, MM2, 0, MM2)))/2))/M_PI + 
       (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
       (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
       (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
       (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)))/pow(m_t,2) + 
  (64.*m_alpha*ME2*M_PI*m_s*(pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - 
       B0i(bb0, m_t, MM2, MM2) + (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, 
         ME2) + (-2.*MM2 + m_t)*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
       2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
       2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 + ME*MM)*
        C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 2.*ME2*C0i(cc11, ME2, m_t, ME2, 0, 
         ME2, ME2) - ME*MM*C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 
       4.*ME2*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) - 
       2.*ME*MM*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) - 
       2.*ME2*C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) - 
       ME*MM*C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
       m_t*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc1, MM2, m_t, MM2, 0, MM2, 
          MM2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, MM2, m_t, MM2, 
          0, MM2, MM2)) + (2.*ME2 - ME*MM)*C0i(cc22, ME2, m_t, ME2, 0, ME2, 
         ME2) + MM2*(-3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, 
           MM2, MM2) - 3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc22, MM2, m_t, MM2, 0, MM2, MM2))) - 
     m_alpha*M_PI*((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - 
            Re(B0i(bb1, ME2, 0, ME2)))/2))/M_PI + 
       (2.*m_alpha*(1./4 + (-Re(B0i(bb0, MM2, 0, MM2)) - 
            Re(B0i(bb1, MM2, 0, MM2)))/2))/M_PI + 
       (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
       (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
       (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
       (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)))/pow(m_t,2) + 
  (64.*m_alpha*MM2*M_PI*m_s*(pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - 
       B0i(bb0, m_t, MM2, MM2) + (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, 
         ME2) + (-2.*MM2 + m_t)*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
       2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
       2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 + ME*MM)*
        C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 2.*ME2*C0i(cc11, ME2, m_t, ME2, 0, 
         ME2, ME2) - ME*MM*C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 
       4.*ME2*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) - 
       2.*ME*MM*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) - 
       2.*ME2*C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) - 
       ME*MM*C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
       m_t*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc1, MM2, m_t, MM2, 0, MM2, 
          MM2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, MM2, m_t, MM2, 
          0, MM2, MM2)) + (2.*ME2 - ME*MM)*C0i(cc22, ME2, m_t, ME2, 0, ME2, 
         ME2) + MM2*(-3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, 
           MM2, MM2) - 3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc22, MM2, m_t, MM2, 0, MM2, MM2))) - 
     m_alpha*M_PI*((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - 
            Re(B0i(bb1, ME2, 0, ME2)))/2))/M_PI + 
       (2.*m_alpha*(1./4 + (-Re(B0i(bb0, MM2, 0, MM2)) - 
            Re(B0i(bb1, MM2, 0, MM2)))/2))/M_PI + 
       (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
       (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
       (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
       (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)))/pow(m_t,2) - 
  (32.*m_alpha*M_PI*pow(m_s,2)*(pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - 
       B0i(bb0, m_t, MM2, MM2) + (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, 
         ME2) + (-2.*MM2 + m_t)*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
       2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
       2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 + ME*MM)*
        C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 2.*ME2*C0i(cc11, ME2, m_t, ME2, 0, 
         ME2, ME2) - ME*MM*C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 
       4.*ME2*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) - 
       2.*ME*MM*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) - 
       2.*ME2*C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) - 
       ME*MM*C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
       m_t*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc1, MM2, m_t, MM2, 0, MM2, 
          MM2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, MM2, m_t, MM2, 
          0, MM2, MM2)) + (2.*ME2 - ME*MM)*C0i(cc22, ME2, m_t, ME2, 0, ME2, 
         ME2) + MM2*(-3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, 
           MM2, MM2) - 3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc22, MM2, m_t, MM2, 0, MM2, MM2))) - 
     m_alpha*M_PI*((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - 
            Re(B0i(bb1, ME2, 0, ME2)))/2))/M_PI + 
       (2.*m_alpha*(1./4 + (-Re(B0i(bb0, MM2, 0, MM2)) - 
            Re(B0i(bb1, MM2, 0, MM2)))/2))/M_PI + 
       (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
       (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
       (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
       (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)))/pow(m_t,2) - 
  (32.*m_alpha*MM2*M_PI*(pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - 
       B0i(bb0, m_t, MM2, MM2) + (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, 
         ME2) + (-2.*MM2 + m_t)*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
       2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
       2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 + ME*MM)*
        C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 2.*ME2*C0i(cc11, ME2, m_t, ME2, 0, 
         ME2, ME2) - ME*MM*C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 
       4.*ME2*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) - 
       2.*ME*MM*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) - 
       2.*ME2*C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) - 
       ME*MM*C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
       m_t*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc1, MM2, m_t, MM2, 0, MM2, 
          MM2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, MM2, m_t, MM2, 
          0, MM2, MM2)) + (2.*ME2 - ME*MM)*C0i(cc22, ME2, m_t, ME2, 0, ME2, 
         ME2) + MM2*(-3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, 
           MM2, MM2) - 3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc22, MM2, m_t, MM2, 0, MM2, MM2))) - 
     m_alpha*M_PI*((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - 
            Re(B0i(bb1, ME2, 0, ME2)))/2))/M_PI + 
       (2.*m_alpha*(1./4 + (-Re(B0i(bb0, MM2, 0, MM2)) - 
            Re(B0i(bb1, MM2, 0, MM2)))/2))/M_PI + 
       (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
       (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
       (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
       (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)))/m_t + 
  (4.*((-32.*m_alpha*ME2*MM2*M_PI)/m_t - (8.*m_alpha*M_PI*pow(ME2 + MM2 - m_s,2))/m_t - 
     (8.*m_alpha*MM2*M_PI*(-2.*ME2 + m_t))/m_t - (8.*m_alpha*ME2*M_PI*(-2.*MM2 + m_t))/m_t)*
    (pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - B0i(bb0, m_t, MM2, MM2) + 
       (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, ME2) + 
       (-2.*MM2 + m_t)*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
       2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
       2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 + ME*MM)*
        C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 2.*ME2*C0i(cc11, ME2, m_t, ME2, 0, 
         ME2, ME2) - ME*MM*C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 
       4.*ME2*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) - 
       2.*ME*MM*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) - 
       2.*ME2*C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) - 
       ME*MM*C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
       m_t*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc1, MM2, m_t, MM2, 0, MM2, 
          MM2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, MM2, m_t, MM2, 
          0, MM2, MM2)) + (2.*ME2 - ME*MM)*C0i(cc22, ME2, m_t, ME2, 0, ME2, 
         ME2) + MM2*(-3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, 
           MM2, MM2) - 3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc22, MM2, m_t, MM2, 0, MM2, MM2))) - 
     m_alpha*M_PI*((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - 
            Re(B0i(bb1, ME2, 0, ME2)))/2))/M_PI + 
       (2.*m_alpha*(1./4 + (-Re(B0i(bb0, MM2, 0, MM2)) - 
            Re(B0i(bb1, MM2, 0, MM2)))/2))/M_PI + 
       (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
       (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
       (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
       (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)))/m_t - 
  (32.*m_alpha*pow(ME2,2)*M_PI*(pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - 
       B0i(bb0, m_t, MM2, MM2) + (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, 
         ME2) + 2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
       2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 - ME*MM - m_t)*
        (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, 
          ME2)) + m_t*(C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, 
          MM2)) + (2.*ME2 + ME*MM)*(C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 
         2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, 
          ME2, ME2)) + MM2*(-2.*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) - 
         3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 0, 
          MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) - 
         3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, 
          MM2, MM2))) - m_alpha*M_PI*
      ((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - 
            Re(B0i(bb1, ME2, 0, ME2)))/2))/M_PI + 
       (2.*m_alpha*(1./4 + (-Re(B0i(bb0, MM2, 0, MM2)) - 
            Re(B0i(bb1, MM2, 0, MM2)))/2))/M_PI + 
       (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
       (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
       (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
       (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)))/pow(m_t,2) + 
  (64.*m_alpha*ME2*MM2*M_PI*(pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - 
       B0i(bb0, m_t, MM2, MM2) + (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, 
         ME2) + 2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
       2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 - ME*MM - m_t)*
        (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, 
          ME2)) + m_t*(C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, 
          MM2)) + (2.*ME2 + ME*MM)*(C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 
         2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, 
          ME2, ME2)) + MM2*(-2.*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) - 
         3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 0, 
          MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) - 
         3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, 
          MM2, MM2))) - m_alpha*M_PI*
      ((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - 
            Re(B0i(bb1, ME2, 0, ME2)))/2))/M_PI + 
       (2.*m_alpha*(1./4 + (-Re(B0i(bb0, MM2, 0, MM2)) - 
            Re(B0i(bb1, MM2, 0, MM2)))/2))/M_PI + 
       (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
       (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
       (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
       (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)))/pow(m_t,2) - 
  (32.*m_alpha*pow(MM2,2)*M_PI*(pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - 
       B0i(bb0, m_t, MM2, MM2) + (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, 
         ME2) + 2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
       2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 - ME*MM - m_t)*
        (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, 
          ME2)) + m_t*(C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, 
          MM2)) + (2.*ME2 + ME*MM)*(C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 
         2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, 
          ME2, ME2)) + MM2*(-2.*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) - 
         3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 0, 
          MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) - 
         3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, 
          MM2, MM2))) - m_alpha*M_PI*
      ((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - 
            Re(B0i(bb1, ME2, 0, ME2)))/2))/M_PI + 
       (2.*m_alpha*(1./4 + (-Re(B0i(bb0, MM2, 0, MM2)) - 
            Re(B0i(bb1, MM2, 0, MM2)))/2))/M_PI + 
       (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
       (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
       (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
       (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)))/pow(m_t,2) - 
  (64.*m_alpha*ME2*M_PI*(pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - 
       B0i(bb0, m_t, MM2, MM2) + (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, 
         ME2) + 2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
       2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 - ME*MM - m_t)*
        (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, 
          ME2)) + m_t*(C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, 
          MM2)) + (2.*ME2 + ME*MM)*(C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 
         2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, 
          ME2, ME2)) + MM2*(-2.*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) - 
         3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 0, 
          MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) - 
         3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, 
          MM2, MM2))) - m_alpha*M_PI*
      ((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - 
            Re(B0i(bb1, ME2, 0, ME2)))/2))/M_PI + 
       (2.*m_alpha*(1./4 + (-Re(B0i(bb0, MM2, 0, MM2)) - 
            Re(B0i(bb1, MM2, 0, MM2)))/2))/M_PI + 
       (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
       (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
       (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
       (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)))/m_t - 
  (64.*m_alpha*MM2*M_PI*(pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - 
       B0i(bb0, m_t, MM2, MM2) + (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, 
         ME2) + 2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
       2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 - ME*MM - m_t)*
        (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, 
          ME2)) + m_t*(C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, 
          MM2)) + (2.*ME2 + ME*MM)*(C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 
         2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, 
          ME2, ME2)) + MM2*(-2.*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) - 
         3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 0, 
          MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) - 
         3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, 
          MM2, MM2))) - m_alpha*M_PI*
      ((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - 
            Re(B0i(bb1, ME2, 0, ME2)))/2))/M_PI + 
       (2.*m_alpha*(1./4 + (-Re(B0i(bb0, MM2, 0, MM2)) - 
            Re(B0i(bb1, MM2, 0, MM2)))/2))/M_PI + 
       (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
       (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
       (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
       (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)))/m_t + 
  (64.*m_alpha*ME2*M_PI*m_u*(pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - 
       B0i(bb0, m_t, MM2, MM2) + (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, 
         ME2) + 2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
       2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 - ME*MM - m_t)*
        (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, 
          ME2)) + m_t*(C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, 
          MM2)) + (2.*ME2 + ME*MM)*(C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 
         2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, 
          ME2, ME2)) + MM2*(-2.*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) - 
         3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 0, 
          MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) - 
         3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, 
          MM2, MM2))) - m_alpha*M_PI*
      ((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - 
            Re(B0i(bb1, ME2, 0, ME2)))/2))/M_PI + 
       (2.*m_alpha*(1./4 + (-Re(B0i(bb0, MM2, 0, MM2)) - 
            Re(B0i(bb1, MM2, 0, MM2)))/2))/M_PI + 
       (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
       (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
       (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
       (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)))/pow(m_t,2) + 
  (64.*m_alpha*MM2*M_PI*m_u*(pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - 
       B0i(bb0, m_t, MM2, MM2) + (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, 
         ME2) + 2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
       2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 - ME*MM - m_t)*
        (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, 
          ME2)) + m_t*(C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, 
          MM2)) + (2.*ME2 + ME*MM)*(C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 
         2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, 
          ME2, ME2)) + MM2*(-2.*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) - 
         3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 0, 
          MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) - 
         3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, 
          MM2, MM2))) - m_alpha*M_PI*
      ((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - 
            Re(B0i(bb1, ME2, 0, ME2)))/2))/M_PI + 
       (2.*m_alpha*(1./4 + (-Re(B0i(bb0, MM2, 0, MM2)) - 
            Re(B0i(bb1, MM2, 0, MM2)))/2))/M_PI + 
       (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
       (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
       (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
       (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)))/pow(m_t,2) - 
  (32.*m_alpha*M_PI*pow(m_u,2)*(pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - 
       B0i(bb0, m_t, MM2, MM2) + (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, 
         ME2) + 2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
       2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 - ME*MM - m_t)*
        (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, 
          ME2)) + m_t*(C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, 
          MM2)) + (2.*ME2 + ME*MM)*(C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 
         2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, 
          ME2, ME2)) + MM2*(-2.*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) - 
         3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 0, 
          MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) - 
         3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, 
          MM2, MM2))) - m_alpha*M_PI*
      ((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - 
            Re(B0i(bb1, ME2, 0, ME2)))/2))/M_PI + 
       (2.*m_alpha*(1./4 + (-Re(B0i(bb0, MM2, 0, MM2)) - 
            Re(B0i(bb1, MM2, 0, MM2)))/2))/M_PI + 
       (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
       (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
       (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
       (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)))/pow(m_t,2) + 
  (8.*m_alpha*M_PI*((8.*pow(m_alpha,2)*ME*ME2*MM*(ME2 + ME*MM)*
       (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc11, ME2, m_t, ME2, 0, ME2, 
         ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
        C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, 
         ME2)))/m_t + (4.*pow(m_alpha,2)*pow(ME2,2)*MM2*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
        C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, 
          ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
        C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)))/m_t + 
     (16.*pow(m_alpha,2)*ME*ME2*MM*MM2*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
        C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, 
          ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
        C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)))/m_t + 
     (8.*pow(m_alpha,2)*ME*MM*(ME2 + ME*MM)*MM2*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
        C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, 
          ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
        C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)))/m_t + 
     (4.*pow(m_alpha,2)*ME2*pow(MM2,2)*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
        C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, 
          ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
        C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)))/m_t - 
     (8.*pow(m_alpha,2)*ME*MM*(ME2 + ME*MM)*m_s*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
        C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, 
          ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
        C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)))/m_t - 
     (4.*pow(m_alpha,2)*ME2*MM2*m_s*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
        C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, 
          ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
        C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)))/m_t - 
     (8.*pow(m_alpha,2)*ME2*(pow(MM2,2) + (MM2 - m_t)*(ME2 - m_u))*
       (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc11, ME2, m_t, ME2, 0, ME2, 
         ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
        C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, 
         ME2)))/m_t - 4.*pow(m_alpha,2)*ME*MM*m_u*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
       C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, 
         ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 
        0, ME2, ME2)) + (16.*pow(m_alpha,2)*ME2*pow(MM2,2)*(C0i(cc1, MM2, m_t, MM2, 0, MM2, 
         MM2) + C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 
        2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, 
         MM2) + C0i(cc22, MM2, m_t, MM2, 0, MM2, MM2)))/m_t - 
     4.*pow(m_alpha,2)*MM2*m_u*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
       C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, 
         MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 
        0, MM2, MM2)) + 4.*ME2*MM*
      (-2.*pow(m_alpha,2)*ME*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
         C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, 
           ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
         C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)) + 
       2.*pow(m_alpha,2)*MM*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, 
           MM2, MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc22, MM2, m_t, MM2, 0, MM2, MM2))) + 
     MM*MM2*(-2.*pow(m_alpha,2)*ME*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
         C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, 
           ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
         C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)) + 
       2.*pow(m_alpha,2)*MM*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, 
           MM2, MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + 
         C0i(cc22, MM2, m_t, MM2, 0, MM2, MM2))) - 
     (2.*ME*ME2*MM2*(-2.*pow(m_alpha,2)*ME*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
          C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, 
            ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
          C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)) + 2.*pow(m_alpha,2)*MM*
         (C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 0, MM2, 
           MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) + 
          C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, MM2, 
           MM2))))/m_t - (10*ME2*MM*MM2*
       (-2.*pow(m_alpha,2)*ME*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
          C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, 
            ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
          C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)) + 2.*pow(m_alpha,2)*MM*
         (C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 0, MM2, 
           MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) + 
          C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, MM2, 
           MM2))))/m_t - 
     (2.*ME*pow(MM2,2)*(-2.*pow(m_alpha,2)*ME*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
          C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, 
            ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
          C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)) + 2.*pow(m_alpha,2)*MM*
         (C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 0, MM2, 
           MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) + 
          C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, MM2, 
           MM2))))/m_t - (2.*ME*MM2*(ME2 + MM2 - m_u)*
       (-2.*pow(m_alpha,2)*ME*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
          C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, 
            ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
          C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)) + 2.*pow(m_alpha,2)*MM*
         (C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 0, MM2, 
           MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) + 
          C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, MM2, 
           MM2))))/m_t + 
     (2.*ME*MM2*m_u*(-2.*pow(m_alpha,2)*ME*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
          C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, 
            ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
          C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)) + 2.*pow(m_alpha,2)*MM*
         (C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 0, MM2, 
           MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) + 
          C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, MM2, 
           MM2))))/m_t - (MM*(-2.*(pow(ME2,2) + ME2*(2.*MM2 - m_s)) + MM2*m_s + 
        (2.*ME2 + MM2 - 2.*m_s)*m_u)*(-2.*pow(m_alpha,2)*ME*(C0i(cc1, ME2, m_t, ME2, 0, ME2, 
           ME2) + C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 
          2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, 
           ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)) + 
        2.*pow(m_alpha,2)*MM*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, 
           MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) + 
          C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, MM2, 
           MM2))))/m_t + (12.*pow(m_alpha,2)*ME*ME2*MM*
       (ME2*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc11, ME2, m_t, ME2, 0, 
           ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
          C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, 
           ME2)) + MM2*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
          C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, 
            MM2, MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + 
          C0i(cc22, MM2, m_t, MM2, 0, MM2, MM2))))/m_t + 
     (12.*pow(m_alpha,2)*ME*MM*MM2*(ME2*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
          C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, 
            ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
          C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)) + 
        MM2*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 0, 
           MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) + 
          C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, MM2, 
           MM2))))/m_t - (12.*pow(m_alpha,2)*ME*MM*m_u*
       (ME2*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc11, ME2, m_t, ME2, 0, 
           ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
          C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, 
           ME2)) + MM2*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + 
          C0i(cc11, MM2, m_t, MM2, 0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, 
            MM2, MM2) + C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + 
          C0i(cc22, MM2, m_t, MM2, 0, MM2, MM2))))/m_t - 
     (4.*pow(m_alpha,2)*ME*ME2*MM*((ME2 - 4.*ME*MM)*(C0i(cc1, ME2, m_t, ME2, 0, ME2, 
           ME2) + C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 
          2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, 
           ME2, ME2) + C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)) + 
        MM2*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 0, 
           MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) + 
          C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, MM2, 
           MM2))))/m_t - (4.*pow(m_alpha,2)*ME*MM*MM2*
       ((ME2 - 4.*ME*MM)*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
          C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, 
            ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
          C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)) + 
        MM2*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 0, 
           MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) + 
          C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, MM2, 
           MM2))))/m_t + (4.*pow(m_alpha,2)*ME*MM*m_u*
       ((ME2 - 4.*ME*MM)*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 
          C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, 
            ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
          C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)) + 
        MM2*(C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 0, 
           MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) + 
          C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, MM2, 
           MM2))))/m_t - 4.*ME2*(pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - 
         B0i(bb0, m_t, MM2, MM2) + (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, 
           ME2) + (-2.*MM2 + m_t)*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
         2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
         2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 + ME*MM)*
          C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 2.*ME2*C0i(cc11, ME2, m_t, ME2, 
           0, ME2, ME2) - ME*MM*C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 
         4.*ME2*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) - 
         2.*ME*MM*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) - 
         2.*ME2*C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) - 
         ME*MM*C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
         m_t*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc1, MM2, m_t, MM2, 0, 
            MM2, MM2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
           C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2)) + (2.*ME2 - ME*MM)*
          C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2) + 
         MM2*(-3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 
            0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) - 
           3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, 
            MM2, MM2))) - m_alpha*M_PI*
        ((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - Re(B0i(bb1, ME2, 
                0, ME2)))/2))/M_PI + (2.*m_alpha*(1./4 + 
            (-Re(B0i(bb0, MM2, 0, MM2)) - Re(B0i(bb1, MM2, 0, MM2)))/2))/M_PI + 
         (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
         (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
         (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
         (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)) + 
     (8.*ME2*MM2*(pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - 
          B0i(bb0, m_t, MM2, MM2) + (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, 
            ME2) + (-2.*MM2 + m_t)*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
          2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
          2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 + ME*MM)*
           C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + 2.*ME2*C0i(cc11, ME2, m_t, ME2, 
            0, ME2, ME2) - ME*MM*C0i(cc11, ME2, m_t, ME2, 0, ME2, ME2) + 
          4.*ME2*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) - 
          2.*ME*MM*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) - 
          2.*ME2*C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) - 
          ME*MM*C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
          m_t*(C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc1, MM2, m_t, MM2, 0, 
             MM2, MM2) + C0i(cc2, ME2, m_t, ME2, 0, ME2, ME2) + 
            C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2)) + (2.*ME2 - ME*MM)*
           C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2) + 
          MM2*(-3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 
             0, MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) - 
            3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, 
             MM2, MM2))) - m_alpha*M_PI*
         ((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - Re(B0i(bb1, ME2, 
                 0, ME2)))/2))/M_PI + (2.*m_alpha*(1./4 + 
             (-Re(B0i(bb0, MM2, 0, MM2)) - Re(B0i(bb1, MM2, 0, MM2)))/2))/
           M_PI + (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
          (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
          (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
          (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)))/m_t - 
     (4.*pow(ME2,2)*(pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - 
          B0i(bb0, m_t, MM2, MM2) + (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, 
            ME2) + 2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
          2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 - ME*MM - m_t)*
           (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, 
             ME2, ME2)) + m_t*(C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
            C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc2, MM2, m_t, MM2, 0, 
             MM2, MM2)) + (2.*ME2 + ME*MM)*(C0i(cc11, ME2, m_t, ME2, 0, ME2, 
             ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
            C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)) + 
          MM2*(-2.*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) - 
            3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 0, 
             MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) - 
            3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, 
             MM2, MM2))) - m_alpha*M_PI*
         ((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - Re(B0i(bb1, ME2, 
                 0, ME2)))/2))/M_PI + (2.*m_alpha*(1./4 + 
             (-Re(B0i(bb0, MM2, 0, MM2)) - Re(B0i(bb1, MM2, 0, MM2)))/2))/
           M_PI + (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
          (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
          (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
          (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)))/m_t - 
     (24.*ME2*MM2*(pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - 
          B0i(bb0, m_t, MM2, MM2) + (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, 
            ME2) + 2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
          2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 - ME*MM - m_t)*
           (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, 
             ME2, ME2)) + m_t*(C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
            C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc2, MM2, m_t, MM2, 0, 
             MM2, MM2)) + (2.*ME2 + ME*MM)*(C0i(cc11, ME2, m_t, ME2, 0, ME2, 
             ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
            C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)) + 
          MM2*(-2.*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) - 
            3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 0, 
             MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) - 
            3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, 
             MM2, MM2))) - m_alpha*M_PI*
         ((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - Re(B0i(bb1, ME2, 
                 0, ME2)))/2))/M_PI + (2.*m_alpha*(1./4 + 
             (-Re(B0i(bb0, MM2, 0, MM2)) - Re(B0i(bb1, MM2, 0, MM2)))/2))/
           M_PI + (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
          (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
          (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
          (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)))/m_t - 
     (4.*pow(MM2,2)*(pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - 
          B0i(bb0, m_t, MM2, MM2) + (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, 
            ME2) + 2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
          2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 - ME*MM - m_t)*
           (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, 
             ME2, ME2)) + m_t*(C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
            C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc2, MM2, m_t, MM2, 0, 
             MM2, MM2)) + (2.*ME2 + ME*MM)*(C0i(cc11, ME2, m_t, ME2, 0, ME2, 
             ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
            C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)) + 
          MM2*(-2.*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) - 
            3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 0, 
             MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) - 
            3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, 
             MM2, MM2))) - m_alpha*M_PI*
         ((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - Re(B0i(bb1, ME2, 
                 0, ME2)))/2))/M_PI + (2.*m_alpha*(1./4 + 
             (-Re(B0i(bb0, MM2, 0, MM2)) - Re(B0i(bb1, MM2, 0, MM2)))/2))/
           M_PI + (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
          (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
          (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
          (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)))/m_t + 
     (8.*ME2*m_u*(pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - 
          B0i(bb0, m_t, MM2, MM2) + (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, 
            ME2) + 2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
          2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 - ME*MM - m_t)*
           (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, 
             ME2, ME2)) + m_t*(C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
            C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc2, MM2, m_t, MM2, 0, 
             MM2, MM2)) + (2.*ME2 + ME*MM)*(C0i(cc11, ME2, m_t, ME2, 0, ME2, 
             ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
            C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)) + 
          MM2*(-2.*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) - 
            3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 0, 
             MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) - 
            3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, 
             MM2, MM2))) - m_alpha*M_PI*
         ((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - Re(B0i(bb1, ME2, 
                 0, ME2)))/2))/M_PI + (2.*m_alpha*(1./4 + 
             (-Re(B0i(bb0, MM2, 0, MM2)) - Re(B0i(bb1, MM2, 0, MM2)))/2))/
           M_PI + (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
          (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
          (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
          (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)))/m_t + 
     (8.*MM2*m_u*(pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - 
          B0i(bb0, m_t, MM2, MM2) + (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, 
            ME2) + 2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
          2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 - ME*MM - m_t)*
           (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, 
             ME2, ME2)) + m_t*(C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
            C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc2, MM2, m_t, MM2, 0, 
             MM2, MM2)) + (2.*ME2 + ME*MM)*(C0i(cc11, ME2, m_t, ME2, 0, ME2, 
             ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
            C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)) + 
          MM2*(-2.*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) - 
            3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 0, 
             MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) - 
            3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, 
             MM2, MM2))) - m_alpha*M_PI*
         ((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - Re(B0i(bb1, ME2, 
                 0, ME2)))/2))/M_PI + (2.*m_alpha*(1./4 + 
             (-Re(B0i(bb0, MM2, 0, MM2)) - Re(B0i(bb1, MM2, 0, MM2)))/2))/
           M_PI + (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
          (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
          (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
          (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)))/m_t - 
     (4.*pow(m_u,2)*(pow(m_alpha,2)*1. + pow(m_alpha,2)*(-B0i(bb0, m_t, ME2, ME2) - 
          B0i(bb0, m_t, MM2, MM2) + (-2.*ME2 + m_t)*C0i(cc0, ME2, m_t, ME2, 0, ME2, 
            ME2) + 2.*C0i(cc00, ME2, m_t, ME2, 0, ME2, ME2) + 
          2.*C0i(cc00, MM2, m_t, MM2, 0, MM2, MM2) - (2.*ME2 - ME*MM - m_t)*
           (C0i(cc1, ME2, m_t, ME2, 0, ME2, ME2) + C0i(cc2, ME2, m_t, ME2, 0, 
             ME2, ME2)) + m_t*(C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) + 
            C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc2, MM2, m_t, MM2, 0, 
             MM2, MM2)) + (2.*ME2 + ME*MM)*(C0i(cc11, ME2, m_t, ME2, 0, ME2, 
             ME2) + 2.*C0i(cc12, ME2, m_t, ME2, 0, ME2, ME2) + 
            C0i(cc22, ME2, m_t, ME2, 0, ME2, ME2)) + 
          MM2*(-2.*C0i(cc0, MM2, m_t, MM2, 0, MM2, MM2) - 
            3.*C0i(cc1, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc11, MM2, m_t, MM2, 0, 
             MM2, MM2) + 2.*C0i(cc12, MM2, m_t, MM2, 0, MM2, MM2) - 
            3.*C0i(cc2, MM2, m_t, MM2, 0, MM2, MM2) + C0i(cc22, MM2, m_t, MM2, 0, 
             MM2, MM2))) - m_alpha*M_PI*
         ((2.*m_alpha*(1./4 + (-Re(B0i(bb0, ME2, 0, ME2)) - Re(B0i(bb1, ME2, 
                 0, ME2)))/2))/M_PI + (2.*m_alpha*(1./4 + 
             (-Re(B0i(bb0, MM2, 0, MM2)) - Re(B0i(bb1, MM2, 0, MM2)))/2))/
           M_PI + (2.*m_alpha*ME2*Re(B0i(dbb0, ME2, 0, ME2)))/M_PI + 
          (2.*m_alpha*MM2*Re(B0i(dbb0, MM2, 0, MM2)))/M_PI - 
          (2.*m_alpha*ME2*Re(B0i(dbb1, ME2, 0, ME2)))/M_PI - 
          (2.*m_alpha*MM2*Re(B0i(dbb1, MM2, 0, MM2)))/M_PI)))/m_t))/m_t;

}


DECLARE_VIRTUALME2_GETTER(EXTRAXS::emu_emu_Virtual,"emu_emu_Virtual")
Virtual_ME2_Base *ATOOLS::Getter
<PHASIC::Virtual_ME2_Base,PHASIC::Process_Info,EXTRAXS::emu_emu_Virtual>::
operator()(const Process_Info &pi) const
{
  if (pi.m_loopgenerator.find("Internal")!=0) return NULL;
  // if (pi.m_fi.m_nlotype==nlo_type::loop) {
  Flavour_Vector fl(pi.ExtractFlavours());
  if (fl.size()!=4) return NULL;
  if ( ( fl[0].Kfcode() == Flavour(kf_e) && fl[1].Kfcode() == Flavour(kf_mu)  ) || 
       ( fl[0].Kfcode() == Flavour(kf_mu) && fl[1].Kfcode() == Flavour(kf_e)  )
    )
  {
    return new emu_emu_Virtual(pi,fl,0.0,0);
  }
  return NULL;
}
