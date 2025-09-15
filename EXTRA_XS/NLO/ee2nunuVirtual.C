#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "clooptools.h"

using namespace PHASIC;
using namespace ATOOLS;
using namespace MODEL;

namespace EXTRAXS {
  class ee2nunuVirtual :  public Virtual_ME2_Base {
  public:
    ee2nunuVirtual(const Process_Info& pi, const Flavour_Vector& flavs,
                  const double& ep2, const double& ep);
    ~ee2nunuVirtual(){};
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    Complex Full();
    double m_me, m_mpi, m_pij[5][5], m_s45;
    double m_betapi, m_s, m_sp;
    double ME2, MP2, m_alpha;
    double S34,T24,T14,S12;
    double S,U,T;
    double SW2, SW, CW, CW2, MZ2;
    Vec4D_Vector m_pmom;
  };
}

using namespace EXTRAXS;


ee2nunuVirtual::ee2nunuVirtual(const Process_Info& pi, const Flavour_Vector& flavs,
                  const double& ep2, const double& ep) :
                  Virtual_ME2_Base(pi, flavs)
{
  PRINT_INFO("initialised XS_ee2nunuVirtual");
  // Flavour_Vector outflavs = args.m_outflavs;
  // p_bornme = dynamic_cast<ME2_Base*>(PHASIC::Tree_ME2_Base::GetME2(bargs));
  // if (!p_bornme) THROW(fatal_error,"no born me found.");
  // const Flavour_Vector fl = args.Flavours();
  double ME = flavs[0].Mass();
  double MP = flavs[3].Mass();
  ME2 = ME*ME;
  MP2 = MP*MP;
  m_mode=0;
  // m_oqcd = 1;
  // m_oew  = 2;
  // m_me  = outflavs[0].Mass();
  // m_mpi = outflavs[2].Mass();
  // m_s = sqr(rpa->gen.Ecms());
  // m_betapi = sqrt(1-m_mpi*m_mpi/m_s);
  SW2 = MODEL::s_model->ComplexConstant("csin2_thetaW").real();
  CW2 = 1.-SW*SW;
  MZ2 = Flavour(23).Mass()*Flavour(23).Mass();
  CW = sqrt(CW2);
  SW = sqrt(SW2);
}


void ee2nunuVirtual::Calc(const ATOOLS::Vec4D_Vector& p)
{
  m_pmom = p;
  m_sp = (p[2]+p[3]).Abs2();
  m_s45 = (p[3]+p[4]).Abs2();
  m_betapi = sqrt(1-m_mpi*m_mpi/m_sp);
  S = m_sp;
  S34 = (p[2]+p[3]).Abs2();
  T14 = (p[0]-p[3]).Abs2();
  T24 = (p[1]-p[3]).Abs2();
  T=(p[0]-p[2]).Abs2();;
  U=(p[0]-p[3]).Abs2();
  m_alpha = (*aqed)(m_sp);
  // PRINT_VAR(Full());
  double res;// = Initial()+Final()+IFI();
  double Q2 = S;
  double s15 = (p[0]+p[4]).Abs2(); 
  double s25 = (p[1]+p[4]).Abs2(); 
  double s35 = (p[2]+p[4]).Abs2();

  Setlambda(0.);
  m_res.Finite() = Full().real();
  Setlambda(-1);
  m_res.IR() = Full().real();

}


Complex ee2nunuVirtual::Full(){
  return -1/16.*(M_PI*(pow(1 - 4*SW2,2)*(2.*ME2 - T - U)*(ME2*ME2 - T*U)*
        (-A0(ME2) + ME2*(-1. - B0i(bb0,1,ME2,ME2) + 2.*B0i(bb0,ME2,0,ME2))) + 
       2.*SW2*((1. - 4*SW2)*(2.*ME2 - T - U)*(ME2*ME2 - T*U)*
           (A0(ME2) + ME2*(1. + B0i(bb0,1,ME2,ME2) - 2.*B0i(bb0,ME2,0,ME2))) \
+ 2.*pow(ME2 - T,2)*(2.*(-1 + 4*ME2)*SW2 + 
             (1 - 4*ME2)*SW2*B0i(bb0,1,ME2,ME2) + 
             2.*(ME2 + 2.*SW2 - 8*ME2*SW2)*(2.*ME2 - T - U)*
              (-B0i(bb0,1,ME2,ME2) + B0i(bb0,ME2,0,ME2)) + 
             4.*(1. - 4.*ME2)*ME2*SW2*C0i(cc0,ME2,1,ME2,0,ME2,ME2) + 
             2.*(-1. + 4*ME2)*SW2*C0i(cc0,ME2,1,ME2,0,ME2,ME2)) - 
          ME2*(4*(1. - 2.*SW2 + ME2*(-3. + 8.*SW2))*(2.*ME2 - T - U)*
              (-B0i(bb0,1,ME2,ME2) + B0i(bb0,ME2,0,ME2)) + 
             (-1. + 4*ME2)*(1 - 2.*SW2)*
              (2. - B0i(bb0,1,ME2,ME2) + 
                2.*(1 - 2.*ME2)*C0i(cc0,ME2,1,ME2,0,ME2,ME2)))) - 
       (1 - 2.*SW2)*((1 - 4*SW2)*(2.*ME2 - T - U)*(ME2*ME2 - T*U)*
           (A0(ME2) + ME2*(1. + B0i(bb0,1,ME2,ME2) - 2.*B0i(bb0,ME2,0,ME2))) + 
          2.*ME2*(2.*(-1. + 4*ME2)*SW2 + (1 - 4*ME2)*SW2*B0i(bb0,1,ME2,ME2) + 
             2.*(ME2 + 2.*SW2 - 8*ME2*SW2)*(2.*ME2 - T - U)*
              (-B0i(bb0,1,ME2,ME2) + B0i(bb0,ME2,0,ME2)) + 
             4*(1 - 4*ME2)*ME2*SW2*C0i(cc0,ME2,1,ME2,0,ME2,ME2) + 
             2.*(-1 + 4*ME2)*SW2*C0i(cc0,ME2,1,ME2,0,ME2,ME2)) - 
          pow(ME2 - U,2)*(4*(1 - 2.*SW2 + ME2*(-3 + 8*SW2))*
              (2.*ME2 - T - U)*(-B0i(bb0,1,ME2,ME2) + B0i(bb0,ME2,0,ME2)) + 
             (-1. + 4*ME2)*(1. - 2.*SW2)*
              (2. - B0i(bb0,1,ME2,ME2) + 
                2.*(1 - 2.*ME2)*C0i(cc0,ME2,1,ME2,0,ME2,ME2))))))/
   (CW2*(-1 + 4*ME2)*pow(-1 + MZ2,2)*(SW2*SW2));
}


DECLARE_VIRTUALME2_GETTER(EXTRAXS::ee2nunuVirtual,"ee2nunuVirtual")
Virtual_ME2_Base *ATOOLS::Getter
<PHASIC::Virtual_ME2_Base,PHASIC::Process_Info,EXTRAXS::ee2nunuVirtual>::
operator()(const Process_Info &pi) const
{
  return NULL;
  const Flavour_Vector fl = pi.ExtractFlavours();
  if(fl.size()!=4) return NULL;
   if ( (fl[0]==Flavour(kf_e) || fl[1]==Flavour(kf_e))  && fl[1]==fl[0].Bar() &&
      (fl[2].IsNeutrino() && fl[2].IsNeutrino()) && fl[3]==fl[2].Bar()) return new ee2nunuVirtual(pi,fl,0.0,0);
    return NULL;
}
