#include "ATOOLS/Math/Vec4.H"
#include "ATOOLS/Math/Vector.H"
#include "MODEL/Main/Coupling_Data.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"

namespace ATOOLS {
  class NLO_subevt;
}

namespace PHASIC {

  class EWVirtKFactor_Setter : public PHASIC::KFactor_Setter_Base {
  private:
    ATOOLS::Vec4D_Vector m_p;
    PHASIC::Virtual_ME2_Base* p_ewloop;
    MODEL::Coupling_Map m_cpls;
    double m_deltaew;

    void CopyMomenta();
    void CopyMomenta(const ATOOLS::NLO_subevt& evt);

    void InitEWVirt();
    void CalcEWCorrection();

  public:

    EWVirtKFactor_Setter(const PHASIC::KFactor_Setter_Arguments &args);
    ~EWVirtKFactor_Setter();

    // Default KFactor method
    double KFactor(const int mode=0);
    // KFactor for Comix, not yet tested or validated
    double KFactor(const ATOOLS::NLO_subevt& evt);
  };
}

#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Phys/Flavour.H"

#include "MODEL/Main/Model_Base.H"

#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"


using namespace PHASIC;
using namespace ATOOLS;

EWVirtKFactor_Setter::EWVirtKFactor_Setter
(const KFactor_Setter_Arguments &args):
  KFactor_Setter_Base(args), p_ewloop(NULL), m_deltaew(0.)
{
  DEBUG_FUNC("");
  InitEWVirt();
}

EWVirtKFactor_Setter::~EWVirtKFactor_Setter()
{
  if (p_ewloop) { delete p_ewloop; p_ewloop=NULL; }
  for (const auto& kv : m_cpls)
    delete kv.second;
}

double EWVirtKFactor_Setter::KFactor(const int mode)
{
  if(!m_on) return 1.;
  CopyMomenta();
  CalcEWCorrection();
  return (1.+m_deltaew);
}

double EWVirtKFactor_Setter::KFactor(const NLO_subevt& evt)
{
  if(!m_on) return 1.;
  CopyMomenta(evt);
  CalcEWCorrection();
  return (1.+m_deltaew);
}

void EWVirtKFactor_Setter::CopyMomenta()
{
  m_p=p_proc->Integrator()->Momenta();
}

void EWVirtKFactor_Setter::CopyMomenta(const NLO_subevt& evt)
{
  Vec4D_Vector p(evt.p_mom, &(evt.p_mom[evt.m_n]));
  for (size_t i(0);i<p_proc->NIn();++i) p[i]=-p[i];
}

void EWVirtKFactor_Setter::InitEWVirt()
{
  Process_Info loop_pi(p_proc->Info());
  loop_pi.m_fi.m_nlotype=nlo_type::loop;
  loop_pi.m_maxcpl=p_proc->MaxOrders();
  loop_pi.m_mincpl=p_proc->MinOrders();
  for (size_t i(0);i<loop_pi.m_fi.m_nlocpl.size();++i)
    loop_pi.m_fi.m_nlocpl[i]=(i==1?1:0);
  for (size_t i(0);i<Min(loop_pi.m_maxcpl.size(),
                         p_proc->Info().m_fi.m_nlocpl.size());++i) {
    loop_pi.m_maxcpl[i]+=loop_pi.m_fi.m_nlocpl[i]-p_proc->Info().m_fi.m_nlocpl[i];
    loop_pi.m_mincpl[i]+=loop_pi.m_fi.m_nlocpl[i]-p_proc->Info().m_fi.m_nlocpl[i];
  }
  msg_Debugging()<<"Load "<<loop_pi.m_loopgenerator
                 <<" process for "<<p_proc->Name()
                 <<" of order "<<loop_pi.m_mincpl<<" .. "<<loop_pi.m_maxcpl
                 <<std::endl;
  p_ewloop=PHASIC::Virtual_ME2_Base::GetME2(loop_pi);
  if (!p_ewloop) {
    THROW(not_implemented,"Couldn't find EW Virtual for "
                          +p_proc->Name());
  }
  MODEL::s_model->GetCouplings(m_cpls);
  p_ewloop->SetCouplings(m_cpls);
  p_ewloop->SetSubType(sbt::qed);
}

void EWVirtKFactor_Setter::CalcEWCorrection()
{
  DEBUG_FUNC("");
  m_deltaew=0.;
  p_ewloop->SetRenScale(p_proc->ScaleSetter()->Scale(stp::ren,1));
  p_ewloop->Calc(m_p);
  // OL returns V/(as/2pi*B)
  double fac(p_ewloop->AlphaQED()/2./M_PI);
  double B(p_ewloop->ME_Born());
  double V(p_ewloop->ME_Finite()*fac*B);
  m_deltaew=V/B;
  if (msg_LevelIsDebugging()) {
    msg_Out()<<"p_T    = "<<(m_p[2]+m_p[3]).PPerp()<<std::endl;
    msg_Out()<<"\\mu_R  = "<<p_proc->ScaleSetter()->Scale(stp::ren,1)<<std::endl;
    msg_Out()<<"cpl    = "<<fac<<std::endl;
    msg_Out()<<"VI_e2  = "<<p_ewloop->ME_E2()*fac*B<<std::endl;
    msg_Out()<<"VI_e1  = "<<p_ewloop->ME_E1()*fac*B<<std::endl;
    msg_Out()<<"VI_fin = "<<V<<std::endl;
    msg_Out()<<"B      = "<<B<<std::endl;
    msg_Out()<<"\\delta = "<<m_deltaew<<std::endl;
  }
}


DECLARE_GETTER(PHASIC::EWVirtKFactor_Setter,"EWVirt",
               PHASIC::KFactor_Setter_Base,
               PHASIC::KFactor_Setter_Arguments);

KFactor_Setter_Base *ATOOLS::Getter<
  PHASIC::KFactor_Setter_Base,
  PHASIC::KFactor_Setter_Arguments,
  PHASIC::EWVirtKFactor_Setter>::operator()
(const KFactor_Setter_Arguments &args) const
{
  msg_Info()<<"Loading EWVirt KFactor for "<<args.p_proc->Name()<<std::endl;
  return new EWVirtKFactor_Setter(args);
}

void ATOOLS::Getter<
  PHASIC::KFactor_Setter_Base,
  PHASIC::KFactor_Setter_Arguments,
  PHASIC::EWVirtKFactor_Setter>::PrintInfo(std::ostream &str,
                                           const size_t width) const
{
  str<<"EWVirt K-Factor\n";
}
