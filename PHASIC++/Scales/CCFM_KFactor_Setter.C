#include "PHASIC++/Scales/KFactor_Setter_Base.H"

#include "PHASIC++/Scales/MEPS_Scale_Setter.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Process/Process_Base.H"
#include "MODEL/Main/Running_AlphaS.H"

#include <map>

namespace PHASIC {

  class CCFM_Sudakov {
  private:

    ATOOLS::Flavour m_fl;
    MODEL::Running_AlphaS *p_as;

    double m_Q2, m_mur2;
    int m_fo;

  public:

    CCFM_Sudakov(const ATOOLS::Flavour &fl);

    double Delta(const double &z,const double &q2,const double &Q2);
    double Delta1(const double &z,const double &q2,const double &Q2,const double &mur2);

  };// end of class CCFM_Sudakov

  class CCFM_KFactor_Setter: public KFactor_Setter_Base {
  private:

    MEPS_Scale_Setter *p_meps;

    std::map<ATOOLS::Flavour,CCFM_Sudakov*> m_suds;

    double KFactor(const ATOOLS::Cluster_Amplitude *ampl,
		   const double &mur2,const int mode);

  public:

    CCFM_KFactor_Setter(const KFactor_Setter_Arguments &args);

    ~CCFM_KFactor_Setter();

    double KFactor(const int mode=0);
    double KFactor(const ATOOLS::NLO_subevt& evt);

  };// end of class CCFM_KFactor_Setter

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(CCFM_KFactor_Setter,"CCFM",
	       KFactor_Setter_Base,KFactor_Setter_Arguments);

KFactor_Setter_Base *ATOOLS::Getter
<KFactor_Setter_Base,KFactor_Setter_Arguments,CCFM_KFactor_Setter>::
operator()(const KFactor_Setter_Arguments &args) const
{
  return new CCFM_KFactor_Setter(args);
}

void ATOOLS::Getter<KFactor_Setter_Base,KFactor_Setter_Arguments,
		    CCFM_KFactor_Setter>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"CCFM kfactor scheme\n";
}

CCFM_KFactor_Setter::CCFM_KFactor_Setter
(const KFactor_Setter_Arguments &args):
  KFactor_Setter_Base(args)
{
  p_meps=dynamic_cast<MEPS_Scale_Setter*>(p_proc->ScaleSetter());
  if (p_meps==NULL) THROW(fatal_error,"Must use METS scale");
  m_suds[Flavour(kf_gluon)] = new CCFM_Sudakov(Flavour(kf_gluon));
  for (size_t i(0);i<=6;++i) {
    m_suds[Flavour(i,0)] = new CCFM_Sudakov(Flavour(i,0));
    m_suds[Flavour(i,1)] = new CCFM_Sudakov(Flavour(i,1));
  }
}

CCFM_KFactor_Setter::~CCFM_KFactor_Setter()
{
  for (std::map<ATOOLS::Flavour,CCFM_Sudakov*>::const_iterator
	 sit(m_suds.begin());sit!=m_suds.end();++sit)
    delete sit->second;
}

double CCFM_KFactor_Setter::KFactor(const ATOOLS::NLO_subevt &evt)
{
  if (!m_on || evt.m_me==0.0 || evt.p_ampl==NULL) return m_weight=1.0;
  DEBUG_FUNC(p_proc->Name());
  msg_Debugging()<<evt<<"\n";
  return KFactor(evt.p_ampl,evt.m_mu2[stp::ren],0);
}

double CCFM_KFactor_Setter::KFactor(const int mode)
{
  if (!m_on || p_meps->Amplitudes().empty()) return m_weight=1.0;
  Cluster_Amplitude *ampl(p_meps->Amplitudes().back());
  if (p_proc->Info().Has(nlo_type::real)) {
    if (p_proc->GetSubevtList()->back()->p_ampl &&
	p_proc->GetSubevtList()->back()->p_ampl->Next())
      ampl=p_proc->GetSubevtList()->back()->p_ampl->Next();
    else {
      msg_Debugging()<<"Real emission amplitude not found\n";
      return m_weight=1.0;
    }
  }
  DEBUG_FUNC(p_proc->Name()<<", mode = "<<mode);
  double mur2(p_proc->ScaleSetter()->Scale(stp::ren));
  return KFactor(ampl,mur2,mode);
}

double CCFM_KFactor_Setter::KFactor
(const ATOOLS::Cluster_Amplitude *ampl,const double &muR2,const int mode)
{
  m_weight=1.0;
  msg_Debugging()<<"\\mu_R = "<<sqrt(muR2)<<"\n";
  double sub(0.0);
  Vec4D_Vector p(ampl->Legs().size());
  for (size_t i(0);i<p.size();++i) p[i]=ampl->Leg(i)->Mom();
  for (;ampl->Next();ampl=ampl->Next()) {
    Cluster_Amplitude *next(ampl->Next());
    msg_Debugging()<<*ampl<<"\n";
    for (size_t i(0);i<next->Legs().size();++i) {
      Cluster_Leg *l(next->Leg(i));
      int is((l->Id()&3)?1:0);
      if (l->K()==0) continue;
      Cluster_Leg *lj(ampl->IdLeg(ampl->IdNew()));
      Vec4D qcur(lj->Mom()), Qold(p[ID(l->Id()).front()]), Qnew(Qold+qcur);
      for (size_t j(0), k(0);k<p.size();++j) {
	if (ampl->Leg(j)==lj) {
	  msg_Debugging()<<"<- p_"<<i<<" = "<<p[i]<<"\n";
	  msg_Debugging()<<"<- p_"<<k<<" = "<<p[k]<<"\n";
	  p[i]+=p[k];
	  msg_Debugging()<<"-> p_"<<i<<" = "<<p[i]<<"\n";
	  p.erase(p.begin()+k);
	  continue;
	}
	++k;
      }
      double z(-Qnew.PPlus()>-Qnew.PMinus()?
	       Qnew.PPlus()/Qold.PPlus():
	       Qnew.PMinus()/Qold.PMinus());
      std::map<ATOOLS::Flavour,CCFM_Sudakov*>::iterator sit(m_suds.find(l->Flav()));
      if (!is || sit==m_suds.end()) continue;
      double gamma[2]={0.0,0.0};
      if ((mode&1) && p_proc->Info().m_fi.m_nlotype!=nlo_type::lo) {
       	gamma[0]=sit->second->Delta1(z,qcur.PPerp2(),Qnew.PPerp2(),muR2);
      }
      double f=1./z/(1./z+1./(1.-z));
      double delta=sit->second->Delta(z,qcur.PPerp2(),Qnew.PPerp2());
      msg_Debugging()<<"Sudakov for "<<ID(l->Id())<<"["<<is
		     <<"] -> \\Delta_{"<<l->Flav()<<"}("<<z<<","<<qcur.PPerp()
		     <<","<<Qnew.PPerp()<<") = "<<delta<<", \\Gamma = "
		     <<gamma[0]<<", f = "<<f<<"\n";
      m_weight*=(1.-f)+f*delta*(1.0-gamma[0]);
    }
  }
  msg_Debugging()<<*ampl<<"\n";
  msg_Debugging()<<"w = "<<m_weight<<"\n";
  return m_weight;
}

CCFM_Sudakov::CCFM_Sudakov(const ATOOLS::Flavour &fl):
  m_fl(fl), p_as(MODEL::as)
{
}

double CCFM_Sudakov::Delta(const double &z,const double &q2,const double &Q2)
{
  double as2pi((*p_as)(q2)/(2.0*M_PI));
  if (m_fl.IsGluon()) {
    return exp(-as2pi*3.0*(2.0*log(1.0/z)*log(Q2/q2/z)));
  }
  return 1;
}

double CCFM_Sudakov::Delta1(const double &z,const double &q2,const double &Q2,const double &mur2)
{
  double as2pi((*p_as)(mur2)/(2.0*M_PI));
  if (m_fl.IsGluon()) {
    return -as2pi*3.0*(2.0*log(1.0/z)*log(Q2/q2/z));
  }
  return 0;
}
