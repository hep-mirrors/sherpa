#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include "PHASIC++/Scales/Tag_Setter.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "PDF/Main/Shower_Base.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

namespace PHASIC {

  class METS_Scale_Setter: public Scale_Setter_Base {
  private:

    std::string m_muf2tag, m_mur2tag;

    ATOOLS::Algebra_Interpreter m_muf2calc, m_mur2calc;

    Tag_Setter m_muf2tagset, m_mur2tagset;

    ATOOLS::Vec4D_Vector   m_p;
    ATOOLS::Flavour_Vector m_f;

  public:

    METS_Scale_Setter(Process_Base *const proc,
		     const std::string &mur2tag,
		     const std::string &muf2tag);

    double CalculateScale(const std::vector<ATOOLS::Vec4D> &p);
    double CalculateScale2(const std::vector<ATOOLS::Vec4D> &p);

    ATOOLS::Vec4D Momentum(const size_t &i) const;

    void SetScale(const std::string &mu2tag,Tag_Setter &mu2tagset,
		  ATOOLS::Algebra_Interpreter &mu2calc);

  };// end of class Scale_Setter_Base

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(METS_Scale_Setter_Getter,"METS",
	       Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *METS_Scale_Setter_Getter::
operator()(const Scale_Setter_Arguments &args) const
{
  return new METS_Scale_Setter(args.p_proc,args.m_ren,args.m_fac);
}

void METS_Scale_Setter_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"mets scale scheme\n";
}

METS_Scale_Setter::METS_Scale_Setter
(Process_Base *const proc,
 const std::string &mur2tag,const std::string &muf2tag): 
  Scale_Setter_Base(proc,true), m_muf2tagset(this), m_mur2tagset(this)
{
  m_p.resize(4);
  SetScale(muf2tag,m_muf2tagset,m_muf2calc);
  SetScale(mur2tag,m_mur2tagset,m_mur2calc);
  m_f=p_proc->Flavours();
  for (size_t i(0);i<p_proc->NIn();++i) m_f[i]=m_f[i].Bar();
}

Vec4D METS_Scale_Setter::Momentum(const size_t &i) const
{
  if (i>m_p.size()) THROW(fatal_error,"Momentum index too large");
  return m_p[i];
}

double METS_Scale_Setter::CalculateScale2(const std::vector<ATOOLS::Vec4D> &momenta) 
{
  p_proc->Integrator()->SwapInOrder();
  double muf2(CalculateScale(momenta));
  p_proc->Integrator()->RestoreInOrder();
  return muf2;
}

double METS_Scale_Setter::CalculateScale
(const std::vector<ATOOLS::Vec4D> &momenta)
{
  DEBUG_FUNC("");
  if (!m_kfkey.Assigned()) {
    std::string kfinfo("O(QCD)="+ToString(p_proc->OrderQCD()));
    msg_Debugging()<<"Assign '"<<p_proc->Name()
		   <<"' '"<<kfinfo<<"'\n";
    m_kfkey.Assign(p_proc->Name(),3,0,p_proc->
		   Integrator()->PSHandler()->GetInfo());
    m_kfkey.SetInfo(kfinfo);
  }
  p_proc->Integrator()->SetMomenta(momenta);
  p_proc->Generator()->SetClusterDefinitions
    (p_proc->Shower()->GetClusterDefinitions());
  Cluster_Amplitude *ampl
    (p_proc->Generator()->ClusterConfiguration(p_proc));
  double kt2max(0.0);
  while (ampl->Next()) {
    kt2max=Max(kt2max,ampl->KT2QCD());
    ampl=ampl->Next();
  }
  msg_Debugging()<<"Core = "<<*ampl<<"\n";
  m_p.resize(ampl->Legs().size());
  ColorID c[4]={ampl->Leg(0)->Col(),ampl->Leg(1)->Col(),
		ampl->Leg(2)->Col(),ampl->Leg(3)->Col()};
  size_t qcd(0);
  for (size_t i(0);i<m_p.size();++i) {
    m_p[i]=ampl->Leg(i)->Mom();
    if (c[i].m_i>0 || c[i].m_j>0) qcd+=1<<i;
  }
  double kt2cmin(std::numeric_limits<double>::max());
  if (qcd!=15) {
    if ((c[0].m_i>0 && c[0].m_i==c[1].m_j) ||
	(c[0].m_j>0 && c[0].m_j==c[1].m_i) ||
	(c[2].m_i>0 && c[2].m_i==c[3].m_j) ||
	(c[2].m_j>0 && c[2].m_j==c[3].m_i)) {
      kt2cmin=Min(kt2cmin,(m_p[0]+m_p[1]).Abs2());
    }
    if ((c[0].m_i>0 && c[0].m_i==c[2].m_j) ||
	(c[0].m_j>0 && c[0].m_j==c[2].m_i) ||
	(c[1].m_i>0 && c[1].m_i==c[3].m_j) ||
	(c[1].m_j>0 && c[1].m_j==c[3].m_i)) {
      kt2cmin=Min(kt2cmin,dabs((m_p[0]+m_p[2]).Abs2()));
    }
    if ((c[0].m_i>0 && c[0].m_i==c[3].m_j) ||
	(c[0].m_j>0 && c[0].m_j==c[3].m_i) ||
	(c[1].m_i>0 && c[1].m_i==c[2].m_j) ||
	(c[1].m_j>0 && c[1].m_j==c[2].m_i)) {
      kt2cmin=Min(kt2cmin,dabs((m_p[0]+m_p[3]).Abs2()));
    }
  }
  if (kt2cmin==std::numeric_limits<double>::max()) {
    kt2cmin=1.0/(m_p[0]+m_p[1]).Abs2();
    kt2cmin+=1.0/(m_p[0]+m_p[2]).Abs2();
    kt2cmin+=1.0/(m_p[0]+m_p[3]).Abs2();
    kt2cmin=-1.0/kt2cmin;
  }
  m_scale[stp::ren]=m_scale[stp::fac]=Max(kt2max,kt2cmin);
  msg_Debugging()<<"QCD scale = "<<sqrt(m_scale[stp::ren])<<"\n";
  m_scale[stp::ren]=m_mur2calc.Calculate()->Get<double>();
  m_scale[stp::fac]=m_muf2calc.Calculate()->Get<double>();
  msg_Debugging()<<"Set \\mu_r = "
		 <<sqrt(m_scale[stp::ren])<<", \\mu_f = "
		 <<sqrt(m_scale[stp::fac])<<"\n";
  ampl->SetMuF2(m_scale[stp::fac]); 
  ampl->SetMuR2(m_scale[stp::ren]);
  while (ampl->Prev()) {
    ampl=ampl->Prev();
    ampl->SetMuF2(m_scale[stp::fac]); 
    ampl->SetMuR2(m_scale[stp::ren]);
  }
  m_kfkey[0]=m_scale[stp::ren];
  m_kfkey[2]=m_kfkey[1]=m_scale[stp::fac];
  return m_scale[stp::fac];
}

void METS_Scale_Setter::SetScale
(const std::string &mu2tag,Tag_Setter &mu2tagset,Algebra_Interpreter &mu2calc)
{ 
  if (mu2tag=="" || mu2tag=="0") THROW(fatal_error,"No scale specified");
  msg_Debugging()<<METHOD<<"(): scale '"<<mu2tag
		 <<"' in '"<<p_proc->Name()<<"' {\n";
  msg_Indent();
  mu2tagset.SetCalculator(&mu2calc);
  mu2calc.SetTagReplacer(&mu2tagset);
  mu2calc.AddTag("MU_F","1.0");
  mu2calc.AddTag("MU_R","1.0");
  mu2calc.AddTag("H_T","1.0");
  mu2calc.AddTag("Q2_CUT","1.0");
  mu2calc.AddTag("Q2_MIN","1.0");
  Process_Integrator *ib(p_proc->Integrator());
  for (size_t i=0;i<ib->NIn()+ib->NOut();++i) 
    mu2calc.AddTag("p["+ToString(i)+"]",ToString(ib->Momenta()[i]));
  mu2calc.Interprete(mu2tag);
  msg_Debugging()<<"}\n";
}
