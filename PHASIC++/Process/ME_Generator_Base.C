#include "PHASIC++/Process/ME_Generator_Base.H"

#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Phys/Flavour.H"
#include "MODEL/Main/Model_Base.H"
#include <algorithm>

#define COMPILE__Getter_Function
#define OBJECT_TYPE PHASIC::ME_Generator_Base
#define PARAMETER_TYPE PHASIC::ME_Generator_Key
#define EXACTMATCH false
#include "ATOOLS/Org/Getter_Function.C"

using namespace PHASIC;
using namespace ATOOLS;

ME_Generator_Base::~ME_Generator_Base()
{
}

void ME_Generator_Base::SetPSMasses(Data_Reader *const dr)
{
  ATOOLS::Flavour_Vector allflavs(MODEL::s_model->IncludedFlavours());
  std::vector<size_t> psmassive,psmassless;
  std::vector<size_t> defpsmassive,defpsmassless;
  dr->VectorFromFile(psmassive,"MASSIVE_PS");
  dr->VectorFromFile(psmassless,"MASSLESS_PS");
  bool respect=(bool)dr->GetValue<int>("RESPECT_MASSIVE_FLAG",0);
  // check consistency
  for (size_t i(0);i<psmassive.size();++i)
    if (std::find(psmassless.begin(),psmassless.end(),psmassive[i])!=
        psmassless.end()) THROW(fatal_error,"Inconsinstent input.");
  for (size_t i(0);i<psmassless.size();++i)
    if (Flavour(psmassless[i]).IsMassive())
      THROW(fatal_error,"Cannot shower massive particle massless.");
  // set defaults
  // respect=0 -> def: dusgy massless, rest massive
  // respect=1 -> def: only massive massive, rest massless
  // TODO: need to fill in those that are massive already?
  if (!respect) {
    defpsmassless.push_back(kf_d);
    defpsmassless.push_back(kf_u);
    defpsmassless.push_back(kf_s);
    defpsmassless.push_back(kf_gluon);
    defpsmassless.push_back(kf_photon);
    for (size_t i(0);i<allflavs.size();++i) {
      if (allflavs[i].IsDummy()) continue;
      size_t kf(allflavs[i].Kfcode());
      bool add(true);
      for (size_t j(0);j<defpsmassive.size();++j)
        if (kf==defpsmassive[j]) { add=false; break; }
      for (size_t j(0);j<defpsmassless.size();++j)
        if (kf==defpsmassless[j]) { add=false; break; }
      if (add)  defpsmassive.push_back(kf);
    }
  }
  else {
    for (size_t i(0);i<allflavs.size();++i) {
      if (allflavs[i].IsDummy()) continue;
      size_t kf(allflavs[i].Kfcode());
      bool add(true);
      for (size_t j(0);j<defpsmassive.size();++j)
        if (kf==defpsmassive[j]) { add=false; break; }
      for (size_t j(0);j<defpsmassless.size();++j)
        if (kf==defpsmassless[j]) { add=false; break; }
      if (add && allflavs[i].IsMassive())  defpsmassive.push_back(kf);
      if (add && !allflavs[i].IsMassive()) defpsmassless.push_back(kf);
    }
  }
  // then remove and add those specified manually
  for (size_t i(0);i<psmassive.size();++i) {
    defpsmassless.erase(std::remove(defpsmassless.begin(),defpsmassless.end(),
                                    psmassive[i]),defpsmassless.end());
    if (std::find(defpsmassive.begin(),defpsmassive.end(),psmassive[i])==
        defpsmassive.end()) defpsmassive.push_back(psmassive[i]);
  }
  for (size_t i(0);i<psmassless.size();++i) {
    defpsmassive.erase(std::remove(defpsmassive.begin(),defpsmassive.end(),
                                   psmassless[i]),defpsmassive.end());
    if (std::find(defpsmassless.begin(),defpsmassless.end(),psmassless[i])==
        defpsmassless.end()) defpsmassless.push_back(psmassless[i]);
  }
  // fill massive ones into m_psmass
  for (size_t i(0);i<defpsmassive.size();++i) {
    Flavour fl(defpsmassive[i],0);
    m_psmass.insert(fl);
    m_psmass.insert(fl.Bar());
    msg_Tracking()<<METHOD<<"(): "<<m_name<<": Using massive PS for "<<fl<<".\n";
  }
  Flavour_Vector mf(m_psmass.begin(),m_psmass.end());
  msg_Info()<<METHOD<<"(): Massive PS flavours for "<<m_name<<": "
                    <<mf<<std::endl;
}

namespace PHASIC {

  class ShiftMasses_Energy: public Function_Base {
  private:
    std::vector<double> m_m2, m_p2;
  public:
    ShiftMasses_Energy(Mass_Selector *const ms,
		    Cluster_Amplitude *const ampl,int mode)
    {
      if (mode<0) {
	for (size_t i(0);i<ampl->NIn();++i) {
	  m_p2.push_back(ampl->Leg(i)->Mom().PSpat2());
	  m_m2.push_back(ms->Mass2(ampl->Leg(i)->Flav()));
	}
      }
      else {
	for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
	  m_p2.push_back(ampl->Leg(i)->Mom().PSpat2());
	  m_m2.push_back(ms->Mass2(ampl->Leg(i)->Flav()));
	}
      }
    }
    virtual double operator()(double x)
    {
      double E=0.0;
      for (size_t i(0);i<m_m2.size();++i)
	E+=sqrt(m_m2[i]+x*x*m_p2[i]);
      return E;
    }
  };// end of class ShiftMasses_Energy

}// end of namespace PHASIC

int ME_Generator_Base::ShiftMasses(Cluster_Amplitude *const ampl)
{
  if (m_psmass.empty()) return 0;
  bool run=false;
  Vec4D cms, pb[2]={ampl->Leg(0)->Mom(),ampl->Leg(1)->Mom()};
  for (size_t i(0);i<ampl->Legs().size();++i) {
    if (i<ampl->NIn()) cms-=ampl->Leg(i)->Mom();
    if (m_psmass.find(ampl->Leg(i)->Flav())!=
	m_psmass.end()) run=true;
  }
  if (!run) return 1;
  DEBUG_FUNC(m_name);
  msg_Debugging()<<"Before shift: "<<*ampl<<"\n";
  Poincare boost(cms);
  boost.Boost(cms);
  for (size_t i(0);i<ampl->Legs().size();++i)
    ampl->Leg(i)->SetMom(boost*ampl->Leg(i)->Mom());
  boost.Invert();
  if (ampl->NIn()>1) {
    ShiftMasses_Energy etot(this,ampl,-1);
    double xi(etot.WDBSolve(cms[0],0.0,1.0));
    if (!IsEqual(etot(xi),cms[0],rpa->gen.Accu())) return -1;
    for (size_t i(0);i<ampl->NIn();++i) {
      Vec4D p(xi*ampl->Leg(i)->Mom());
      p[0]=-sqrt(Mass2(ampl->Leg(i)->Flav())+p.PSpat2());
      ampl->Leg(i)->SetMom(boost*p);
    }
  }
  ShiftMasses_Energy etot(this,ampl,1);
  double xi(etot.WDBSolve(cms[0],0.0,1.0));
  if (!IsEqual(etot(xi),cms[0],rpa->gen.Accu())) return -1;
  for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
    Vec4D p(xi*ampl->Leg(i)->Mom());
    p[0]=sqrt(Mass2(ampl->Leg(i)->Flav())+p.PSpat2());
    ampl->Leg(i)->SetMom(boost*p);
  }
  msg_Debugging()<<"After shift: "<<*ampl<<"\n";
  if (ampl->NIn()==2) {
    if (ampl->Leg(0)->Mom()[3]*pb[0][3]<0.0) return -1;
    if (ampl->Leg(1)->Mom()[3]*pb[1][3]<0.0) return -1;
  }
  return 1;
}

double ME_Generator_Base::Mass(const ATOOLS::Flavour &fl) const
{
  if (m_massmode==0) return fl.Mass();
  if (m_psmass.find(fl)!=m_psmass.end()) return fl.Mass(true);
  return fl.Mass();
}

void ME_Generator_Base::SetClusterDefinitions
(PDF::Cluster_Definitions_Base *const defs)
{
}

void ME_Generator_Base::PreCluster
(Process_Base *const proc,const ATOOLS::Vec4D_Vector &p)
{
}

Cluster_Amplitude *ME_Generator_Base::ClusterConfiguration
(Process_Base *const proc,const Vec4D_Vector &p,const size_t &mode)
{
  return NULL;
}

void ME_Generator_Base::ShowSyntax(const int mode)
{
  if (!msg_LevelIsInfo() || mode==0) return;
  msg_Out()<<METHOD<<"(): {\n\n";
  ME_Generator_Getter::PrintGetterInfo(msg->Out(),15);
  msg_Out()<<"\n}"<<std::endl;
}
