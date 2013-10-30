#include "PHASIC++/Process/ME_Generator_Base.H"

#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Phys/Flavour.H"
#include "MODEL/Main/Model_Base.H"

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
  std::vector<size_t> psmassive;
  std::vector<size_t> psmassless;
  dr->VectorFromFile(psmassive,"MASSIVE_PS");
  dr->VectorFromFile(psmassless,"MASSLESS_PS");
  bool massive=!((bool)dr->GetValue<int>("RESPECT_MASSIVE_FLAG",0));
  if (massive) {
    Flavour_Set tempflavs;
    for (size_t i(0);i<psmassless.size();++i) {
      Flavour fl(psmassless[i],0);
      tempflavs.insert(fl);
      tempflavs.insert(fl.Bar());
      msg_Info()<<METHOD<<"(): "<<m_name<<": Using massless PS for "<<fl<<".\n";
    }
    // find m_psmass=allflavs-tempflavs
    for (size_t i(0);i<allflavs.size();++i) {
      if (allflavs[i].IsDummy()) continue;
      bool enter(true);
      for (Flavour_Set::const_iterator it(tempflavs.begin());
           it!=tempflavs.end();++it) {
        if (allflavs[i].Kfcode()==it->Kfcode()) { enter=false;break; }
      }
      if (enter) m_psmass.insert(allflavs[i]);
      if (enter) msg_Tracking()<<METHOD<<"(): "<<m_name
                               <<": Using massive PS for "<<allflavs[i]<<".\n";
    }
  }
  else {
    for (size_t i(0);i<psmassive.size();++i) {
      Flavour fl(psmassive[i],0);
      m_psmass.insert(fl);
      m_psmass.insert(fl.Bar());
      msg_Info()<<METHOD<<"(): "<<m_name<<": Using massive PS for "<<fl<<".\n";
    }
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

void ME_Generator_Base::ShowSyntax(const int mode)
{
  if (!msg_LevelIsInfo() || mode==0) return;
  msg_Out()<<METHOD<<"(): {\n\n";
  ME_Generator_Getter::PrintGetterInfo(msg->Out(),15);
  msg_Out()<<"\n}"<<std::endl;
}
