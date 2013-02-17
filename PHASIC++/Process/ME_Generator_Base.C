#include "PHASIC++/Process/ME_Generator_Base.H"

#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/Poincare.H"

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
  std::vector<double> psmass;
  dr->VectorFromFile(psmass,"MASSIVE_PS");
  for (size_t i(0);i<psmass.size();++i) {
    Flavour fl((int)psmass[i],0);
    m_psmass.insert(fl);
    m_psmass.insert(fl.Bar());
    msg_Info()<<m_name<<": Using massive PS for "<<fl<<".\n";
  }
}

bool ME_Generator_Base::ShiftMasses(Cluster_Amplitude *const ampl)
{
  if (m_psmass.empty()) return true;
  DEBUG_FUNC(m_name);
  msg_Debugging()<<"Before shift: "<<*ampl<<"\n";
  for (size_t i(0);i<ampl->Legs().size();++i) {
    Cluster_Leg *li(ampl->Leg(i));
    if (m_psmass.find(li->Flav())==m_psmass.end()) continue;
    Vec4D pk, pi(li->Mom());
    if (i<ampl->NIn()) pk=ampl->Leg(1-i)->Mom();
    else {
      for (size_t j(ampl->NIn());j<ampl->Legs().size();++j)
	if (i!=j) pk+=ampl->Leg(j)->Mom();
    }
    Vec4D Q(pk+pi);
    double sk(pk.Abs2()), si(li->Mom().Abs2());
    double mi2(Mass2(li->Flav())), Q2(Q.Abs2());
    double po(sqr(Q2-si-sk)-4.0*si*sk);
    double pn(sqr(Q2-mi2-sk)-4.0*mi2*sk);
    if (pn<0.0 ^ po<0.0) return false;
    Vec4D npk(sqrt(pn/po)*(pk-(Q*pk)/Q2*Q)+(Q2+sk-mi2)/(2.0*Q2)*Q);
    li->SetMom(Q-npk);
    if (i<ampl->NIn()) ampl->Leg(1-i)->SetMom(npk);
    else {
      if (ampl->Legs().size()==ampl->NIn()+2) {
	for (size_t j(ampl->NIn());j<ampl->Legs().size();++j)
	  if (i!=j) ampl->Leg(j)->SetMom(npk);
      }
      else {
	Poincare ocms(pk), ncms(npk);
	ncms.Invert();
	for (size_t j(ampl->NIn());j<ampl->Legs().size();++j)
	  if (i!=j) ampl->Leg(j)->SetMom
		      (ncms*(ocms*ampl->Leg(j)->Mom()));
      }
    }
    msg_Debugging()<<"After shifting "<<i<<": "<<*ampl<<"\n";
  }
  return true;
}

double ME_Generator_Base::Mass(const ATOOLS::Flavour &fl) const
{
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
