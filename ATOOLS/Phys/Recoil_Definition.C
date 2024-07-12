#include "ATOOLS/Phys/Recoil_Definition.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE ATOOLS::Recoil_Definition
#define PARAMETER_TYPE ATOOLS::RecoilDefinition_Key
#define EXACTMATCH true
#include "ATOOLS/Org/Getter_Function.C"

using namespace ATOOLS;

void Recoil_Definition::ShowSyntax(const int mode)
{
  if (!msg_LevelIsInfo() || mode==0) return;
  msg_Out()<<METHOD<<"(): {\n\n";
  RecoilDefinition_Getter::PrintGetterInfo(msg->Out(),25);
  msg_Out()<<"\n}"<<std::endl;
}

namespace ATOOLS {

  class Recoil_InitialState: public Recoil_Definition {
  public:

    Vec4D Recoil(const Cluster_Amplitude *ampl)
    {
      Vec4D rec;
      for (size_t i(0);i<ampl->NIn();++i)
	rec+=ampl->Leg(i)->Mom();
      return rec;
    }

    std::vector<int> RecoilTags(const Cluster_Amplitude *ampl)
    {
      if (!ampl) return std::vector<int>(4);
      std::vector<int> tags(ampl->Legs().size(),0);
      for (size_t i(0);i<ampl->NIn();++i) tags[i]=3;
      return tags;
    }

  };// end of class Recoil_FinalState

  
  class Recoil_FinalState: public Recoil_Definition {
  public:

    Vec4D Recoil(const Cluster_Amplitude *ampl)
    {
      Vec4D rec;
      for (size_t i(ampl->NIn());i<ampl->Legs().size();++i)
	rec+=ampl->Leg(i)->Mom();
      return rec;
    }

    std::vector<int> RecoilTags(const Cluster_Amplitude *ampl)
    {
      if (!ampl) return std::vector<int>(4);
      std::vector<int> tags(ampl->Legs().size(),0);
      for (size_t i(2);i<ampl->Legs().size();++i) tags[i]=3;
      return tags;
    }

  };// end of class Recoil_FinalState

  class Recoil_EWFinalState: public Recoil_Definition {
  public:

    Vec4D Recoil(const Cluster_Amplitude *ampl)
    {
      Vec4D rec;
      for (size_t i(ampl->NIn());i<ampl->Legs().size();++i)
	if (!ampl->Leg(i)->Flav().Strong()) rec+=ampl->Leg(i)->Mom();
      return rec;
    }

    std::vector<int> RecoilTags(const Cluster_Amplitude *ampl)
    {
      std::vector<int> tags(ampl->Legs().size(),0);
      for (size_t i(ampl->NIn());i<ampl->Legs().size();++i)
	if (!ampl->Leg(i)->Flav().Strong()) tags[i]=3;
      return tags;
    }

  };// end of class Recoil_EWFinalState

}

DECLARE_GETTER(Recoil_InitialState,"InitialState",
	       Recoil_Definition,RecoilDefinition_Key);
Recoil_Definition *ATOOLS::Getter
<Recoil_Definition,RecoilDefinition_Key,Recoil_InitialState>::
operator()(const RecoilDefinition_Key &key) const
{ return new Recoil_InitialState(); }

void Getter<Recoil_Definition,RecoilDefinition_Key,Recoil_InitialState>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"complete initial state"; }


DECLARE_GETTER(Recoil_FinalState,"FinalState",
	       Recoil_Definition,RecoilDefinition_Key);
Recoil_Definition *ATOOLS::Getter
<Recoil_Definition,RecoilDefinition_Key,Recoil_FinalState>::
operator()(const RecoilDefinition_Key &key) const
{ return new Recoil_FinalState(); }

void Getter<Recoil_Definition,RecoilDefinition_Key,Recoil_FinalState>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"complete final state"; }

DECLARE_GETTER(Recoil_EWFinalState,"EWFinalState",
	       Recoil_Definition,RecoilDefinition_Key);
Recoil_Definition *ATOOLS::Getter
<Recoil_Definition,RecoilDefinition_Key,Recoil_EWFinalState>::
operator()(const RecoilDefinition_Key &key) const
{ return new Recoil_EWFinalState(); }

void Getter<Recoil_Definition,RecoilDefinition_Key,Recoil_EWFinalState>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"electroweak final state"; }
