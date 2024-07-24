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

  class Recoil_ColorPartner: public Recoil_Definition {
  public:

    Vec4D Recoil(const Cluster_Amplitude *ampl,int,int,int cspect)
    {
      return ampl->Legs()[cspect]->Mom();
    }

    std::vector<int> RecoilTags(const Cluster_Amplitude *ampl,int,int,int cspect)
    {
      if (!ampl) return std::vector<int>(4);
      std::vector<int> tags(ampl->Legs().size(),0);
      tags[cspect] = 3;
      return tags;
    }

  };// end of class Recoil_PassiveFinalState

  class Recoil_PassiveFinalState: public Recoil_Definition {
  public:

    Vec4D Recoil(const Cluster_Amplitude *ampl,int split,int em,int)
    {
      Vec4D rec;
      for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
        if(i==split || i==em) continue;
	rec+=ampl->Leg(i)->Mom();
      }
      return rec;
    }

    std::vector<int> RecoilTags(const Cluster_Amplitude *ampl,int split,int em,int)
    {
      if (!ampl) return std::vector<int>(4);
      std::vector<int> tags(ampl->Legs().size(),0);
      for (size_t i(2);i<ampl->Legs().size();++i) {
        if(i==split || i==em) continue;
        tags[i]=3;
      }
      return tags;
    }

  };// end of class Recoil_PassiveFinalState


  class Recoil_FinalState: public Recoil_Definition {
  public:

    Vec4D Recoil(const Cluster_Amplitude *ampl,int,int,int)
    {
      Vec4D rec;
      for (size_t i(ampl->NIn());i<ampl->Legs().size();++i)
	rec+=ampl->Leg(i)->Mom();
      return rec;
    }

    std::vector<int> RecoilTags(const Cluster_Amplitude *ampl,int,int,int)
    {
      if (!ampl) return std::vector<int>(4);
      std::vector<int> tags(ampl->Legs().size(),0);
      for (size_t i(2);i<ampl->Legs().size();++i) tags[i]=3;
      return tags;
    }

  };// end of class Recoil_FinalState

  class Recoil_EWFinalState: public Recoil_Definition {
  public:

    Vec4D Recoil(const Cluster_Amplitude *ampl,int,int,int)
    {
      Vec4D rec;
      for (size_t i(ampl->NIn());i<ampl->Legs().size();++i)
	if (!ampl->Leg(i)->Flav().Strong()) rec+=ampl->Leg(i)->Mom();
      return rec;
    }

    std::vector<int> RecoilTags(const Cluster_Amplitude *ampl,int,int,int)
    {
      std::vector<int> tags(ampl->Legs().size(),0);
      for (size_t i(ampl->NIn());i<ampl->Legs().size();++i)
	if (!ampl->Leg(i)->Flav().Strong()) tags[i]=3;
      return tags;
    }

  };// end of class Recoil_EWFinalState

}

DECLARE_GETTER(Recoil_ColorPartner,"ColorPartner",
	       Recoil_Definition,RecoilDefinition_Key);
Recoil_Definition *ATOOLS::Getter
<Recoil_Definition,RecoilDefinition_Key,Recoil_ColorPartner>::
operator()(const RecoilDefinition_Key &key) const
{ return new Recoil_ColorPartner(); }

void Getter<Recoil_Definition,RecoilDefinition_Key,Recoil_ColorPartner>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"color partner"; }


DECLARE_GETTER(Recoil_PassiveFinalState,"PassiveFinalState",
	       Recoil_Definition,RecoilDefinition_Key);
Recoil_Definition *ATOOLS::Getter
<Recoil_Definition,RecoilDefinition_Key,Recoil_PassiveFinalState>::
operator()(const RecoilDefinition_Key &key) const
{ return new Recoil_PassiveFinalState(); }

void Getter<Recoil_Definition,RecoilDefinition_Key,Recoil_PassiveFinalState>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"final state except splitter (and emission)"; }


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
