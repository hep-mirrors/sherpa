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

    Vec4D Recoil(const Cluster_Amplitude *ampl,size_t splits,size_t specs,int)
    {
      Vec4D rec;
      for (size_t i(ampl->NIn());i<ampl->Legs().size();++i)
        if(ampl->Leg(i)->Id()&specs) rec+=ampl->Leg(i)->Mom();
      return rec;
    }

    std::vector<int> RecoilTags
      (const Cluster_Amplitude *ampl,size_t splits,size_t specs,int)
        {
          if (!ampl) return std::vector<int>(4);
          std::vector<int> tags(ampl->Legs().size(),0);
          for (size_t i(ampl->NIn());i<ampl->Legs().size();++i)
            if(ampl->Leg(i)->Id()&specs) tags[i]=3;
          return tags;
        }

    int Mode(const Cluster_Amplitude *ampl,int split) const { return -1; }
  };// end of class Recoil_ColorPartner

  class Recoil_PassiveFinalState: public Recoil_Definition {
  public:


      Vec4D Recoil(const Cluster_Amplitude *ampl,size_t splits,size_t specs,int)
      {
        Vec4D rec;
        for (size_t i(ampl->NIn());i<ampl->Legs().size();++i)
          if((ampl->Leg(i)->Id()&splits)==0) rec+=ampl->Leg(i)->Mom();
        return rec;
      }

      std::vector<int> RecoilTags
        (const Cluster_Amplitude *ampl,size_t splits,size_t specs,int)
          {
            if (!ampl) return std::vector<int>(4);
            std::vector<int> tags(ampl->Legs().size(),0);
            for (size_t i(2);i<ampl->Legs().size();++i)
              if((ampl->Leg(i)->Id()&splits)==0) tags[i]=3;
            return tags;
          }

      int Mode(const Cluster_Amplitude *ampl,int split) const { return -1; }
  };// end of class Recoil_PassiveFinalState

  class Recoil_FinalState: public Recoil_Definition {
  public:

    Recoil_FinalState() {
      m_prefac = -1;
    } 
    
    Vec4D Recoil(const Cluster_Amplitude *ampl,size_t splits,size_t specs,int)
      {
        Vec4D rec;
        for (size_t i(ampl->NIn());i<ampl->Legs().size();++i)
          rec+=ampl->Leg(i)->Mom();
        return rec;
      }

    std::vector<int> RecoilTags
        (const Cluster_Amplitude *ampl,size_t splits,size_t specs,int)
          {
            if (!ampl) return std::vector<int>(4);
            std::vector<int> tags(ampl->Legs().size(),0);
            for (size_t i(2);i<ampl->Legs().size();++i) tags[i]=3;
            return tags;
          }

    int Mode(const Cluster_Amplitude *ampl,int split) const { return -1; }
  };// end of class Recoil_FinalState

  class Recoil_EWFinalState: public Recoil_Definition {
  public:


      Vec4D Recoil(const Cluster_Amplitude *ampl,size_t splits,size_t specs,int)
      {
        Vec4D rec;
        for (size_t i(ampl->NIn());i<ampl->Legs().size();++i)
          if (!ampl->Leg(i)->Flav().Strong()) rec+=ampl->Leg(i)->Mom();
        return rec;
      }
        
      std::vector<int> RecoilTags
        (const Cluster_Amplitude *ampl,size_t splits,size_t specs,int)
          {
            std::vector<int> tags(ampl->Legs().size(),0);
            for (size_t i(ampl->NIn());i<ampl->Legs().size();++i)
              if (!ampl->Leg(i)->Flav().Strong()) tags[i]=3;
            return tags;
          }

    int Mode(const Cluster_Amplitude *ampl,int split) const { return -1; }
  };// end of class Recoil_EWFinalState

  class Recoil_EWChargedMode: public Recoil_Definition {
    enum {
      all = -1,
      pos = 0,
      neg = 1
    };

    Vec4D Recoil(const Cluster_Amplitude *ampl,size_t,size_t,int recmode)
    {
      Vec4D rec;
      for (size_t i(ampl->NIn());i<ampl->Legs().size();++i)
	if (!ampl->Leg(i)->Flav().Strong() &&
            (recmode==pos && ampl->Leg(i)->Flav().Charge() > 0 ||
             recmode==neg && ampl->Leg(i)->Flav().Charge() < 0 ||
             recmode==all))
          rec+=ampl->Leg(i)->Mom();
      return rec;
    }

    std::vector<int> RecoilTags(const Cluster_Amplitude *ampl,size_t,size_t,int recmode)
      {
        std::vector<int> tags(ampl->Legs().size(),0);
        for (size_t i(ampl->NIn());i<ampl->Legs().size();++i)
          if (!ampl->Leg(i)->Flav().Strong() &&
              (recmode==pos && ampl->Leg(i)->Flav().Charge() > 0 ||
               recmode==neg && ampl->Leg(i)->Flav().Charge() < 0 ||
               recmode==all)) {
            tags[i]=3;
          }
        return tags;
      }

    int Mode(const Cluster_Amplitude *ampl,int splits) const
    {
      double charge = 0;
      for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
        if((ampl->Leg(i)->Id()&splits)==0) charge += ampl->Leg(i)->Flav().Charge();
      }
      if(charge == 0) {
        //THROW(fatal_error,"No unique charge assignment. Recoil definition invalid.");
        return all;
      }
      return charge>0?neg:pos;
    }
    
  };// end of class Recoil_EWChargeMode


  class Recoil_StrongFinalState: public Recoil_Definition {
  public:

    Vec4D Recoil(const Cluster_Amplitude *ampl,size_t splits,size_t specs,int)
    {
      Vec4D rec;
      for (size_t i(ampl->NIn());i<ampl->Legs().size();++i)
	if (ampl->Leg(i)->Flav().Strong()) rec+=ampl->Leg(i)->Mom();
      return rec;
    }

    std::vector<int> RecoilTags
    (const Cluster_Amplitude *ampl,size_t splits,size_t specs,int)
    {
      std::vector<int> tags(ampl->Legs().size(),0);
      for (size_t i(ampl->NIn());i<ampl->Legs().size();++i)
	if (ampl->Leg(i)->Flav().Strong()) tags[i]=3;
      return tags;
    }

    int Mode(const Cluster_Amplitude *ampl,int split) const { return -1; }
  };// end of class Recoil_StrongFinalState

  class Recoil_PassiveStrongFinalState: public Recoil_Definition {
  public:

    Vec4D Recoil(const Cluster_Amplitude *ampl,size_t splits,size_t specs,int)
    {
      Vec4D rec;
      for (size_t i(ampl->NIn());i<ampl->Legs().size();++i)
        if((ampl->Leg(i)->Id()&splits)==0 && ampl->Leg(i)->Flav().Strong()) rec+=ampl->Leg(i)->Mom();
      return rec;
    }

    std::vector<int> RecoilTags
    (const Cluster_Amplitude *ampl,size_t splits,size_t specs,int)
    {
      if (!ampl) return std::vector<int>(4);
      std::vector<int> tags(ampl->Legs().size(),0);
      for (size_t i(2);i<ampl->Legs().size();++i)
        if((ampl->Leg(i)->Id()&splits)==0 && ampl->Leg(i)->Flav().Strong()) tags[i]=3;
      return tags;
    }

    int Mode(const Cluster_Amplitude *ampl,int split) const { return -1; }
  };// end of class Recoil_PassiveFinalState

  
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

DECLARE_GETTER(Recoil_EWChargedMode,"EWChargedMode",
	       Recoil_Definition,RecoilDefinition_Key);
Recoil_Definition *ATOOLS::Getter
<Recoil_Definition,RecoilDefinition_Key,Recoil_EWChargedMode>::
operator()(const RecoilDefinition_Key &key) const
{ return new Recoil_EWChargedMode(); }

void Getter<Recoil_Definition,RecoilDefinition_Key,Recoil_EWChargedMode>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"electroweak final state"; }
DECLARE_GETTER(Recoil_StrongFinalState,"StrongFinalState",
	       Recoil_Definition,RecoilDefinition_Key);
Recoil_Definition *ATOOLS::Getter
<Recoil_Definition,RecoilDefinition_Key,Recoil_StrongFinalState>::
operator()(const RecoilDefinition_Key &key) const
{ return new Recoil_StrongFinalState(); }

void Getter<Recoil_Definition,RecoilDefinition_Key,Recoil_StrongFinalState>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"strongly interacting final state"; }

DECLARE_GETTER(Recoil_PassiveStrongFinalState,"PassiveStrongFinalState",
	       Recoil_Definition,RecoilDefinition_Key);
Recoil_Definition *ATOOLS::Getter
<Recoil_Definition,RecoilDefinition_Key,Recoil_PassiveStrongFinalState>::
operator()(const RecoilDefinition_Key &key) const
{ return new Recoil_PassiveStrongFinalState(); }

void Getter<Recoil_Definition,RecoilDefinition_Key,Recoil_PassiveStrongFinalState>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"strongly interacting final state except splitter (and emission)"; }

