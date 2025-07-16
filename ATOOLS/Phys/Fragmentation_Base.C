#include "ATOOLS/Phys/Fragmentation_Base.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE ATOOLS::Fragmentation_Base
#define PARAMETER_TYPE ATOOLS::Fragmentation_Getter_Parameters
#define EXACTMATCH false
#include "ATOOLS/Org/Getter_Function.C"

using namespace ATOOLS;
using namespace std;

Fragmentation_Base::Fragmentation_Base()
{
  Settings& s = Settings::GetMainSettings();
  m_shrink = s["COMPRESS_PARTONIC_DECAYS"].SetDefault(true).Get<bool>();
  m_flagpartonics = s["FLAG_PARTONIC_DECAYS"].SetDefault(true).Get<bool>();
}
   
Fragmentation_Base::~Fragmentation_Base() 
{
}
void Fragmentation_Base::Shrink(Blob_List * bloblist) {
  list<Blob *> deleteblobs;
  Particle_Vector * parts;
  for (Blob_List::reverse_iterator blit=bloblist->rbegin();
       blit!=bloblist->rend();++blit) {
    Blob * blob = (*blit);
    if (blob->Type()==btp::Fragmentation) {
      Blob * showerblob(blob->InParticle(0)->ProductionBlob());
      Blob * decblob(showerblob->InParticle(0)->ProductionBlob());
      if (decblob->Type()!=btp::Hadron_Decay) continue;
      showerblob->DeleteInParticles(0);
      showerblob->DeleteOutParticles(0);
      deleteblobs.push_back(blob);
      deleteblobs.push_back(showerblob);
      while (!blob->GetOutParticles().empty()) {
	Particle * part = 
	  blob->RemoveOutParticle(blob->GetOutParticles().front());
	decblob->AddToOutParticles(part);
      }
      decblob->SetStatus(blob_status::needs_hadrondecays);
      decblob->AddData("Partonic",new Blob_Data<int>(m_flagpartonics));
    }
  }
  for (list<Blob *>::iterator blit=deleteblobs.begin();
       blit!=deleteblobs.end();blit++) bloblist->Delete((*blit));
}
void Fragmentation_Base::ReadMassParameters()
{
  auto s = Settings::GetMainSettings()["AHADIC"];
  double mglue =  s["M_GLUE"].SetDefault(0.00).Get<double>();
  double mud =   s["M_UP_DOWN"].SetDefault(0.30).Get<double>();
  double ms =    s["M_STRANGE"].SetDefault(0.40).Get<double>();
  double mc =    s["M_CHARM"].SetDefault(1.80).Get<double>();
  double mb =    s["M_BOTTOM"].SetDefault(5.10).Get<double>();
  double mdiq =  s["M_DIQUARK_OFFSET"].SetDefault(0.30).Get<double>();
  double bind0 = s["M_BIND_0"].SetDefault(0.12).Get<double>();
  double bind1 = s["M_BIND_1"].SetDefault(0.50).Get<double>();
  Flavour(kf_gluon).SetHadMass(mglue);
  Flavour(kf_d).SetHadMass(mud);
  Flavour(kf_u).SetHadMass(mud);
  Flavour(kf_s).SetHadMass(ms);
  Flavour(kf_c).SetHadMass(mc);
  Flavour(kf_b).SetHadMass(mb);
  Flavour(kf_ud_0).SetHadMass((2.*mud+mdiq)*(1.+bind0));
  Flavour(kf_uu_1).SetHadMass((2.*mud+mdiq)*(1.+bind1));
  Flavour(kf_ud_1).SetHadMass((2.*mud+mdiq)*(1.+bind1));
  Flavour(kf_dd_1).SetHadMass((2.*mud+mdiq)*(1.+bind1));
  Flavour(kf_su_0).SetHadMass((ms+mud+mdiq)*(1.+bind0));
  Flavour(kf_sd_0).SetHadMass((ms+mud+mdiq)*(1.+bind0));
  Flavour(kf_su_1).SetHadMass((ms+mud+mdiq)*(1.+bind1));
  Flavour(kf_sd_1).SetHadMass((ms+mud+mdiq)*(1.+bind1));
  Flavour(kf_ss_1).SetHadMass((2.*ms+mdiq)*(1.+bind1));
}


namespace ATOOLS {
  class No_Fragmentation : public Fragmentation_Base {
  public:
    No_Fragmentation(const std::string& shower) {}
    ~No_Fragmentation() {}
    
    Return_Value::code Hadronize(ATOOLS::Blob_List *) {
      return Return_Value::Nothing;
    }
  };
}

DEFINE_FRAGMENTATION_GETTER(No_Fragmentation, "None");
