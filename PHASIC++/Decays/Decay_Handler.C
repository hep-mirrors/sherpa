#include "PHASIC++/Decays/Decay_Handler.H"

#include "PHASIC++/Decays/Decay_Channel.H"
#include "PHASIC++/Decays/Decay_Table.H"
#include "PHASIC++/Decays/Decay_Map.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

namespace METOOLS {
  class Spin_Density;
  class Amplitude2_Tensor;
}

/// TODO: move ms to Decay_Handler instead of Decay_Map?
Decay_Handler::Decay_Handler()
{
}

Decay_Handler::~Decay_Handler()
{
  if (p_decaymap) delete p_decaymap; p_decaymap=NULL;
}

METOOLS::Amplitude2_Tensor* Decay_Handler::FillOnshellDecay(Blob *blob,
                                                           METOOLS::Spin_Density* sigma)
{
  DEBUG_FUNC("");
  Decay_Channel* dc(NULL);
  Blob_Data_Base* data = (*blob)["dc"];
  if (data) {
    dc=data->Get<Decay_Channel*>();
  }
  else {
    Decay_Table* table=p_decaymap->FindDecay(blob->InParticle(0)->Flav());
    if (table==NULL) {
      msg_Error()<<METHOD<<"Error: Did not find "<<blob->InParticle(0)->Flav()
                 <<" in decay tables."<<endl;
      throw Return_Value::Retry_Event;
    }
    dc=table->Select();
    blob->AddData("dc",new Blob_Data<Decay_Channel*>(dc));
  }
  if (!dc) THROW(fatal_error,"No decay channel found for "
                             +blob->InParticle(0)->Flav().IDName()+".");
  msg_Debugging()<<*dc<<std::endl;

  Particle* inpart=blob->InParticle(0);
  inpart->SetStatus(part_status::decayed);
  Flavour flav; Particle* particle=NULL;
  for (size_t i=1; i<dc->Flavs().size(); ++i) {
    flav=dc->Flavs()[i];
    if (inpart->Flav().IsAnti()!=dc->GetDecaying().IsAnti()) flav = flav.Bar();
    particle = new Particle(0, flav);
    particle->SetFinalMass(Mass(flav));
    particle->SetStatus(part_status::active);
    particle->SetNumber();
    particle->SetInfo('D');
    blob->AddToOutParticles( particle );
  }
  
  size_t n=1+blob->NOutP();
  vector<Vec4D> moms(n);
  moms[0]=inpart->Momentum();

  METOOLS::Amplitude2_Tensor* amps(NULL);
  if (sigma) {
    std::vector<Particle*> parts;
    parts.insert(parts.end(), blob->InParticle(0));
    Particle_Vector outparts=blob->GetOutParticles();
    parts.insert(parts.end(), outparts.begin(), outparts.end());
    dc->GenerateKinematics(moms,
                           inpart->Flav().IsAnti()!=dc->GetDecaying().IsAnti(),
                           sigma,parts);
    amps=dc->Amps();
  }
  else {
    dc->GenerateKinematics(moms,
                           inpart->Flav().IsAnti()!=dc->GetDecaying().IsAnti());
  }
  for (size_t i=1; i<n; i++) blob->GetOutParticles()[i-1]->SetMomentum(moms[i]);
  msg_Debugging()<<*blob<<std::endl;
  return amps;
}

double Decay_Handler::DecayWidth(const ATOOLS::Flavour& flav)
{
  if (p_decaymap) {
    Decay_Table* table=p_decaymap->FindDecay(flav);
    if (table) return table->TotalWidth();
    else THROW(fatal_error, "Internal Error 1");
  }
  else return flav.Width();
}

bool Decay_Handler::DiceMass(ATOOLS::Particle* p, double max)
{
  DEBUG_FUNC(p->Flav()<<"  m<"<<max);
  Blob* decayblob=p->DecayBlob();
  Blob_Data_Base* data = (*decayblob)["dc"];
  if (data) {
    Decay_Channel* dc = data->Get<Decay_Channel*>();
    if (!dc) THROW(fatal_error,"Missing decay channel for "
                               +decayblob->ShortProcessName()+".");
    if (dc->NOut()==1 && p->FinalMass()<max) return true;
    double width = p_decaymap->FindDecay(p->Flav())->TotalWidth();
    double mass=dc->GenerateMass(max, width);
    if (mass>0.0) p->SetFinalMass(mass);
    else return false;
  }
  return true;
}

void Decay_Handler::CleanUp()
{
  if (p_decaymap) p_decaymap->ResetCounters();
}


#include "EXTRA_XS/One2Three/Comix1to3.H"
ATOOLS::Flavour Decay_Handler::PropFlav(ATOOLS::Blob* blob)
{
  // TODO: replace hack
  Decay_Channel* dc(NULL);
  Blob_Data_Base* data = (*blob)["dc"];
  if (!data) return Flavour(kf_none);

  dc=data->Get<Decay_Channel*>();
  DEBUG_VAR(*dc);
  EXTRAXS::Comix1to3* amp=dynamic_cast<EXTRAXS::Comix1to3*>(dc->GetDiagrams()[0]);
  if (!amp) THROW(fatal_error, "Internal error.");
  return amp->Prop();
}







