#include "HADRONS++/Main/Hadron_Decay_Handler.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Read_Write_Base.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "HADRONS++/Main/Hadron_Decay_Map.H"
#include "HADRONS++/Main/Hadron_Decay_Table.H"
#include "HADRONS++/Main/Hadron_Decay_Channel.H"
#include "HADRONS++/Main/Mixing_Handler.H"
#include "METOOLS/SpinCorrelations/Spin_Density.H"
#include "METOOLS/SpinCorrelations/Decay_Matrix.H"
#include <algorithm>

using namespace ATOOLS;
using namespace PHASIC;
using namespace std;
using namespace HADRONS;
using namespace METOOLS;

Hadron_Decay_Handler::Hadron_Decay_Handler() :
  Decay_Handler()
{
  Settings& s = Settings::GetMainSettings();
  RegisterSettings();
  const double maxproperlifetime{ s["MAX_PROPER_LIFETIME"].Get<double>() };
  if(maxproperlifetime > 0.0) {
    for(KFCode_ParticleInfo_Map::const_iterator kfit(s_kftable.begin());
	kfit!=s_kftable.end();++kfit) {
      Flavour flav(kfit->first);
      if (flav.IsOn() && flav.IsHadron() && !flav.IsStable() &&
          0.197e-12>maxproperlifetime*flav.Width() && flav.Kfcode()!=kf_K)
      {
        flav.SetStable(true);
      }
    }
  }

  string decaydatazip = s["DECAYDATA"].Get<std::string>();
  string decaydata = decaydatazip.replace(decaydatazip.length()-4,4,"/");
  string decayfile = s["DECAYFILE"].Get<std::string>();
  string decayconstfile = s["DECAYCONSTFILE"].Get<std::string>();
  string bdecayfile = s["B_DECAYFILE"].Get<std::string>();
  string cdecayfile = s["C_DECAYFILE"].Get<std::string>();
  string aliasfile = s["HADRONALIASESFILE"].Get<std::string>();
  string aliasdecayfile = s["ALIASDECAYFILE"].Get<std::string>();
  My_In_File::OpenDB(decaydatazip);
  Hadron_Decay_Map * dmap = new Hadron_Decay_Map(this);
  dmap->ReadInConstants(decaydata, decayconstfile);
  dmap->ReadInPartonicDecays(Flavour(kf_b),decaydata,bdecayfile);
  dmap->ReadInPartonicDecays(Flavour(kf_c),decaydata,cdecayfile);
  dmap->ReadHadronAliases(decaydata, aliasfile);
  dmap->Read(decaydata, decayfile, true);
  dmap->Read(decaydata, aliasdecayfile);
  dmap->Initialise();
  dmap->ReadFixedTables(decaydata, "FixedDecays.dat");
  p_decaymap=dmap;

  int decay_tau_hard=s["DECAY_TAU_HARD"].Get<int>();
  for (Decay_Map::iterator it=p_decaymap->begin(); it!=p_decaymap->end(); ++it) {
    Flavour flav(it->first);
    if (flav.Kfcode()==kf_tau && decay_tau_hard) continue;
    flav.SetDecayHandler(this);
  }
  
  p_mixinghandler = new Mixing_Handler();
  p_mixinghandler->SetModel(dmap->StartModel());
  dmap->SetMixingHandler(p_mixinghandler);
  My_In_File::CloseDB(decaydatazip);
}

Hadron_Decay_Handler::~Hadron_Decay_Handler()
{
  Hadron_Decay_Map* dmap=dynamic_cast<Hadron_Decay_Map*>(p_decaymap);
  delete dmap; p_decaymap=NULL;
  delete p_mixinghandler; p_mixinghandler=NULL;
}

void Hadron_Decay_Handler::RegisterSettings()
{
  Settings& s = Settings::GetMainSettings();
  s["HADRON_DECAYS_QED_CORRECTIONS"].SetDefault(1);
  s["MAX_PROPER_LIFETIME"].SetDefault(-1.0);
  const auto path = s["DECAYPATH"]
    .SetDefault(rpa->gen.Variable("SHERPA_SHARE_PATH") + "/")
    .Get<std::string>();
  s["DECAYDATA"].SetDefault(path + "Decaydata.zip");
  s["DECAYFILE"].SetDefault("HadronDecays.dat");
  s["DECAYCONSTFILE"].SetDefault("HadronConstants.dat");
  s["B_DECAYFILE"].SetDefault("Partonic_b/Decays.dat");
  s["C_DECAYFILE"].SetDefault("Partonic_c/Decays.dat");
  s["HADRONALIASESFILE"].SetDefault("HadronAliases.dat");
  s["ALIASDECAYFILE"].SetDefault("AliasDecays.dat");
  s["SOFT_MASS_SMEARING"].SetDefault(1);
  s["DECAY_TAU_HARD"].SetDefault(0);
}

bool Hadron_Decay_Handler::VetoDecayAndPrepForNew(ATOOLS::Blob* blob)
{
  if (blob->Type()==btp::Hadron_Decay && (*blob)["Partonic"]!=NULL) {
    if (RejectExclusiveChannelsFromFragmentation(blob)) {
      return true;
    }
  }
  return false;
}

void Hadron_Decay_Handler::AfterTreatInitialBlob(Blob* initialblob, Blob_List* bloblist)
{
  MakeShowerBlobsForPartonicDecays(initialblob, bloblist);
}


void Hadron_Decay_Handler::BeforeFillDecayTree(Blob * blob)
{
  p_mixinghandler->PerformMixing(blob->InParticle(0));
}

Amplitude2_Tensor*
Hadron_Decay_Handler::FillOnshellDecay(Blob *blob, Spin_Density* sigma)
{
  Amplitude2_Tensor* amps=Decay_Handler::FillOnshellDecay(blob,sigma);
  Decay_Channel* dc=(*blob)["dc"]->Get<Decay_Channel*>();
  Hadron_Decay_Channel* hdc=dynamic_cast<Hadron_Decay_Channel*>(dc);
  if(hdc && !hdc->SetColorFlow(blob)) {
    msg_Error()<<METHOD<<" failed to set color flow, retrying event."<<endl
               <<*blob<<endl;
    throw Return_Value::Retry_Event;
  }
  return amps;
}

void Hadron_Decay_Handler::CreateDecayBlob(Blob_List* bloblist, Particle* inpart)
{
  DEBUG_FUNC(inpart->RefFlav());
  if(inpart->DecayBlob()) THROW(fatal_error,"Decay blob already exists.");
  if(inpart->Time()==0.0) inpart->SetTime();
  Blob* blob = bloblist->AddBlob(btp::Hadron_Decay);
  blob->SetStatus(blob_status::needs_extraQED);
  blob->AddToInParticles(inpart);
  SetPosition(blob);
  blob->SetTypeSpec("Sherpa");
  Decay_Table* table=p_decaymap->FindDecay(blob->InParticle(0)->Flav());
  if (table==NULL) {
    msg_Error()<<METHOD<<" decay table not found, retrying event."<<endl
               <<*blob<<endl;
    throw Return_Value::Retry_Event;
  }
  blob->AddData("dc",new Blob_Data<Decay_Channel*>(table->Select()));
  blob->AddData("p_onshell",new Blob_Data<Vec4D>(inpart->Momentum()));
  DEBUG_VAR(inpart->Momentum());
}

void Hadron_Decay_Handler::SetPosition(ATOOLS::Blob* blob)
{
  Particle* inpart = blob->InParticle(0);
  if(inpart->Flav().Kfcode()==kf_K) {
    blob->SetPosition(inpart->XProd());
    return;
  }
  
  // boost lifetime into lab
  double gamma = 1./rpa->gen.Accu();
  if (inpart->Flav().HadMass()>rpa->gen.Accu()) {
    gamma = inpart->E()/inpart->Flav().HadMass();
  }
  else {
    double q2    = dabs(inpart->Momentum().Abs2());
    if (q2>rpa->gen.Accu()) gamma = inpart->E()/sqrt(q2);
  }
  double lifetime_boosted = gamma * inpart->Time();
  
  Vec3D      spatial = inpart->Distance( lifetime_boosted );
  Vec4D     position = Vec4D( lifetime_boosted*rpa->c(), spatial );
  blob->SetPosition( inpart->XProd() + position ); // in mm
}

bool Hadron_Decay_Handler::RejectExclusiveChannelsFromFragmentation(Blob* decblob)
{
  if (decblob->Type()!=btp::Hadron_Decay || (*decblob)["Partonic"]==NULL)
    return false;

  DEBUG_FUNC(decblob->InParticle(0)->Flav());
  Return_Value::IncCall(METHOD);
        
  bool anti=false;
  Decay_Map::iterator dt=
    p_decaymap->find(decblob->InParticle(0)->Flav());
  if (dt==p_decaymap->end()) {
    anti=true;
    dt=p_decaymap->find(decblob->InParticle(0)->Flav().Bar());
    if (dt==p_decaymap->end() || dt->second.size()!=1) {
      PRINT_INFO("Internal error.");
      throw Return_Value::Retry_Event;
    }
  }
  Flavour_Vector tmp, tmpno;
  for(int i=0;i<decblob->NOutP();i++) {
    Flavour flav(anti?decblob->OutParticle(i)->Flav().Bar():
                 decblob->OutParticle(i)->Flav());
    tmp.push_back(flav);
    if (!flav.IsPhoton()) tmpno.push_back(flav);
  }
  
  std::sort(tmp.begin(), tmp.end(), Decay_Channel::FlavourSort);
  Flavour_Vector compflavs(tmp.size()+1), compflavsno(tmpno.size()+1);
  compflavs[0]=(anti?decblob->InParticle(0)->Flav().Bar():decblob->InParticle(0)->Flav());
  for (size_t i=0; i<tmp.size(); ++i) compflavs[i+1]=tmp[i];
  if (tmpno.size()!=tmp.size()) {
    std::sort(tmpno.begin(), tmpno.end(), Decay_Channel::FlavourSort);
    compflavsno[0]=(anti?decblob->InParticle(0)->Flav().Bar():
                    decblob->InParticle(0)->Flav());
    for (size_t i=0; i<tmpno.size(); ++i) compflavsno[i+1]=tmpno[i];
  }
  
  if (dt->second[0]->GetDecayChannel(compflavs) ||
      (tmpno.size()!=tmp.size() && 
       dt->second[0]->GetDecayChannel(compflavsno))) {
    Return_Value::IncRetryPhase(METHOD);
    DEBUG_INFO("rejected. retrying decay.");
    decblob->DeleteOutParticles();
    decblob->InParticle(0)->SetStatus(part_status::active);
    decblob->AddStatus(blob_status::internal_flag);
    Flavour infl(decblob->InParticle(0)->Flav());
    Hadron_Decay_Table * hdt = 
      dynamic_cast<Hadron_Decay_Table*>(p_decaymap->FindDecay(infl));
    hdt->Select(decblob);
    return true;
  }
  else {
    DEBUG_INFO("not found as exclusive. accepted.");
    return false;
  }
}

void Hadron_Decay_Handler::MakeShowerBlobsForPartonicDecays(Blob* initialblob, Blob_List* bloblist)
{
  DEBUG_FUNC(bloblist->size());
  for (auto initialblob: *bloblist) {
    if (initialblob->Type()!=btp::Hadron_Decay) continue;
    bool partonic=false;
    for (auto part: initialblob->GetOutParticles())
      if (part->Flav().IsQCD() && part->DecayBlob()==NULL) partonic=true;
    if (!partonic) continue;
    DEBUG_INFO("found partonic decay: "<<*initialblob);
    Blob* showerblob = bloblist->AddBlob(btp::Shower);
    showerblob->AddStatus(blob_status::needs_showers);

    Cluster_Amplitude* ampl = Cluster_Amplitude::New();
    ampl->SetMS(this);
    for (int i(0);i<initialblob->NOutP();++i) {
      Particle *p(initialblob->OutParticle(i));
      if (p->GetFlow(1)==0 && p->GetFlow(2)==0) continue;
      showerblob->AddToInParticles(p);
      ColorID col(p->GetFlow(1),p->GetFlow(2));
      ampl->CreateLeg(p->Momentum(),p->Flav(),col,1<<(i+ampl->NIn()));
    }
    double mu2=initialblob->InParticle(0)->Momentum().Abs2();
    ampl->SetMuF2(mu2);
    ampl->SetKT2(mu2);
    ampl->SetMuQ2(mu2);
    showerblob->AddData("ClusterAmplitude",new Blob_Data<Cluster_Amplitude*>(ampl));
    DEBUG_VAR(*showerblob);
  }
}
