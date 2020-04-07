#include"SHERPA/PerturbativePhysics/Hard_Decay_Handler.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Phys/Weight_Info.H"
#include "ATOOLS/Math/Random.H"
#include "PHASIC++/Decays/Decay_Map.H"
#include "PHASIC++/Decays/Decay_Table.H"
#include "PHASIC++/Decays/Decay_Channel.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Single_Vertex.H"
#include "MODEL/Main/Color_Function.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "PHASIC++/Channels/Rambo.H"
#include "PHASIC++/Channels/Decay_Dalitz.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "SHERPA/SoftPhysics/Soft_Photon_Handler.H"
#include "ATOOLS/Phys/Variations.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "EXTRA_XS/One2Two/Comix1to2.H"
#include "EXTRA_XS/One2Three/Comix1to3.H"

#include <iostream>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <assert.h>

using namespace SHERPA;
using namespace ATOOLS;
using namespace MODEL;
using namespace EXTRAXS;
using namespace PHASIC;
using namespace METOOLS;
using namespace std;

Hard_Decay_Handler::Hard_Decay_Handler() :
  Decay_Handler(),
  p_newsublist(NULL), m_resultdir(""), m_offshell(""),
  m_set_widths(false),
  m_br_weights(true), m_usemass(true), m_min_prop_width(0.0)
{
  auto& s = Settings::GetMainSettings();
  auto ds = s["HARD_DECAYS"];
  /*
    TODO: Writing out a 1->3 channel which might have a large width in one
    resolved configuration, and a small one in another?
  */
  m_store_results   = ds["Store_Results"].SetDefault(0).Get<int>();
  m_br_weights      = ds["Apply_Branching_Ratios"].SetDefault(true).Get<bool>();
  int decay_tau       = ds["Decay_Tau"].SetDefault(false).Get<bool>();
  m_set_widths      = ds["Set_Widths"].SetDefault(false).Get<bool>();
  m_min_prop_width  = ds["Min_Prop_Width"].SetDefault(0.0).Get<double>();
  m_int_accuracy    = ds["Int_Accuracy"].SetDefault(0.01).Get<double>();
  m_int_niter       = ds["Int_NIter"].SetDefault(2500).Get<int>();
  m_int_target_mode = ds["Int_Target_Mode"].SetDefault(0).Get<int>();
  m_offshell        = ds["Resolve_Decays"]
    .SetDefault("Threshold")
    .UseNoneReplacements()
    .Get<std::string>();
  m_resultdir       = ds["Result_Directory"]
    .SetDefault(s["RESULT_DIRECTORY"].Get<std::string>() + "/Decays")
    .Get<std::string>();
  if (m_store_results)
    MakeDir(m_resultdir, true);

  // also need to tell shower whats massive now
  // TODO: need to use the same mass-selector
  // for now, implement alike
  SetDecayMasses();

  DEBUG_FUNC("");
  p_decaymap = new Decay_Map(this);
  KFCode_ParticleInfo_Map::const_iterator it;
  for (it=s_kftable.begin();it!=s_kftable.end();it++) {
    Flavour flav(it->first);
    if (flav.IsHadron()) continue;
    if (flav.Kfcode()==kf_tau && !decay_tau) continue;
    if (flav.IsOn() && !flav.IsStable()) {
      flav.SetDecayHandler(this);
      Decay_Table* dt=new Decay_Table(flav, this);
      vector<Decay_Table*> decaytables;
      decaytables.push_back(dt);
      p_decaymap->insert(make_pair(flav,decaytables));
      ReadDecayTable(flav);
      if (flav!=flav.Bar() && !flav.Bar().IsStable()) {
        Decay_Table* dt=new Decay_Table(flav.Bar(), this);
        vector<Decay_Table*> decaytables;
        decaytables.push_back(dt);
        p_decaymap->insert(make_pair(flav.Bar(),decaytables));
        ReadDecayTable(flav.Bar());
      }
    }
  }
  
  // initialize them sorted by masses:
  Decay_Map::iterator dmit;
  msg_Debugging()<<"Initialising hard decay tables: two-body decays.\n";
  for (dmit=p_decaymap->begin(); dmit!=p_decaymap->end(); ++dmit) {
    InitializeDirectDecays(dmit->second.at(0));
  }
  msg_Debugging()<<"Initialising hard decay tables: three-body decays.\n";
  for (dmit=p_decaymap->begin(); dmit!=p_decaymap->end(); ++dmit) {
    InitializeOffshellDecays(dmit->second.at(0));
  }
  msg_Debugging()<<"Initialising hard decay tables: customizing decay tables.\n";
  CustomizeDecayTables();

  if (m_set_widths)
    for (dmit=p_decaymap->begin(); dmit!=p_decaymap->end(); ++dmit) {
      dmit->second.at(0)->Flav().SetWidth(dmit->second.at(0)->TotalWidth());
    }

  if (p_decaymap->size()) msg_Info()<<endl<<*p_decaymap<<endl;
  WriteDecayTables();
}

Hard_Decay_Handler::~Hard_Decay_Handler()
{
  if (p_newsublist) {
    for (size_t i=0; i<p_newsublist->size(); ++i) delete (*p_newsublist)[i];
    p_newsublist->clear();
    delete p_newsublist;
  }
}

void Hard_Decay_Handler::SetDecayMasses()
{
  Settings& s = Settings::GetMainSettings();
  ATOOLS::Flavour_Vector allflavs(MODEL::s_model->IncludedFlavours());
  std::vector<size_t> defpsmassive,defpsmassless;
  const std::vector<size_t> psmassive  { s["MASSIVE_PS"].GetVector<size_t>()  };
  const std::vector<size_t> psmassless { s["MASSLESS_PS"].GetVector<size_t>() };
  const bool respect{ s["RESPECT_MASSIVE_FLAG"].Get<bool>() };
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
  // fill massive ones into m_decmass
  for (size_t i(0);i<defpsmassive.size();++i) {
    Flavour fl(defpsmassive[i],0);
    m_decmass.insert(fl);
    m_decmass.insert(fl.Bar());
  }
  Flavour_Vector mf;
  for (Flavour_Set::iterator fit(m_decmass.begin());fit!=m_decmass.end();++fit)
    if (fit->Mass(true)!=fit->Mass(false)) mf.push_back(*fit);
  msg_Info()<<METHOD<<"(): Massive decay flavours: "<<mf<<std::endl;
}


void Hard_Decay_Handler::InitializeDirectDecays(Decay_Table* dt)
{
  DEBUG_FUNC(dt->Flav());
  Flavour inflav=dt->Flav();
  Vertex_Table::const_iterator vlit=s_model->VertexTable()->find(inflav);
  const Vertex_List& vlist=(vlit!=s_model->VertexTable()->end()) ? (vlit->second) : Vertex_List();
  // temporary hack:
  // prepare reduced vertex list, not containing duplicates
  // in the sense of having identical external flavours
  Vertex_List reduced_vlist;
  for (Vertex_List::const_iterator vit=vlist.begin();vit!=vlist.end();++vit){
    Vertex_List::const_iterator vjt=reduced_vlist.begin();
    for(;vjt!=reduced_vlist.end(); ++vjt)
      if((*vjt)->in == (*vit)->in) break;
    if(vjt==reduced_vlist.end()) reduced_vlist.push_back(*vit);
  }

  msg_Debugging()<<"Vertices:"<<std::endl;
  for (size_t i=0;i<reduced_vlist.size();i++) {
    Single_Vertex* sv=reduced_vlist[i];
    if (!ProperVertex(sv)) continue;
    msg_Debugging()<<"  "<<i<<": "<<*sv<<std::endl;
    Decay_Channel* dc=new Decay_Channel(inflav, this);
    for (int j=1; j<sv->NLegs(); ++j) dc->AddDecayProduct(sv->in[j]);

    Comix1to2* diagram=new Comix1to2(dc->Flavs());
    dc->AddDiagram(diagram);

    dc->SetChannels(new Multi_Channel(""));
    dc->Channels()->SetNin(1);
    dc->Channels()->SetNout(dc->NOut());
    Rambo* rambo = new Rambo(1,dc->NOut(),&dc->Flavs().front(),this);
    dc->Channels()->Add(rambo);
    dc->Channels()->Reset();

    if (CalculateWidth(dc)) dt->AddDecayChannel(dc);
    else delete dc;
  }
  dt->UpdateWidth();
  if (m_set_widths) dt->Flav().SetWidth(dt->TotalWidth());
}

void Hard_Decay_Handler::InitializeOffshellDecays(Decay_Table* dt) {
  DEBUG_FUNC(dt->Flav()<<" "<<dt->size());
  size_t dtsize=dt->size();
  for (size_t i=0;i<dtsize;++i) {
    Decay_Channel* dc=dt->at(i);
    vector<Decay_Channel*> new_dcs=ResolveDecay(dc);
    if (TriggerOffshell(dc, new_dcs)) {
      dc->SetActive(-1);
      for (size_t j=0; j<new_dcs.size(); ++j) {
        // check for duplicates
        Decay_Channel* dup=dt->GetDecayChannel(new_dcs[j]->Flavs());
        if (dup && dup->Active()>=0) {
          DEBUG_INFO("Adding new diagram to "<<*dup);
          for (size_t k=0; k<new_dcs[j]->GetDiagrams().size(); ++k) {
            dup->AddDiagram(new_dcs[j]->GetDiagrams()[k]);
          }
          for (size_t k=0; k<new_dcs[j]->Channels()->Channels().size(); ++k) {
            dup->AddChannel(new_dcs[j]->Channels()->Channels()[k]);
          }
          new_dcs[j]->ResetDiagrams();
          new_dcs[j]->ResetChannels();
          delete new_dcs[j];
          new_dcs[j]=NULL;
          CalculateWidth(dup);
        }
        else {
          DEBUG_INFO("Adding "<<new_dcs[j]->Name());
          dt->AddDecayChannel(new_dcs[j]);
        }
      }
    }
    else {
      DEBUG_INFO("Keeping factorised.");
      for (size_t j=0; j<new_dcs.size(); ++j) {
        if (new_dcs[j]) new_dcs[j]->SetActive(-1);
      }
    }
  }
  dt->UpdateWidth();
  if (m_set_widths) dt->Flav().SetWidth(dt->TotalWidth());
}


void Hard_Decay_Handler::CustomizeDecayTables()
{
  auto s = Settings::GetMainSettings()["HARD_DECAYS"]["Channels"];
  DEBUG_FUNC(s.GetKeys().size());
  for (const auto& decay : s.GetKeys()) {
    DEBUG_VAR(decay);

    // obtain decay channel (creating it if appropriate)
    pair<Decay_Table*, Decay_Channel*> match=p_decaymap->FindDecayChannel(decay, true);
    Decay_Channel* dc = match.second;
    if (!match.first || !dc) {
      PRINT_INFO("Ignoring unknown decay channel: " << decay);
      continue;
    }

    DEBUG_VAR(dc->Name());

    // update properties
    for (const auto& propname : s[decay].GetKeys()) {
      auto propsetting = s[decay][propname];
      DEBUG_VAR(propname);
      if (propname == "Status") {
        match.first->SetChannelStatus
          (dc,propsetting.SetDefault(dc->Active()).Get<int>());
      }
      else if (propname == "Width") {
        dc->SetWidth(propsetting.SetDefault(dc->Width()).Get<double>());
        dc->SetDeltaWidth(0.0);
      }
      else {
        THROW(fatal_error,
              "Unknown HARD_DECAYS:Channels property '" + propname + "'");
      }
    }
  }
  for (Decay_Map::iterator dmit=p_decaymap->begin(); dmit!=p_decaymap->end(); ++dmit) {
    dmit->second.at(0)->UpdateWidth();
  }
}


bool Hard_Decay_Handler::TriggerOffshell(Decay_Channel* dc, vector<Decay_Channel*> new_dcs) {
  DEBUG_FUNC(dc->Name()<<"... "<<new_dcs.size());

  if (m_offshell=="Threshold") {
    double outmass=0.0;
    for (size_t j=1; j<dc->Flavs().size(); ++j)
      outmass+=this->Mass(dc->Flavs()[j]);
    DEBUG_INFO(this->Mass(dc->Flavs()[0])<<" vs "<<outmass);
    return (this->Mass(dc->Flavs()[0])<outmass);
  }
  else if (m_offshell=="ByWidth") {
    double sum_resolved_widths=0.0;
    for (size_t i=0; i<new_dcs.size(); ++i) {
      if (new_dcs[i]) sum_resolved_widths+=new_dcs[i]->Width();
    }
    return (sum_resolved_widths>dc->Width());
  }
  else if (m_offshell=="None") {
    return false;
  } else {
    THROW(fatal_error,
          "Parameter HARD_DECAYS:Resolve_Decays set to wrong value.");
  }
  return false;
}

vector<Decay_Channel*> Hard_Decay_Handler::ResolveDecay(Decay_Channel* dc1)
{
  DEBUG_FUNC(dc1->Name());
  vector<Decay_Channel*> new_dcs;
  const std::vector<ATOOLS::Flavour> flavs1(dc1->Flavs());
  for (size_t j=1;j<flavs1.size();++j) {
    bool ignore=false;
    if (flavs1[j].Width()<m_min_prop_width) continue;
    for (size_t k=1; k<j; ++k) {
      // TODO Do we really have to avoid double counting e.g. in h -> Z Z?
      // Further iterations: W+ -> b t -> b W b -> b b .. ?
      if (flavs1[j]==flavs1[k]) ignore=true;
    }
    if (ignore) continue;
    Vertex_Table::const_iterator it=s_model->VertexTable()->find(flavs1[j]);
    const Vertex_List& vlist(it->second);
    // temporary hack:
    // prepare reduced vertex list, not containing duplicates
    // in the sense of having identical external flavours
    Vertex_List reduced_vlist;
    for (Vertex_List::const_iterator vit=vlist.begin();vit!=vlist.end();++vit){
      Vertex_List::const_iterator vjt=reduced_vlist.begin();
      for(;vjt!=reduced_vlist.end(); ++vjt)
        if((*vjt)->in == (*vit)->in) break;
      if(vjt==reduced_vlist.end()) reduced_vlist.push_back(*vit);
    }

    for (size_t k=0;k<reduced_vlist.size();k++) {
      Single_Vertex* sv = reduced_vlist[k];
      if (!ProperVertex(sv)) continue;
      // TODO so far special case 1->3 only
      Decay_Channel* dc=new Decay_Channel(flavs1[0], this);
      size_t nonprop(0), propi(0), propj(0);
      dc->AddDecayProduct(flavs1[3-j]);
      dc->AddDecayProduct(sv->in[1]);
      dc->AddDecayProduct(sv->in[2]);
      DEBUG_FUNC("trying "<<dc->Name());
      // TODO what about W+ -> b t -> b W b -> b b ... two diagrams, factor 2, ...?
      // TODO what about identical particles like W' -> b t -> b W b ... two diagrams, factor 2, ...?
      // TODO what about W -> W gamma -> l v gamma?
      for (size_t l=1; l<4; ++l) {
        if (dc->Flavs()[l]==flavs1[3-j]) nonprop=l;
      }
      for (size_t l=1; l<4; ++l) {
        if (l!=nonprop && propi>0) propj=l;
        else if (l!=nonprop && propi==0) propi=l;
      }

      assert(dc1->GetDiagrams().size()==1);
      DEBUG_VAR(dc->Flavs());
      DEBUG_VAR(flavs1[j]);
      Comix1to3* diagram=new Comix1to3(dc->Flavs(),flavs1[j],
                                       nonprop, propi, propj);

      dc->AddDiagram(diagram);

      dc->SetChannels(new Multi_Channel(""));
      dc->Channels()->SetNin(1);
      dc->Channels()->SetNout(dc->NOut());
      dc->Channels()->Add(new Rambo(1,dc->NOut(),&dc->Flavs().front(),this));
      if (flavs1[j].Width()>0.0) {
        dc->Channels()->Add(new Decay_Dalitz(&dc->Flavs().front(),
                                             flavs1[j].Mass(), flavs1[j].Width(),
                                             nonprop, propi, propj, this));
      }
      dc->Channels()->Reset();

      if (CalculateWidth(dc)) new_dcs.push_back(dc);
      else delete dc;
    }
  }

  return new_dcs;
  /* TODO:
     What about cases where a 1->3 was resolved from one 1->2, but would have
     been possible from another 1->2 where it was decided not to be resolved?
  */
}

bool Hard_Decay_Handler::CalculateWidth(Decay_Channel* dc)
{
  // Integrate or use results read from decay table file
  double outmass=0.0;
  for (size_t l=1; l<dc->Flavs().size(); ++l)
    outmass+=this->Mass(dc->Flavs()[l]);
  if (this->Mass(dc->Flavs()[0])>outmass) {
    DEBUG_FUNC("Starting calculation of "<<dc->Name()<<" width now.");
    if (m_store_results && m_read.find(dc->Flavs()[0])!=m_read.end()) {
      if (m_read[dc->Flavs()[0]].find(dc->IDCode())!=
          m_read[dc->Flavs()[0]].end()) {
        const vector<double>& results(m_read[dc->Flavs()[0]][dc->IDCode()]);
        dc->SetIWidth(results[0]);
        dc->SetIDeltaWidth(results[1]);
        dc->SetMax(results[2]);
      }
      else {
        msg_Tracking()<<"    Integrating "<<dc->Name()<<endl;
        dc->CalculateWidth(m_int_accuracy,
                           m_int_target_mode==0 ? dc->Flavs()[0].Width() : 0.0,
                           m_int_niter);
      }
    }
    else {
      msg_Tracking()<<"    Integrating "<<dc->Name()<<endl;
      dc->CalculateWidth(m_int_accuracy,
                         m_int_target_mode==0 ? dc->Flavs()[0].Width() : 0.0,
                         m_int_niter);
    }
  }
  else {
    dc->SetActive(-1);
    dc->SetIWidth(0.0);
    dc->SetIDeltaWidth(0.0);
    dc->SetMax(0.0);
  }
  dc->SetWidth(dc->IWidth());
  dc->SetDeltaWidth(dc->IDeltaWidth());
  return true;
}

bool Hard_Decay_Handler::ProperVertex(MODEL::Single_Vertex* sv)
{
  if (sv->dec) return false;

  for (int i(0); i<sv->NLegs(); ++i)
    if (sv->in[i].IsDummy()) return false;

  if (sv->NLegs()!=3) return false; // TODO

  // TODO: ignore radiation graphs. should we?
  for (int i=1; i<sv->NLegs(); ++i) {
    if (sv->in[i].Kfcode()==sv->in[0].Kfcode()) {
      return false;
    }
  }

  // what about extra particles like Z4 if Z stable?

  return true;
}


void Hard_Decay_Handler::CreateDecayBlob(Blob_List* bloblist, ATOOLS::Particle* inpart)
{
  DEBUG_FUNC(inpart->Flav());
  Blob* blob = bloblist->AddBlob(btp::Hard_Decay);
  blob->AddStatus(blob_status::needs_extraQED);
  blob->AddToInParticles(inpart);
  blob->SetTypeSpec("Sherpa");
  Decay_Table* table=p_decaymap->FindDecay(blob->InParticle(0)->Flav());
  if (table==NULL) {
    msg_Error()<<METHOD<<" decay table not found, retrying event."<<endl
               <<*blob<<endl;
    throw Return_Value::Retry_Event;
  }
  blob->AddData("dc",new Blob_Data<Decay_Channel*>(table->Select()));

  DEBUG_INFO("p_onshell="<<inpart->Momentum());
  blob->AddData("p_onshell",new Blob_Data<Vec4D>(inpart->Momentum()));
  DEBUG_INFO("succeeded.");
}

void Hard_Decay_Handler::FindDecayProducts(Particle* decayer,
                                           list<Particle*>& decayprods)
{
  if (decayer->DecayBlob()==NULL) {
    decayprods.push_back(decayer);
  }
  else {
    for (size_t i=0; i<decayer->DecayBlob()->NOutP(); ++i) {
      FindDecayProducts(decayer->DecayBlob()->OutParticle(i), decayprods);
    }
  }
}

double Hard_Decay_Handler::BRFactor(ATOOLS::Blob* blob) const
{
  double brfactor=1.0;
  for (size_t i=0; i<blob->NOutP(); ++i) {
    Particle* part=blob->OutParticle(i);
    Decay_Table* dt=p_decaymap->FindDecay(part->RefFlav());
    if (dt) {
      brfactor*=dt->ActiveWidth()/dt->TotalWidth();
      if (part->DecayBlob() && part->DecayBlob()->Type()==btp::Hard_Decay)
        brfactor*=BRFactor(part->DecayBlob());
    }
  }
  return brfactor;
}

void Hard_Decay_Handler::AfterTreatInitialBlob(Blob* blob, Blob_List* bloblist)
{
  if (blob->Type()!=btp::Signal_Process) return;

  double brfactor=m_br_weights ? BRFactor(blob) : 1.0;
  DEBUG_VAR(brfactor);
  Blob_Data_Base * bdbmeweight((*blob)["MEWeight"]);
  if (bdbmeweight) {
    bdbmeweight->Set<double>(brfactor*bdbmeweight->Get<double>());
  }
  Blob_Data_Base * wgtinfo((*blob)["MEWeightInfo"]);
  if (wgtinfo) *wgtinfo->Get<ME_Weight_Info*>()*=brfactor;

  Blob_Data_Base * weights((*blob)["Weights"]);
  if (weights) weights->Get<Event_Weights>()*=brfactor;

  NLO_subevtlist* sublist(NULL);
  Blob_Data_Base * bdb((*blob)["NLO_subeventlist"]);
  if (bdb) sublist=bdb->Get<NLO_subevtlist*>();
  if (sublist) {
    // If the blob contains a NLO_subeventlist, we have to attach decays
    // in the sub-events as well. The decay has to be identical for infrared
    // cancellations, so we simply take each decay and boost it to the sub
    // kinematics to replace the particle in the subevent

    DEBUG_FUNC("");

    vector<list<Particle*> > decayprods(blob->NOutP());
    size_t newn(2);
    for (size_t i=0; i<blob->NOutP()-1; ++i) {
      // iterate over out-particles excluding real emission parton
      list<Particle*> decayprods_i;
      FindDecayProducts(blob->OutParticle(i), decayprods_i);
      DEBUG_VAR(blob->OutParticle(i)->Flav());
      list<Particle*>::const_iterator it;
      for (it=decayprods_i.begin(); it!=decayprods_i.end(); ++it) {
        DEBUG_VAR((*it)->Flav());
      }
      decayprods[i]=decayprods_i;
      newn+=decayprods_i.size();
    }
    DEBUG_VAR(newn);

    if (p_newsublist) {
      for (size_t i=0; i<p_newsublist->size(); ++i) delete (*p_newsublist)[i];
      p_newsublist->clear();
    }
    else p_newsublist=new NLO_subevtlist();
    for (size_t i=0; i<sublist->size(); ++i) {
      // iterate over sub events and replace decayed particles
      NLO_subevt* sub((*sublist)[i]);
      DEBUG_VAR(*(*sublist)[i]);

      if (sub->IsReal()) newn+=1;
      Flavour* newfls = new Flavour[newn];
      Vec4D* newmoms = new Vec4D[newn];
      size_t* newid = new size_t[newn];
      for (size_t n=0; n<newn; ++n) newid[n]=0;
      NLO_subevt* newsub=new NLO_subevt(*sub);
      newsub->m_n=newn;
      newsub->p_id=newid;
      newsub->p_fl=newfls;
      newsub->p_mom=newmoms;
      newsub->m_delete=true;
      p_newsublist->push_back(newsub);

      int nin=blob->NInP();
      for (size_t j=0, jnew=0; j<nin+blob->NOutP()-1; ++j, ++jnew) {
        if (j<2 || decayprods[j-nin].size()==1) {
          newfls[jnew]=sub->p_fl[j];
          newmoms[jnew]=sub->p_mom[j];
          continue;
        }
        if (sub->p_fl[j]!=blob->OutParticle(j-nin)->Flav()) {
          THROW(fatal_error, "Internal Error 1");
        }

        Vec4D oldmom=blob->OutParticle(j-nin)->Momentum();
        Vec4D newmom=sub->p_mom[j];
        Poincare cms(oldmom);
        Poincare newframe(newmom);
        newframe.Invert();

        list<Particle*>::const_iterator it;
        for (it=decayprods[j-nin].begin(); it!=decayprods[j-nin].end(); ++it) {
          newfls[jnew]=(*it)->Flav();
          newmoms[jnew]=Vec4D(newframe*(cms*(*it)->Momentum()));
          jnew++;
        }
        jnew--; // because we replaced one particle
      }
      if (sub->IsReal()) {
        // Add remaining parton for real event
        newfls[newn-1]=sub->p_fl[sub->m_n-1];
        newmoms[newn-1]=sub->p_mom[sub->m_n-1];
      }
    }
    bdb->Set<NLO_subevtlist*>(p_newsublist);
    DEBUG_INFO("New subevts:");
    for (size_t i=0;i<p_newsublist->size();++i) {
      (*p_newsublist)[i]->m_result*=brfactor;
      (*p_newsublist)[i]->m_results*=brfactor;
      (*p_newsublist)[i]->m_me*=brfactor;
      (*p_newsublist)[i]->m_mewgt*=brfactor;
      DEBUG_VAR(*(*p_newsublist)[i]);
    }
  }
}

void Hard_Decay_Handler::ReadDecayTable(Flavour decayer)
{
  DEBUG_FUNC(decayer);
  if (!m_store_results) return;
  Data_Reader reader = Data_Reader("|",";","!");
  reader.AddComment("#");
  reader.AddComment("//");
  reader.AddWordSeparator("\t");
  reader.SetInputPath(m_resultdir);
  reader.SetInputFile(decayer.ShellName());

  vector<vector<string> > file;
  if(reader.MatrixFromFile(file)) {
    for (size_t iline=0; iline<file.size(); ++iline) {
      if (file[iline].size()==4) {
        string decaychannel=file[iline][0];
        vector<double> results(3);
        for (size_t i=0; i<3; ++i) results[i]=ToType<double>(file[iline][i+1]);
        m_read[decayer].insert(make_pair(decaychannel, results));
      }
      else {
        PRINT_INFO("Wrong format in decay table in "<<m_resultdir);
      }
    }
  }
}

void Hard_Decay_Handler::WriteDecayTables()
{
  if (!(m_store_results & 1)) return;

  Decay_Map::iterator dmit;
  for (dmit=p_decaymap->begin(); dmit!=p_decaymap->end(); ++dmit) {
    ofstream ostr((m_resultdir+"/"+dmit->first.ShellName()).c_str());
    ostr<<"# Decay table for "<<dmit->first<<endl;
    ostr<<"# IDCode                   \tWidth     \tDeltaWidth \tMaximum"<<endl<<endl;
    Decay_Table::iterator dtit;
    for (dtit=dmit->second[0]->begin(); dtit!=dmit->second[0]->end(); ++dtit) {
      ostr<<setw(25)<<left<<(*dtit)->IDCode()<<"\t"
          <<setw(12)<<left<<(*dtit)->IWidth()<<"\t"
          <<setw(12)<<left<<(*dtit)->IDeltaWidth()<<"\t"
          <<setw(12)<<left<<(*dtit)->Max()<<endl;
    }
    ostr.close();
  }
}

double Hard_Decay_Handler::Mass(const ATOOLS::Flavour &fl) const
{
  if (m_usemass==0) return fl.Mass();
  if (m_decmass.find(fl)!=m_decmass.end()) return fl.Mass(true);
  return fl.Mass();
}
