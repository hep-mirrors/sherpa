#include"SHERPA/PerturbativePhysics/Hard_Decay_Handler.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Math/Random.H"
#include "PHASIC++/Decays/Decay_Map.H"
#include "PHASIC++/Decays/Decay_Table.H"
#include "PHASIC++/Decays/Decay_Channel.H"

#include "MODEL/Main/Model_Base.H"
#include "MODEL/Interaction_Models/Single_Vertex.H"
#include "MODEL/Interaction_Models/Color_Function.H"
#include "PHASIC++/Decays/Color_Function_Decay.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "PHASIC++/Channels/Rambo.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "SHERPA/SoftPhysics/Soft_Photon_Handler.H"

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

class ParticlePairFirstEnergySort {
public:
  bool operator()(const ParticlePair& a,const ParticlePair& b)
  { return (a.first->Momentum()[0]<b.first->Momentum()[0]); }
};

Hard_Decay_Handler::Hard_Decay_Handler(std::string path, std::string file) :
  p_newsublist(NULL), m_path(""), m_file(""), m_resultdir(""), m_offshell(""),
  m_store_results(false), m_decay_tau(false), m_set_widths(false),
  m_br_weights(true), m_usemass(true)
{
  Data_Reader dr(" ",";","!","=");
  dr.AddWordSeparator("\t");
  dr.SetInputPath(path);
  dr.SetInputFile(file);
  m_mass_smearing=dr.GetValue<int>("HARD_MASS_SMEARING",1);
  m_spincorr=rpa->gen.HardSC();
  /*
    TODO: Writing out a 1->3 channel which might have a large width in one
    resolved configuration, and a small one in another?
  */
  m_store_results=dr.GetValue<int>("STORE_DECAY_RESULTS",0);
  m_br_weights=dr.GetValue<int>("HDH_BR_WEIGHTS",1);
  m_decay_tau=dr.GetValue<int>("DECAY_TAU_HARD",0);
  m_set_widths=dr.GetValue<int>("HDH_SET_WIDTHS",0);
  // also need to tell shower whats massive now
  // TODO: need to use the same mass-selector
  // for now, implement alike
  SetDecayMasses(&dr);
  m_qedmode=dr.GetValue<int>("HDH_QED_CORRECTIONS",1);
  if (m_qedmode && m_spincorr) m_qedmode=2;
  if (m_qedmode>0 && !m_usemass) {
    THROW(fatal_error,std::string("QED corrections to hard decays only ")
                      +std::string("available in massive mode."));
  }
  std::string r=dr.GetValue<std::string>("DECAY_RESULT_DIRECTORY","");
  m_resultdir=(r==""?dr.GetValue<std::string>("RESULT_DIRECTORY","Results"):r);
  if (m_store_results) {
    MakeDir(m_resultdir+"/Decays/", true);
  }
  m_offshell=dr.GetValue<std::string>("RESOLVE_DECAYS", "Threshold");

  Data_Reader nodecayreader("|", ";", "#");
  nodecayreader.SetString(dr.GetValue<std::string>("HDH_NO_DECAY", ""));
  vector<string> disabled_channels;
  nodecayreader.VectorFromString(disabled_channels);
  for (size_t i=0; i<disabled_channels.size(); ++i)
    m_disabled_channels.insert(disabled_channels[i]);
  nodecayreader.SetString(dr.GetValue<std::string>("HDH_ONLY_DECAY", ""));
  vector<string> forced_channels;
  nodecayreader.VectorFromString(forced_channels);
  for (size_t i=0; i<forced_channels.size(); ++i) {
    string fc(forced_channels[i]);
    size_t bracket(fc.find("{")), comma(fc.find(","));
    int kfc(atoi((fc.substr(bracket+1,comma-bracket-1)).c_str()));
    m_forced_channels[Flavour(abs(kfc), kfc<0)].insert(fc);
  }

  DEBUG_FUNC("");
  p_decaymap = new Decay_Map(this);
  KFCode_ParticleInfo_Map::const_iterator it;
  for (it=s_kftable.begin();it!=s_kftable.end();it++) {
    Flavour flav(it->first);
    if (Decays(flav)) {
      Decay_Table* dt=new Decay_Table(flav, this);
      vector<Decay_Table*> decaytables;
      decaytables.push_back(dt);
      p_decaymap->insert(make_pair(flav,decaytables));
      ReadDecayTable(flav);
    }
    if (flav!=flav.Bar() && Decays(flav.Bar())) {
      Decay_Table* dt=new Decay_Table(flav.Bar(), this);
      vector<Decay_Table*> decaytables;
      decaytables.push_back(dt);
      p_decaymap->insert(make_pair(flav.Bar(),decaytables));
      ReadDecayTable(flav.Bar());
    }
  }
  
  // initialize them sorted by masses:
  Decay_Map::iterator dmit;
  Vertex_Table offshell;
  msg_Info()<<"Initialising hard decay tables."<<endl;
  size_t i(0);
  for (dmit=p_decaymap->begin(); dmit!=p_decaymap->end(); ++dmit) {
    offshell.insert(make_pair(dmit->first,Vertex_List()));
    msg_Tracking()<<"  Initialising two-body decays. Step "
		  <<++i<<"/"<<p_decaymap->size()<<" ("<<dmit->first<<")            "
		  <<endl;
    InitializeDirectDecays(dmit->second.at(0));
  }
  i=0;
  if (p_decaymap->size()) msg_Tracking()<<endl;
  for (dmit=p_decaymap->begin(); dmit!=p_decaymap->end(); ++dmit) {
    msg_Tracking()<<"  Initialising three-body decays. Step "
		  <<++i<<"/"<<p_decaymap->size()<<" ("<<dmit->first<<")            "
		  <<endl;
    if (m_offshell=="None") {
      PRINT_INFO("Warning: Ignoring offshell decays as requested.");
    }
    else if (m_offshell=="Threshold") {
      RefineDecaysThreshold(dmit->second.at(0));
    }
    else if (m_offshell=="ByWidth") {
      RefineDecaysByWidth(dmit->second.at(0));
    }
    else THROW(fatal_error, "Parameter RESOLVE_DECAYS set to wrong value.");

    // force/disable specified decay channels
    for (size_t i=0;i<dmit->second.at(0)->size();++i) {
      Decay_Channel* dc=dmit->second.at(0)->at(i);
      if (dc->Active()<1) continue;
      string idc(dc->IDCode());
      size_t bracket(idc.find("{")), comma(idc.find(","));
      int kfc(atoi((idc.substr(bracket+1,comma-bracket-1)).c_str()));
      Flavour dec(abs(kfc), kfc<0);
      if (m_disabled_channels.count(idc) ||
	  (m_forced_channels.find(dec)!=m_forced_channels.end() &&
           m_forced_channels[dec].size() &&
	   !m_forced_channels[dec].count(idc))) 
	dc->SetActive(0);
    }

    dmit->second.at(0)->UpdateWidth();
    if (m_set_widths)
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

void Hard_Decay_Handler::SetDecayMasses(Data_Reader *const dr)
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
  const Vertex_List& vertexlist(vlit->second);

  msg_Debugging()<<"Vertices:"<<std::endl;
  for (size_t i=0;i<vertexlist.size();i++) {
    Single_Vertex* sv=vertexlist[i];
    if (!ProperVertex(sv)) continue;
    msg_Debugging()<<"  "<<i<<": "<<*sv<<std::endl;
    Decay_Channel* dc=new Decay_Channel(inflav, this);
    for (int j=1; j<sv->nleg; ++j) dc->AddDecayProduct(sv->in[j]);

    assert(sv->Color.size()==1);
    Comix1to2* diagram=new Comix1to2(dc->Flavs());
    dc->AddDiagram(diagram,new Color_Function_Decay(sv->Color[0]));

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

void Hard_Decay_Handler::RefineDecaysThreshold(Decay_Table* dt) {
  DEBUG_FUNC(dt->Flav()<<" "<<dt->size());
  size_t dtsize=dt->size();
  for (size_t i=0;i<dtsize;++i) {
    Decay_Channel* dc=dt->at(i);
    if (!dc) THROW(fatal_error,"Could not retrieve decay table.");
    msg_Debugging()<<*dc<<std::endl;
    double outmass=0.0;
    for (size_t j=1; j<dc->Flavs().size(); ++j)
      outmass+=this->Mass(dc->Flavs()[j]);
    msg_Debugging()<<this->Mass(dt->Flav())<<" vs "<<outmass<<std::endl;
    if (this->Mass(dt->Flav())>outmass) continue;

    vector<Decay_Channel*> new_dcs=ResolveDecay(dc);
    for (size_t j=0; j<new_dcs.size(); ++j) {
      dt->AddDecayChannel(new_dcs[j]);
    }
    if (new_dcs.size()>0) {
      dc->SetActive(-1);
      for (size_t j=0; j<new_dcs.size(); ++j) {
        DEBUG_INFO("Adding "<<*new_dcs[j]);
      }
    }
    else {
      for (size_t j=0; j<new_dcs.size(); ++j) {
        new_dcs[j]->SetActive(-1);
      }
    }
  }
}

void Hard_Decay_Handler::RefineDecaysByWidth(Decay_Table* dt) {
  DEBUG_FUNC(dt->Flav()<<" "<<dt->size());
  size_t dtsize=dt->size();
  for (size_t i=0;i<dtsize;++i) {
    Decay_Channel* dc=dt->at(i);
    DEBUG_VAR(*dc);
    vector<Decay_Channel*> new_dcs=ResolveDecay(dc);
    double sum_resolved_widths=0.0;
    for (size_t j=0; j<new_dcs.size(); ++j) {
      Decay_Channel* dup=dt->GetDecayChannel(new_dcs[j]->Flavs());
      if (dup) {
        DEBUG_INFO("Duplicate for "<<*dup);
        for (DiagColVec::const_iterator it=new_dcs[j]->GetDiagrams().begin();
             it!=new_dcs[j]->GetDiagrams().end(); ++it) {
          dup->AddDiagram(it->first, it->second);
        }
        for (vector<Single_Channel*>::const_iterator it=new_dcs[j]->Channels()->Channels().begin();
             it!=new_dcs[j]->Channels()->Channels().end(); ++it) {
          dup->AddChannel(*it);
        }
        new_dcs[j]->ResetDiagrams();
        new_dcs[j]->ResetChannels();
        delete new_dcs[j];
        new_dcs[j]=NULL;
        sum_resolved_widths-=dup->Width();
        CalculateWidth(dup);
        sum_resolved_widths+=dup->Width();
      }
      else {
        sum_resolved_widths+=new_dcs[j]->Width();
        dt->AddDecayChannel(new_dcs[j]);
      }
    }
    DEBUG_INFO("resolved="<<sum_resolved_widths<<" vs. "<<dc->Width());

    if (sum_resolved_widths>dc->Width()) {
      dc->SetActive(-1);
    }
    else {
      for (size_t j=0; j<new_dcs.size(); ++j) {
        if (new_dcs[j]) new_dcs[j]->SetActive(-1);
      }
    }
  }
}

vector<Decay_Channel*> Hard_Decay_Handler::ResolveDecay(Decay_Channel* dc1)
{
  DEBUG_FUNC(*dc1);
  vector<Decay_Channel*> new_dcs;
  const std::vector<ATOOLS::Flavour> flavs1(dc1->Flavs());
  for (size_t j=1;j<flavs1.size();++j) {
    bool ignore=false;
    for (size_t k=1; k<j; ++k) {
      // TODO Do we really have to avoid double counting e.g. in h -> Z Z?
      // Further iterations: W+ -> b t -> b W b -> b b .. ?
      if (flavs1[j]==flavs1[k]) ignore=true;
    }
    if (ignore) continue;
    Vertex_Table::const_iterator it=s_model->VertexTable()->find(flavs1[j]);
    const Vertex_List& vertexlist(it->second);
    for (size_t k=0;k<vertexlist.size();k++) {
      Single_Vertex* sv = vertexlist[k];
      if (!ProperVertex(sv)) continue;
      // TODO so far special case 1->3 only
      Decay_Channel* dc=new Decay_Channel(flavs1[0], this);
      size_t nonprop(0), propi(0), propj(0);
      dc->AddDecayProduct(flavs1[3-j]);
      dc->AddDecayProduct(sv->in[1]);
      dc->AddDecayProduct(sv->in[2]);
      DEBUG_FUNC("trying "<<*dc);
      // TODO what about W+ -> b t -> b W b -> b b ... two diagrams, factor 2, ...?
      // TODO what about W' -> b t -> b W b ... two diagrams, factor 2, ...?
      for (size_t l=1; l<4; ++l) {
        if (dc->Flavs()[l]==flavs1[3-j]) nonprop=l;
      }
      for (size_t l=1; l<4; ++l) {
        if (l!=nonprop && propi>0) propj=l;
        else if (l!=nonprop && propi==0) propi=l;
      }

      assert(dc1->GetDiagrams().size()==1);
      assert(sv->Color.size()==1);
      DEBUG_VAR(dc->Flavs());
      DEBUG_VAR(flavs1[j]);
      Comix1to3* diagram=new Comix1to3(dc->Flavs(),flavs1[j],
                                       nonprop, propi, propj);
      Color_Function_Decay* col1=new Color_Function_Decay(*dc1->GetDiagrams()[0].second);
      DEBUG_VAR(*col1);
      Color_Function_Decay col2(sv->Color[0]);
      DEBUG_VAR(col2);
      vector<int> bumps = col1->Multiply(col2);
      DEBUG_VAR(*col1);
      DEBUG_INFO("Contracting "<<0+bumps[0]<<" with "<<j);
      col1->Contract(0+bumps[0],j);
      DEBUG_VAR(*col1);
      dc->AddDiagram(diagram, col1);

      dc->SetChannels(new Multi_Channel(""));
      dc->Channels()->SetNin(1);
      dc->Channels()->SetNout(dc->NOut());
      Rambo* rambo = new Rambo(1,dc->NOut(),&dc->Flavs().front(),this);
      dc->Channels()->Add(rambo);
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
    DEBUG_INFO("Starting calculation of widths now.");
    DEBUG_INFO(*dc);
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
        dc->CalculateWidth();
      }
    }
    else {
      msg_Tracking()<<"    Integrating "<<dc->Name()<<endl;
      dc->CalculateWidth();
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
  if (!sv->on || sv->dec) return false;

  for (int i(0); i<sv->nleg; ++i)
    if (sv->in[i].IsDummy()) return false;

  if (sv->nleg!=3) return false; // TODO

  // ignore radiation graphs. should we?
  for (int i=1; i<sv->nleg; ++i) {
    if (sv->in[i].Kfcode()==sv->in[0].Kfcode()) {
      return false;
    }
  }

  // what about extra particles like Z4 if Z stable?

  return true;
}


void Hard_Decay_Handler::CreateDecayBlob(ATOOLS::Particle* inpart)
{
  DEBUG_FUNC(inpart->Flav());
  if(inpart->DecayBlob()) THROW(fatal_error,"Decay blob already exists.");
  if(!Decays(inpart->Flav())) return;
  Blob* blob = p_bloblist->AddBlob(btp::Hard_Decay);
  blob->SetStatus(blob_status::needs_showers);
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

void Hard_Decay_Handler::TreatInitialBlob(ATOOLS::Blob* blob,
                                          METOOLS::Amplitude2_Tensor* amps,
                                          const Particle_Vector& origparts)
{
  Decay_Handler_Base::TreatInitialBlob(blob, amps, origparts);

  double brfactor=m_br_weights ? BRFactor(blob) : 1.0;
  DEBUG_VAR(brfactor);
  Blob_Data_Base * bdbweight((*blob)["Weight"]);
  if (bdbweight) bdbweight->Set<double>(brfactor*bdbweight->Get<double>());
  Blob_Data_Base * bdbmeweight((*blob)["MEWeight"]);
  if (bdbmeweight) {
    // msg_Out()<<METHOD<<"(ME = "<<bdbmeweight->Get<double>()<<", "
    // 	     <<"BR = "<<brfactor<<") for\n"<<"   "
    // 	     <<blob->InParticle(0)->Flav()<<" "
    // 	     <<blob->InParticle(1)->Flav()<<" -->";
    // for (size_t i=0;i<blob->NOutP();i++) 
    //   msg_Out()<<" "<<blob->OutParticle(i)->Flav();
    // msg_Out()<<".\n";
    bdbmeweight->Set<double>(brfactor*bdbmeweight->Get<double>());
  }
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
      (*p_newsublist)[i]->m_me*=brfactor;
      (*p_newsublist)[i]->m_mewgt*=brfactor;
      DEBUG_VAR(*(*p_newsublist)[i]);
    }
  }
}

bool Hard_Decay_Handler::DefineInitialConditions(Cluster_Amplitude* ampl,
                                                 Blob* initial_blob)
{
  DEBUG_FUNC(this);
  DEBUG_VAR(*ampl);
  for (int i=0; i<initial_blob->NOutP(); ++i) {
    ampl->Leg(initial_blob->NInP()+i)->SetMom
      (initial_blob->OutParticle(i)->Momentum());
  }
  if (p_clus->ReCluster(ampl)<0) return false;
  if (ampl->NIn()==2) {
    for (Cluster_Amplitude *campl(ampl);
	 campl;campl=campl->Next()) {
      if (-campl->Leg(0)->Mom()[0]>rpa->gen.PBeam(0)[0] ||
	  -campl->Leg(1)->Mom()[0]>rpa->gen.PBeam(1)[0])
	return false;
    }
  }
  size_t imax=ampl->Legs().size()-1;
  for (int i=0; i<initial_blob->NOutP(); ++i) {
    if (initial_blob->OutParticle(i)->DecayBlob()) {
      AddDecayClustering(ampl, initial_blob->OutParticle(i)->DecayBlob(),
                         imax, 1<<(initial_blob->NInP()+i));
    }
  }
  return true;
}

void Hard_Decay_Handler::AddDecayClustering(ATOOLS::Cluster_Amplitude*& ampl,
                                            Blob* blob,
                                            size_t& imax,
                                            size_t idmother)
{
  DEBUG_FUNC("blob->Id()="<<blob->Id()<<" idmother="<<ID(idmother));
  DEBUG_VAR(*blob);
  Particle_Vector daughters;
  ParticlePair_Vector photons;
  for (size_t i(0);i<blob->GetOutParticles().size();++i) {
    Particle * p(blob->OutParticle(i));
    if (p->Info()=='S') photons.push_back(make_pair(p,p));
    else                daughters.push_back(p);
  }
  std::sort(photons.begin(),photons.end(),ParticlePairFirstEnergySort());
  msg_Debugging()<<"daughters: ";
  for (size_t i(0);i<daughters.size();++i)
    msg_Debugging()<<daughters[i]->Flav().IDName()<<" ";
  msg_Debugging()<<" +  "<<photons.size()<<" soft photon(s)"<<std::endl;
  AssignPhotons(daughters,photons);
  if (daughters.size()==2) {
    msg_Debugging()<<"1 to 2 case"<<std::endl;
    Cluster_Amplitude* copy=ampl->InitPrev();
    copy->CopyFrom(ampl);
    copy->SetMS(ampl->MS());
    Cluster_Leg *lij(ampl->IdLeg(idmother));
    copy->SetKT2(lij->Mom().Abs2());
    for (size_t i=0; i<ampl->Legs().size(); ++i)
      ampl->Leg(i)->SetStat(ampl->Leg(i)->Stat()|1);
    lij->SetStat(1|2|4);
    size_t idk(0);
    for (size_t i=0; i<copy->Legs().size(); ++i) {
      copy->Leg(i)->SetK(0);
      if (copy->Leg(i)->Id()!=idmother &&
          (copy->Leg(i)->Col().m_i==lij->Col().m_j ||
           copy->Leg(i)->Col().m_j==lij->Col().m_i))
        idk=copy->Leg(i)->Id();
    }
    if (lij->Col().m_i==0 && lij->Col().m_j==0) {
      // Ad hoc EW partner
      size_t ampl_nout=ampl->Legs().size()-ampl->NIn();
      if (ampl_nout==1) idk=ampl->Leg(0)->Id();
      else {
        size_t select=ampl->Legs().size();
        do {
          select=ampl->NIn()+floor(ran->Get()*ampl_nout);
        } while (ampl->Leg(select)->Id()&idmother ||
                 select>ampl->Legs().size()-1);
        idk=ampl->Leg(select)->Id();
      }
    }
    if (idk==0) THROW(fatal_error,"Colour partner not found");
    lij->SetK(idk);
    Cluster_Leg *d1(copy->IdLeg(idmother));
    size_t stat1(0), stat2(0);
    d1->SetMom(RecombinedMomentum(daughters[0],photons,stat1));
    d1->SetStat(stat1);
    d1->SetFlav(daughters[0]->Flav());
    copy->CreateLeg(RecombinedMomentum(daughters[1],photons,stat2),
                    daughters[1]->RefFlav());
    size_t idnew=1<<(++imax);
    copy->Legs().back()->SetId(idnew);
    copy->Legs().back()->SetStat(stat2);
    Cluster_Amplitude::SetColours(ampl->IdLeg(idmother),
                                  copy->IdLeg(idmother),
                                  copy->Legs().back());
    DEBUG_VAR(*copy);
    Cluster_Amplitude* tmp=copy;
    while (tmp->Next()) {
      tmp=tmp->Next();
      for (size_t i=0; i<tmp->Legs().size(); ++i) {
        if (tmp->Leg(i)->Id()&idmother) {
          tmp->Leg(i)->SetId(tmp->Leg(i)->Id()|idnew);
	}
        if (tmp->Leg(i)->K()&idmother) {
          tmp->Leg(i)->SetK(tmp->Leg(i)->K()|idnew);
        }
      }
      DEBUG_VAR(*tmp);
    }
    std::vector<size_t> ids;
    ids.push_back(idmother);
    ids.push_back(idnew);
    while (photons.size())
      AddPhotonsClustering(copy, daughters, photons, imax, ids);
    if (daughters[0]->DecayBlob())
      AddDecayClustering(copy, daughters[0]->DecayBlob(), imax, idmother);
    if (daughters[1]->DecayBlob())
      AddDecayClustering(copy, daughters[1]->DecayBlob(), imax, idnew);
    ampl=copy;
  }
  else if (daughters.size()==3) {
    msg_Debugging()<<"1 to 3 case"<<std::endl;
    // structure m -> 0 P[->1 2]
    // propagator always combines daughters 1+2
    Cluster_Amplitude* step1=ampl->InitPrev();
    step1->CopyFrom(ampl);
    step1->SetMS(ampl->MS());
    Cluster_Leg *lij(ampl->IdLeg(idmother));
    step1->SetKT2(lij->Mom().Abs2());
    if (!lij) THROW(fatal_error,"Cluster leg of id "+ToString(idmother)
                                +" not found.");
    for (size_t i=0; i<ampl->Legs().size(); ++i)
      ampl->Leg(i)->SetStat(ampl->Leg(i)->Stat()|1);
    lij->SetStat(1|2|4);
    size_t idk(0);
    for (size_t i=0; i<step1->Legs().size(); ++i) {
      step1->Leg(i)->SetK(0);
      if (step1->Leg(i)->Id()!=idmother)
	if (step1->Leg(i)->Col().m_i==lij->Col().m_j ||
	    step1->Leg(i)->Col().m_j==lij->Col().m_i) 
	  idk=step1->Leg(i)->Id();
    }
    if (lij->Col().m_i==0 && lij->Col().m_j==0) {
      // Ad hoc EW partner
      size_t ampl_nout=ampl->Legs().size()-ampl->NIn();
      if (ampl_nout==1) idk=ampl->Leg(0)->Id();
      else {
        size_t select=ampl->Legs().size();
        do {
          select=ampl->NIn()+floor(ran->Get()*ampl_nout);
        } while (ampl->Leg(select)->Id()&idmother || select>ampl->Legs().size()-1);
        idk=ampl->Leg(select)->Id();
      }
    }
    if (idk==0) THROW(fatal_error,"Colour partner not found");
    lij->SetK(idk);
    Cluster_Leg *d1(step1->IdLeg(idmother));
    size_t stat1(0),stat2(0),stat3(0);
    d1->SetMom(RecombinedMomentum(daughters[0],photons,stat1));
    d1->SetStat(stat1);
    d1->SetFlav(daughters[0]->Flav());
    // todo: 1->2 qcd shower with ew fs recoil partner
    // d1->SetK(idmother);// not that simple: w->qq' has color connection in fs
    Decay_Channel* dc(NULL);
    Blob_Data_Base* data = (*blob)["dc"];
    if (data) {
      dc=data->Get<Decay_Channel*>();
      DEBUG_VAR(*dc);
    }
    else THROW(fatal_error, "Internal error.");
    Comix1to3* amp=dynamic_cast<Comix1to3*>(dc->GetDiagrams()[0].first);
    if (!amp) THROW(fatal_error, "Internal error.");
    Flavour prop_flav=amp->Prop();
    Vec4D momd2=RecombinedMomentum(daughters[1],photons,stat2);
    Vec4D momd3=RecombinedMomentum(daughters[2],photons,stat3);
    Vec4D prop_mom=momd2+momd3;
    step1->CreateLeg(prop_mom, prop_flav);
    size_t idnew1=1<<(++imax);
    step1->Legs().back()->SetId(idnew1);
    step1->Legs().back()->SetStat(0);
    Cluster_Amplitude::SetColours(ampl->IdLeg(idmother),
                                  step1->IdLeg(idmother),
                                  step1->Legs().back());
    DEBUG_VAR(*step1);
    Cluster_Amplitude* tmp=step1;
    while (tmp->Next()) {
      tmp=tmp->Next();
      for (size_t i=0; i<tmp->Legs().size(); ++i) {
        if (tmp->Leg(i)->Id()&idmother) {
          tmp->Leg(i)->SetId(tmp->Leg(i)->Id()|idnew1);
	}
        if (tmp->Leg(i)->K()&idmother) {
          tmp->Leg(i)->SetK(tmp->Leg(i)->K()|idnew1);
        }
      }
      DEBUG_VAR(*tmp);
    }

    
    Cluster_Amplitude* step2=step1->InitPrev();
    step2->CopyFrom(step1);
    step2->SetMS(step1->MS());
    for (size_t i=0; i<step1->Legs().size(); ++i)
      step1->Leg(i)->SetStat(step1->Leg(i)->Stat()|1);
    step1->IdLeg(idnew1)->SetStat(1|4);
    step1->IdLeg(idnew1)->SetK(idk);
    for (size_t i=0; i<step2->Legs().size(); ++i) step2->Leg(i)->SetK(0);
    Cluster_Leg *d2(step2->IdLeg(idnew1));
    d2->SetMom(momd2);
    d2->SetStat(stat2);
    d2->SetFlav(daughters[1]->Flav());
    step2->CreateLeg(momd3, daughters[2]->Flav());
    size_t idnew2=1<<(++imax);
    step2->Legs().back()->SetId(idnew2);
    step2->Legs().back()->SetStat(stat3);
    Cluster_Amplitude::SetColours(step1->IdLeg(idnew1),
                                  step2->IdLeg(idnew1),
                                  step2->Legs().back());
    DEBUG_VAR(*step2);
    tmp=step2;
    while (tmp->Next()) {
      tmp=tmp->Next();
      for (size_t i=0; i<tmp->Legs().size(); ++i) {
        if (tmp->Leg(i)->Id()&idnew1) {
          tmp->Leg(i)->SetId(tmp->Leg(i)->Id()|idnew2);
	}
        if (tmp->Leg(i)->K()&idnew1) {
          tmp->Leg(i)->SetK(tmp->Leg(i)->K()|idnew2);
        }
      }
      DEBUG_VAR(*tmp);
    }

    std::vector<size_t> ids;
    ids.push_back(idmother);
    ids.push_back(idnew1);
    ids.push_back(idnew2);
    while (photons.size())
      AddPhotonsClustering(step2,daughters,photons,imax,ids);
    if (daughters[0]->DecayBlob())
      AddDecayClustering(step2, daughters[0]->DecayBlob(), imax, idmother);
    if (daughters[1]->DecayBlob())
      AddDecayClustering(step2, daughters[1]->DecayBlob(), imax, idnew1);
    if (daughters[2]->DecayBlob())
      AddDecayClustering(step2, daughters[2]->DecayBlob(), imax, idnew2);
    ampl=step2;
  }
  else {
    PRINT_VAR(*blob);
    THROW(fatal_error, "1 -> n not implemented yet.");
  }
}

void Hard_Decay_Handler::AddPhotonsClustering(Cluster_Amplitude*& ampl,
                                              const Particle_Vector daughters,
                                              ParticlePair_Vector& photons,
                                              size_t& imax,
                                              const std::vector<size_t>& ids)
{
  DEBUG_FUNC(photons.size()<<" photons to be clustered");
  Particle * photon(photons.back().first);
  Particle * daughter(photons.back().second);
  photons.pop_back();
  size_t idmother(0);
  if      (daughter==daughters[0]) idmother=ids[0];
  else if (daughter==daughters[1]) idmother=ids[1];
  else if (daughter==daughters[2]) idmother=ids[2];
  else THROW(fatal_error,"Did not find id for "+daughter->Flav().IDName());
  msg_Debugging()<<"Cluster "<<photon->Flav()<<" "<<photon->Momentum()
                <<" with "<<daughter->Flav()<<" "<<ID(idmother)<<std::endl;
  Cluster_Amplitude* copy=ampl->InitPrev();
  copy->CopyFrom(ampl);
  copy->SetMS(ampl->MS());
  Cluster_Leg *lij(ampl->IdLeg(idmother));
  for (size_t i=0; i<ampl->Legs().size(); ++i)
    ampl->Leg(i)->SetStat(ampl->Leg(i)->Stat()|1);
  lij->SetStat(1|2|4);
  size_t idk(0);
  for (size_t i=0; i<copy->Legs().size(); ++i) {
    copy->Leg(i)->SetK(0);
    if (copy->Leg(i)->Id()!=idmother &&
        (copy->Leg(i)->Col().m_i==lij->Col().m_j ||
         copy->Leg(i)->Col().m_j==lij->Col().m_i))
      idk=copy->Leg(i)->Id();
  }
  if (lij->Col().m_i==0 && lij->Col().m_j==0) {
    // Ad hoc QED partner, must not be another soft photon
    size_t ampl_nout=ampl->Legs().size()-ampl->NIn();
    if (ampl_nout==1) idk=ampl->Leg(0)->Id();
    else {
      size_t select=ampl->Legs().size();
      size_t nvalid=0;
      for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
        if (!(ampl->Leg(i)->Id()&idmother || i>ampl->Legs().size()-1 ||
              ampl->Leg(i)->Flav().Kfcode()==kf_photon)) {
          nvalid++;
        }
      }
      if (nvalid==0) select=0;
      else {
        do {
          select=ampl->NIn()+floor(ran->Get()*ampl_nout);
        } while (ampl->Leg(select)->Id()&idmother ||
                 select>ampl->Legs().size()-1 ||
                 ampl->Leg(select)->Flav().Kfcode()==kf_photon);
      }
      msg_Debugging()<<"choose ("<<ID(ampl->Leg(select)->Id())<<") "
                     <<ampl->Leg(select)->Flav()<<std::endl;
      idk=ampl->Leg(select)->Id();
    }
  }
  else THROW(fatal_error,"Adding QED to coloured particle.");
  if (idk==0) THROW(fatal_error,"Colour partner not found");
  lij->SetK(idk);
  Cluster_Leg *d1(copy->IdLeg(idmother));
  size_t stat1(0), stat2(0);
  d1->SetMom(RecombinedMomentum(daughter,photons,stat1));
  d1->SetStat(stat1);
  d1->SetFlav(daughter->Flav());
  copy->CreateLeg(photon->Momentum(),photon->RefFlav());
  size_t idnew=1<<(++imax);
  copy->Legs().back()->SetId(idnew);
  copy->Legs().back()->SetStat(stat2);
  Cluster_Amplitude::SetColours(ampl->IdLeg(idmother),
                                copy->IdLeg(idmother),
                                copy->Legs().back());
  DEBUG_VAR(*copy);
  Cluster_Amplitude* tmp=copy;
  while (tmp->Next()) {
    tmp=tmp->Next();
    for (size_t i=0; i<tmp->Legs().size(); ++i) {
      if (tmp->Leg(i)->Id()&idmother) {
        tmp->Leg(i)->SetId(tmp->Leg(i)->Id()|idnew);
      }
      if (tmp->Leg(i)->K()&idmother) {
        tmp->Leg(i)->SetK(tmp->Leg(i)->K()|idnew);
      }
    }
    DEBUG_VAR(*tmp);
  }
  ampl=copy;
}

void Hard_Decay_Handler::AssignPhotons(const Particle_Vector& daughters,
                                       ParticlePair_Vector& photons)
{
  // for every photon, find charged particle that's closest
  // ignore radiation off charged resonance for now
  if (photons.size()) {
    Particle_Vector cdaughters;
    for (size_t i(0);i<daughters.size();++i)
      if (daughters[i]->Flav().Charge()) cdaughters.push_back(daughters[i]);
    if (cdaughters.size()==1) {
      for (size_t i(0);i<photons.size();++i)
        photons[i].second=cdaughters[0];
    }
    else {
      Vec4D cmom(0.,0.,0.,0.);
      Vec4D_Vector cmoms;
      for (size_t i(0);i<cdaughters.size();++i) {
        cmoms.push_back(cdaughters[i]->Momentum());
        cmom+=cmoms[i];
      }
      Poincare ccms(cmom);
      for (size_t i(0);i<cdaughters.size();++i) ccms.Boost(cmoms[i]);
      for (size_t i(0);i<photons.size();++i){
        Vec4D pmom(photons[i].first->Momentum());
        ccms.Boost(pmom);
        size_t id(0);
        double dR(pmom.DR(cmoms[0]));
        for (size_t j(1);j<cmoms.size();++j) {
          double dRj(pmom.DR(cmoms[j]));
          if (dRj<dR) { id=j; dR=dRj; }
        }
        photons[i].second=cdaughters[id];
      }
    }
    for (size_t i(0);i<photons.size();++i) {
      if (photons[i].first==photons[i].second)
        THROW(fatal_error,"Photon has not been assigned.");
      msg_Debugging()<<photons[i].first->Flav()<<" "
                     <<photons[i].first->Momentum()
                     <<" assigned to "<<photons[i].second->Flav()<<std::endl;
    }
  }
}

Vec4D Hard_Decay_Handler::RecombinedMomentum(const Particle * daughter,
                                             const ParticlePair_Vector& photons,
                                             size_t& stat)
{
  Vec4D mom(0.,0.,0.,0.);
  for (size_t i(0);i<photons.size();++i) {
    if (photons[i].second==daughter) {
      mom+=photons[i].first->Momentum();
      stat|=2|4;
    }
  }
  msg_Debugging()<<daughter->Flav()<<": "<<mom<<" "<<stat<<std::endl;
  return mom+daughter->Momentum();
}


void Hard_Decay_Handler::ReadDecayTable(Flavour decayer)
{
  DEBUG_FUNC(decayer);
  if (!m_store_results) return;
  Data_Reader reader = Data_Reader("|",";","!");
  reader.SetAddCommandLine(false);
  reader.AddComment("#");
  reader.AddComment("//");
  reader.AddWordSeparator("\t");
  reader.SetInputPath(m_resultdir+"/Decays/");
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
  if (!m_store_results) return;
  
  Decay_Map::iterator dmit;
  for (dmit=p_decaymap->begin(); dmit!=p_decaymap->end(); ++dmit) {
    ofstream ostr((m_resultdir+"/Decays/"+dmit->first.ShellName()).c_str());
    ostr<<"# Decay table for "<<dmit->first<<endl<<endl;
    Decay_Table::iterator dtit;
    for (dtit=dmit->second[0]->begin(); dtit!=dmit->second[0]->end(); ++dtit) {
      ostr<<(*dtit)->IDCode()<<"\t"<<(*dtit)->IWidth()<<"\t"
          <<(*dtit)->IDeltaWidth()<<"\t"<<(*dtit)->Max()<<endl;
    }
    ostr.close();
  }
}

bool Hard_Decay_Handler::Decays(const ATOOLS::Flavour& flav)
{
  if (flav.IsHadron()) return false;
  if (flav.Kfcode()==kf_tau && !m_decay_tau) return false;
  if (flav.IsStable()) return false;
  return true;
}

double Hard_Decay_Handler::Mass(const ATOOLS::Flavour &fl) const
{
  if (m_usemass==0) return fl.Mass();
  if (m_decmass.find(fl)!=m_decmass.end()) return fl.Mass(true);
  return fl.Mass();
}

