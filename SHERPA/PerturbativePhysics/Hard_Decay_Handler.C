#include"SHERPA/PerturbativePhysics/Hard_Decay_Handler.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
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
#include "METOOLS/SpinCorrelations/Polarized_CrossSections_Handler.H"
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

class ParticlePairFirstEnergySort {
public:
  bool operator()(const ParticlePair& a,const ParticlePair& b)
  { return (a.first->Momentum()[0]<b.first->Momentum()[0]); }
};

Hard_Decay_Handler::Hard_Decay_Handler() :
  p_newsublist(NULL), m_resultdir(""), m_offshell(""),
  m_decay_tau(false), m_set_widths(false), m_use_ho_sm_widths(true),
  m_br_weights(true), m_usemass(true), m_min_prop_width(0.0)
{
  auto& s = Settings::GetMainSettings();
  auto ds = s["HARD_DECAYS"];
  m_mass_smearing   = ds["Mass_Smearing"].SetDefault(1).Get<int>();
  /*
    TODO: Writing out a 1->3 channel which might have a large width in one
    resolved configuration, and a small one in another?
  */
  m_store_results   = ds["Store_Results"].SetDefault(0).Get<int>();
  m_br_weights      = ds["Apply_Branching_Ratios"].SetDefault(true).Get<bool>();
  m_decay_tau       = ds["Decay_Tau"].SetDefault(false).Get<bool>();
  m_set_widths      = ds["Set_Widths"].SetDefault(false).Get<bool>();
  m_use_ho_sm_widths= ds["Use_HO_SM_Widths"].SetDefault(true).Get<bool>();
  if (m_use_ho_sm_widths &&
      (Flavour(kf_h0).Mass()<125.07 || Flavour(kf_h0).Mass()>125.11 ||
       Flavour(kf_t).Mass()<172.4 || Flavour(kf_t).Mass()>172.6))
    THROW(fatal_error, "Use_HO_SM_Widths specified, but particle masses not in SM range.");
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

  m_spincorr=rpa->gen.HardSC();
  m_polcrosssec = ds["Pol_Cross_Section"]["Enabled"].SetDefault(false).Get<bool>();
  // if also final state polarization should be enabled, initialize polarization handler in a separate method in the
  // Initialization_Handler analogous to the Hard_Decay_Handler
  if (m_polcrosssec){
    if (!m_spincorr){
      THROW(fatal_error, "Calculation of polarized cross sections only possible together with spin correlations")
    }
    p_polarization_handler = new METOOLS::Polarized_CrossSections_Handler();
  }


  // also need to tell shower whats massive now
  // TODO: need to use the same mass-selector
  // for now, implement alike
  SetDecayMasses();
  m_qedmode = ds["QED_Corrections"].SetDefault(1).Get<size_t>();
  if (m_qedmode && m_spincorr) m_qedmode=2;
  if (m_qedmode>0 && !m_usemass) {
    THROW(fatal_error,std::string("QED corrections to hard decays only ")
                      +std::string("available in massive mode."));
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
  msg_Debugging()<<"Initialising hard decay tables: two-body decays.\n";
  for (dmit=p_decaymap->begin(); dmit!=p_decaymap->end(); ++dmit) {
    InitializeDirectDecays(dmit->second.at(0));
  }
  msg_Debugging()<<"Initialising hard decay tables: three-body decays.\n";
  for (dmit=p_decaymap->begin(); dmit!=p_decaymap->end(); ++dmit) {
    InitializeOffshellDecays(dmit->second.at(0));
  }
  msg_Debugging()<<"Initialising hard decay tables: customizing decay tables.\n";
  ApplySMWidths();
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


void Hard_Decay_Handler::ApplySMWidths()
{
  if (!m_use_ho_sm_widths) return;

  std::map<std::string, double> ho_sm_widths;
  // Higgs WG BRs 2022
  ho_sm_widths["25,5,-5"]    = 2.382E-03;
  ho_sm_widths["25,15,-15"]  = 2.565E-04;
  ho_sm_widths["25,13,-13"]  = 8.901E-07;
  ho_sm_widths["25,4,-4"]    = 1.182E-04;
  ho_sm_widths["25,3,-3"]    = 1E-06;
  ho_sm_widths["25,21,21"]   = 3.354E-04;
  ho_sm_widths["25,22,22"]   = 9.307E-06;
  ho_sm_widths["25,23,22"]   = 6.318E-06;
  ho_sm_widths["24,2,-1"]    = 0.7041;
  ho_sm_widths["24,4,-3"]    = 0.7041;
  ho_sm_widths["24,12,-11"]  = 0.2256;
  ho_sm_widths["24,14,-13"]  = 0.2256;
  ho_sm_widths["24,16,-15"]  = 0.2256;
  ho_sm_widths["-24,-2,1"]   = 0.7041;
  ho_sm_widths["-24,-4,3"]   = 0.7041;
  ho_sm_widths["-24,-12,11"] = 0.2256;
  ho_sm_widths["-24,-14,13"] = 0.2256;
  ho_sm_widths["-24,-16,15"] = 0.2256;
  ho_sm_widths["23,1,-1"]    = 0.3828;
  ho_sm_widths["23,2,-2"]    = 0.2980;
  ho_sm_widths["23,3,-3"]    = 0.3828;
  ho_sm_widths["23,4,-4"]    = 0.2980;
  ho_sm_widths["23,5,-5"]    = 0.3828;
  ho_sm_widths["23,11,-11"]  = 0.0840;
  ho_sm_widths["23,12,-12"]  = 0.1663;
  ho_sm_widths["23,13,-13"]  = 0.0840;
  ho_sm_widths["23,14,-14"]  = 0.1663;
  ho_sm_widths["23,15,-15"]  = 0.0840;
  ho_sm_widths["23,16,-16"]  = 0.1663;
  ho_sm_widths["6,24,5"]     = 1.32;
  ho_sm_widths["-6,-24,-5"]  = 1.32;

  for (const auto& channel_width : ho_sm_widths) {
    std::string channelid = channel_width.first;
    double width = channel_width.second;
    pair<Decay_Table*, Decay_Channel*> match=p_decaymap->FindDecayChannel(channelid, true);
    Decay_Channel* dc = match.second;
    if (!match.first || !dc) {
      PRINT_INFO("Ignoring unknown decay channel: " << channelid);
      continue;
    }
    dc->SetWidth(width);
  }
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
  Blob_Data_Base * wgtinfo((*blob)["MEWeightInfo"]);
  if (wgtinfo) *wgtinfo->Get<ME_Weight_Info*>()*=brfactor;

  Blob_Data_Base * wgtmap_bdb((*blob)["WeightsMap"]);
  if (wgtmap_bdb) wgtmap_bdb->Get<Weights_Map>()["BR"]=brfactor;

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
      (*p_newsublist)[i]->m_results["BR"]*=brfactor;
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
    copy->SetNLO(0);
    copy->SetFlag(1);
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
        size_t select(0);
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
    d1->SetFromDec(true);
    copy->CreateLeg(RecombinedMomentum(daughters[1],photons,stat2),
                    daughters[1]->RefFlav());
    copy->Legs().back()->SetFromDec(true);
    size_t idnew=1<<(++imax);
    copy->Legs().back()->SetId(idnew);
    copy->Legs().back()->SetStat(stat2);
    Cluster_Amplitude::SetColours(ampl->IdLeg(idmother),
                                  copy->IdLeg(idmother),
                                  copy->Legs().back());
    copy->SetIdNew(idnew);
    DEBUG_VAR(*copy);
    Cluster_Amplitude* tmp=copy;
    while (tmp->Next()) {
      tmp=tmp->Next();
      if (tmp->IdNew()&idmother) tmp->SetIdNew(tmp->IdNew()|idnew);
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
    step1->SetNLO(0);
    step1->SetFlag(1);
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
        size_t select(0);
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
    d1->SetFromDec(true);
    // todo: 1->2 qcd shower with ew fs recoil partner
    // d1->SetK(idmother);// not that simple: w->qq' has color connection in fs
    Decay_Channel* dc(NULL);
    Blob_Data_Base* data = (*blob)["dc"];
    if (data) {
      dc=data->Get<Decay_Channel*>();
      DEBUG_VAR(*dc);
    }
    else THROW(fatal_error, "Internal error.");
    Comix1to3* amp=dynamic_cast<Comix1to3*>(dc->GetDiagrams()[0]);
    if (!amp) THROW(fatal_error, "Internal error.");
    Flavour prop_flav=amp->Prop();
    Vec4D momd2=RecombinedMomentum(daughters[1],photons,stat2);
    Vec4D momd3=RecombinedMomentum(daughters[2],photons,stat3);
    Vec4D prop_mom=momd2+momd3;
    step1->CreateLeg(prop_mom, prop_flav);
    size_t idnew1=1<<(++imax);
    step1->Legs().back()->SetId(idnew1);
    step1->Legs().back()->SetStat(0);
    step1->Legs().back()->SetFromDec(true);
    Cluster_Amplitude::SetColours(ampl->IdLeg(idmother),
                                  step1->IdLeg(idmother),
                                  step1->Legs().back());
    step1->SetIdNew(idnew1);
    DEBUG_VAR(*step1);
    Cluster_Amplitude* tmp=step1;
    while (tmp->Next()) {
      tmp=tmp->Next();
      if (tmp->IdNew()&idmother) tmp->SetIdNew(tmp->IdNew()|idnew1);
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
    step2->SetNLO(0);
    step2->SetFlag(1);
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
    step2->Legs().back()->SetFromDec(true);
    Cluster_Amplitude::SetColours(step1->IdLeg(idnew1),
                                  step2->IdLeg(idnew1),
                                  step2->Legs().back());
    step2->SetIdNew(idnew2);
    DEBUG_VAR(*step2);
    tmp=step2;
    while (tmp->Next()) {
      tmp=tmp->Next();
      if (tmp->IdNew()&(idmother|idnew1))
	tmp->SetIdNew(tmp->IdNew()|idnew2);
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
  copy->SetNLO(0);
  copy->SetFlag(1);
  copy->SetMS(ampl->MS());
  Cluster_Leg *lij(ampl->IdLeg(idmother));
  for (size_t i=0; i<ampl->Legs().size(); ++i)
    ampl->Leg(i)->SetStat(ampl->Leg(i)->Stat()|1);
  lij->SetStat(1|2|4);
  size_t idk(0);
  for (size_t i=0; i<copy->Legs().size(); ++i) {
    copy->Leg(i)->SetK(0);
  }
  if (lij->Col().m_i!=0 || lij->Col().m_j!=0) {
    THROW(fatal_error,"Adding QED to coloured particle.");
  }
  // Ad hoc QED partner, must not be another soft photon
  size_t ampl_nout=ampl->Legs().size()-ampl->NIn();
  if (ampl_nout==1) idk=ampl->Leg(0)->Id();
  else {
    size_t select(0);
    size_t nvalid(0);
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
  copy->SetIdNew(idnew);
  DEBUG_VAR(*copy);
  Cluster_Amplitude* tmp=copy;
  while (tmp->Next()) {
    tmp=tmp->Next();
    if (tmp->IdNew()&idmother) tmp->SetIdNew(tmp->IdNew()|idnew);
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

bool Hard_Decay_Handler::Decays(const ATOOLS::Flavour& flav)
{
  if (flav.IsHadron()) return false;
  if (flav.Kfcode()==kf_tau && !m_decay_tau) return false;
  if (!flav.IsOn() || flav.IsStable()) return false;
  return true;
}

double Hard_Decay_Handler::Mass(const ATOOLS::Flavour &fl) const
{
  if (m_usemass==0) return fl.Mass();
  if (m_decmass.find(fl)!=m_decmass.end()) return fl.Mass(true);
  return fl.Mass();
}

