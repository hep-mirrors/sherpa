#include"SHERPA/PerturbativePhysics/Hard_Decay_Handler.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Org/Shell_Tools.H"
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
#include "EXTRA_XS/One2Three/H_to_bbg_Real.H"
#include "EXTRA_XS/One2Three/H_to_bb_Virtual.H"

#include <iostream>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <assert.h>
#include <limits>

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

class ParticlePairPairFirstEnergySort {
public:
  bool operator()(const ParticlePairPair& a,const ParticlePairPair& b)
  { return (a.first.first->Momentum()[0]+a.first.second->Momentum()[0]
            <b.first.first->Momentum()[0]+b.first.second->Momentum()[0]); }
};

Hard_Decay_Handler::Hard_Decay_Handler() :
/* Sets the variable values. It often checkes if there is a value specified, and if not, it sets a default value.
*/

  p_newsublist(NULL), m_resultdir(""), m_offshell(""),
  m_decay_tau(false), m_set_widths(false),
  m_br_weights(true), m_usemass(true), m_min_prop_width(0.0)   // Variables are set
{
  /* This code manages how particle decays are set up and processed in the simulation. In summary:

Setting Decay Masses:
The code first determines which particle flavours should be treated as massive by checking settings and the particles’ properties. 
It fills a set (m_decmass) with these flavours (and their antiparticles) and then logs the list of flavours considered “massive.”

Initializing Decay Tables and Channels:
For each flavour (and its antiparticle) that is allowed to decay (determined by the Decays() method), a new decay table (dt) is created, 
inserted into the global decay map (p_decaymap), and any pre-existing decay data is read from file via ReadDecayTable().

Then, the decay table entries are initialized in two stages:

Two-body decays (direct decays): Using InitializeDirectDecays(), the vertex information is retrieved from the model’s vertex table (via s_model->VertexTable()), 
any duplicate vertices are filtered out, and new decay channels are created for each valid vertex.

Three-body (offshell) decays: Using InitializeOffshellDecays(), additional channels arising from offshell intermediate states are set up and the 
overall channel statuses are updated.

Width Overwriting and Final Adjustments:
Finally, the code checks for user or heavy‐object–specific overrides for partial widths (through the HO SM widths settings) and sets the decay 
widths accordingly. It then writes the decay tables to files if the configuration requires storing results.

In essence, this code builds, initializes, and adjusts the decay channels for various particle flavours based on model input, 
user settings, and computed integration results, ensuring that the decays are treated consistently during the simulation.
  */
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
  m_min_prop_width  = ds["Min_Prop_Width"].SetDefault(std::numeric_limits<double>::min()).Get<double>();
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
  /* This code initializes the decay map by iterating over all particle flavours stored in s_kftable. 
  For each flavour, it creates a Flavour object from the KF code and then checks with the Decays() 
  method whether that flavour is eligible for decay. If it is, a new Decay_Table is created for the flavour, 
  inserted into the decay map (p_decaymap), and an attempt is made to read existing decay data from file via ReadDecayTable. 
  The same procedure is repeated for the corresponding antiparticle (flav.Bar()) if it is different from the particle itself and also allowed to decay.
  */
  p_decaymap = new Decay_Map(this);
  KFCode_ParticleInfo_Map::const_iterator it;
  for (it=s_kftable.begin();it!=s_kftable.end();it++) {
    Flavour flav(it->first);
    if (Decays(flav)) {
      Decay_Table* dt=new Decay_Table(flav, this);
      p_decaymap->insert(make_pair(flav,dt));
      ReadDecayTable(flav);    // this method essentially does this: m_read[decayer].insert(make_pair(decaychannel, results))
    }
    if (flav!=flav.Bar() && Decays(flav.Bar())) {
      Decay_Table* dt=new Decay_Table(flav.Bar(), this);
      p_decaymap->insert(make_pair(flav.Bar(),dt));
      ReadDecayTable(flav.Bar());
    }
  }
  
  // initialize them sorted by masses:
  /* This code iterates over all entries in the decay map (p_decaymap), which contains decay tables for different particle flavours. 
  For each decay table, it calls InitializeDirectDecays to set up two-body/ three-body decay channels, and it logs a debugging message indicating 
  that the initialization of hard decay tables (specifically two-body decays) is starting. */
  Decay_Map::iterator dmit;
  msg_Debugging()<<"Initialising hard decay tables: two-body decays.\n";
  for (dmit=p_decaymap->begin(); dmit!=p_decaymap->end(); ++dmit) {
    InitializeDirectDecays(dmit->second); // dmit->second retrieves the value of the p_decaymap that is currently being iterated over (flavour, decay table dt)
                                          // dt is a pointer to a Decay_Table object, which is essentially a container holding pointers to Decay_Channel objects. 
                                          // It behaves like a vector (or similar container) where you can iterate over its elements (each element being a Decay_Channel*)
  }
  msg_Debugging()<<"Initialising hard decay tables: three-body decays.\n";
  for (dmit=p_decaymap->begin(); dmit!=p_decaymap->end(); ++dmit) {
    InitializeOffshellDecays(dmit->second);
    dmit->second->UpdateChannelStatuses();
  }


  // overwrite partial widths from HO SM or user input
  if (ds["Use_HO_SM_Widths"].SetDefault(true).Get<bool>()) { // set default value true if not already defined
    if (Flavour(kf_h0).Mass()<125.07 || Flavour(kf_h0).Mass()>125.11 ||
        Flavour(kf_t).Mass()<172.4 || Flavour(kf_t).Mass()>172.6)  // check if the mass of the Higgs boson and top quark are in the expected SM range
      THROW(fatal_error, "Use_HO_SM_Widths specified, but particle masses not in SM range.");
    SetHOSMWidths(ds); // set default values
  }
  for (dmit=p_decaymap->begin(); dmit!=p_decaymap->end(); ++dmit) { // iterate over entries in decay map
    for (auto dc: *(dmit->second)) { // loop over decay channels
      auto s = ds["Channels"][dc->IDCode()]["Width"]; // setting for the width of the decay channel
      if (!s.HasDefault()) s.SetDefault(dc->IWidth()); // if no default value is set, use the integrated width value of the decay channel
      dc->SetWidth(s.Get<double>()); // assign value to the decay channel
    }
    dmit->second->UpdateWidth();
  }

  if (m_set_widths)
    for (dmit=p_decaymap->begin(); dmit!=p_decaymap->end(); ++dmit) {
      dmit->second->Flav().SetWidth(dmit->second->TotalWidth());
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
/*SetDecayMasses() retrieves the configuration for massive and massless partons from 
the global settings and checks for any inconsistencies between these lists. 
Based on the RESPECT_MASSIVE_FLAG, it then assigns default labels to partons, 
ensuring that partons explicitly marked as massive or massless are correctly handled. 
Finally, it populates the m_decmass member with the flavours (and their antiparticles) 
that are to be treated as massive and logs this configuration for further processing in decays.
*/
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
      size_t kf(allflavs[i].Kfcode()); // KF code: unique integer identifier for a particle flavor based on the PDG standard
      bool add(true);
      /*This loop iterates over all flavours in the vector "allflavs". For each flavour that is not marked as a dummy, 
      it retrieves its KF code and then checks whether that code already exists in either the "defpsmassive" or "defpsmassless" lists. 
      If the KF code is not present in either list (i.e. "add" remains true), it adds the code to "defpsmassive". Essentially, 
      this ensures that any flavour not explicitly designated as massless is by default treated as massive.
      */
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
  /*then remove and add those specified manually. This code synchronizes the default massive and massless lists with the manually specified ones. 
  For each element in the massive list (psmassive), it removes that element from the massless defaults (defpsmassless) and then ensures it is 
  present in the massive defaults (defpsmassive). Similarly, for each element in the massless list (psmassless), it removes that element from 
  the massive defaults and then adds it to the massless defaults if it isn’t already present. This ensures that the manually specified settings 
  override any previous assignments.
  */
  for (size_t i(0);i<psmassive.size();++i) {
    defpsmassless.erase(std::remove(defpsmassless.begin(),defpsmassless.end(),
                                    psmassive[i]),defpsmassless.end());
    if (std::find(defpsmassive.begin(),defpsmassive.end(),psmassive[i])==
        defpsmassive.end()) defpsmassive.push_back(psmassive[i]);
  }
  for (size_t i(0);i<psmassless.size();++i) {
    defpsmassive.erase(std::remove(defpsmassive.begin(),defpsmassive.end(),
                                   psmassless[i]),defpsmassive.end());
    // This code checks if the current element from the "psmassive" vector is already present in the "defpsmassive" list. 
    // The function std::find searches the "defpsmassive" container for the given element (psmassive[i]). If std::find 
    // returns defpsmassive.end(), it means the element was not found, so the element is added to "defpsmassive" via push_back.
    if (std::find(defpsmassless.begin(),defpsmassless.end(),psmassless[i])==
        defpsmassless.end()) defpsmassless.push_back(psmassless[i]);
  }
  // fill massive ones into m_decmass
  for (size_t i(0);i<defpsmassive.size();++i) {
    Flavour fl(defpsmassive[i],0);
    m_decmass.insert(fl);
    m_decmass.insert(fl.Bar());  // Antiparticle
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
  Vertex_Table::const_iterator vlit=s_model->VertexTable()->find(inflav);  // vertex list iterator; s_model is a global pointer (or reference) to the current 
                                                                          // model instance, typically of type Model_Base or a similar class. It contains all 
                                                                          // the model-specific information—such as the vertex table, the list of included particle 
                                                                          // flavours, and other interaction data
  // Check if the iterator "vlit" is valid (i.e. not equal to the end of the vertex table); 
  // if so, assign "vlit->second" (the vertex list associated with the found key) to "vlist,"
  const Vertex_List& vlist=(vlit!=s_model->VertexTable()->end()) ? (vlit->second) : Vertex_List(); 
  // temporary hack:
  // prepare reduced vertex list, not containing duplicates
  // in the sense of having identical external flavours
  Vertex_List reduced_vlist;
  for (Vertex_List::const_iterator vit=vlist.begin();vit!=vlist.end();++vit){  // iterate over the vertex list
    Vertex_List::const_iterator vjt=reduced_vlist.begin();
    for(;vjt!=reduced_vlist.end(); ++vjt)
      if((*vjt)->in == (*vit)->in) break;  // filtering out duplicates
    if(vjt==reduced_vlist.end()) reduced_vlist.push_back(*vit);
  }

  msg_Debugging()<<"Vertices:"<<std::endl;
  /* This loop processes each unique vertex in the reduced vertex list, and for every proper vertex it creates a new decay channel 
  for the flavor under consideration. It then adds the decay products (from the vertex), constructs a two-body decay diagram with a 
  Comix1to2 object, sets up the integration channels (using a Multi_Channel object and a Rambo phase-space generator), applies 
  configuration settings (such as the channel's active status), calculates the decay width, and if successful, adds the channel to 
  the decay table; otherwise, it deletes the channel.
  */
  for (size_t i=0;i<reduced_vlist.size();i++) { // iterate over the reduced vertex list
    Single_Vertex* sv=reduced_vlist[i];
    // Proper Vertex, if: 1. The vertex is not marked as already decayed (sv->dec is false).
    // 2. None of the incoming particles (legs) are dummy particles.
    // 3. The vertex has exactly three legs (necessary for direct decay).
    // 4. No external leg (except the first) has the same KF code as the first leg (avoiding self-couplings or unwanted configurations).
    if (!ProperVertex(sv)) continue; 
    msg_Debugging()<<"  "<<i<<": "<<*sv<<std::endl;
    Decay_Channel* dc=new Decay_Channel(inflav, this);
    for (int j=1; j<sv->NLegs(); ++j) dc->AddDecayProduct(sv->in[j]);

    if (dc->Flavs()[0].IDName() == "h0" && dc->Flavs()[1].IDName() == "b" && dc->Flavs()[2].IDName() == "bb"){
          std::cout << "flavs1[0].IDName(): " << dc->Flavs()[0].IDName() << "  to  " << dc->Flavs()[1].IDName() << dc->Flavs()[2].IDName() << std::endl;
        }

    Comix1to2* diagram=new Comix1to2(dc->Flavs());
    dc->AddDiagram(diagram);

    dc->SetChannels(new Multi_Channel(""));
    dc->Channels()->SetNin(1); // One incoming particle
    dc->Channels()->SetNout(dc->NOut()); // Number of outgoing particles

    // create a new instance of a Rambo phase-space generator. The constructor is given the number of incoming particles (1), the number of 
    // outgoing particles (dc->NOut()), and a pointer to the list of particle flavors (starting at &dc->Flavs().front()), along with the 
    //current Hard_Decay_Handler context. 
    Rambo* rambo = new Rambo(1,dc->NOut(),&dc->Flavs().front(),this); 
    dc->Channels()->Add(rambo);
    dc->Channels()->Reset();

    auto s = Settings::GetMainSettings()["HARD_DECAYS"]["Channels"][dc->IDCode()];
    dc->SetActive(s["Status"].SetDefault(dc->Active()).GetVector<int>());

    if (CalculateWidth(dc)) dt->AddDecayChannel(dc);
    else delete dc;
  }
  dt->UpdateWidth(); // recalculate sum of all individual decay channels' widths stored in dt
  if (m_set_widths) dt->Flav().SetWidth(dt->TotalWidth()); // set the particle’s width to newly computed total
}

void Hard_Decay_Handler::InitializeOffshellDecays(Decay_Table* dt) {
  DEBUG_FUNC(dt->Flav()<<" "<<dt->size());
  size_t dtsize=dt->size();
  for (size_t i=0;i<dtsize;++i) {
    Decay_Channel* dc=dt->at(i); // get the decay channel pointer stored at the i-th index
    // ResolveDecay attempts to generate additional (offshell or three-body) decay channel configurations that may arise 
    // from the original two-body decay represented by dc, and returns them in a vector of Decay_Channel pointers.
    vector<Decay_Channel*> new_dcs=ResolveDecay(dc); 
    if (TriggerOffshell(dc, new_dcs)) { // checks if the decay channel is offshell
      dc->SetActiveAll(-1);             // deactivate the original decay channel
      for (size_t j=0; j<new_dcs.size(); ++j) {
        // check for duplicates
        Decay_Channel* dup=dt->GetDecayChannel(new_dcs[j]->Flavs());
        if (dup && dup->Active(0)>=0) { // check if duplicate channel is already active
          DEBUG_INFO("Adding new diagram to "<<*dup);
          for (size_t k=0; k<new_dcs[j]->GetDiagrams().size(); ++k) {
            dup->AddDiagram(new_dcs[j]->GetDiagrams()[k]); // if so, add the new diagrams to the existing decay diagrams
          }
          for (size_t k=0; k<new_dcs[j]->Channels()->Channels().size(); ++k) { // add also new channels
            dup->AddChannel(new_dcs[j]->Channels()->Channels()[k]);
          }
          new_dcs[j]->ResetDiagrams(); // ResetDiagrams and ResetChannels are called to clear the diagrams and channels of the new decay channel
          new_dcs[j]->ResetChannels();
          delete new_dcs[j];
          new_dcs[j]=NULL;
          CalculateWidth(dup);
        }
        else { 
          DEBUG_INFO("Adding "<<new_dcs[j]->Name());
          auto s = Settings::GetMainSettings()["HARD_DECAYS"]["Channels"][new_dcs[j]->IDCode()];
          new_dcs[j]->SetActive(s["Status"].SetDefault(new_dcs[j]->Active()).GetVector<int>());
          dt->AddDecayChannel(new_dcs[j]);
        }
      }
    }
    else { // offshell decay not triggered
      DEBUG_INFO("Keeping factorised.");
      for (size_t j=0; j<new_dcs.size(); ++j) {
        if (new_dcs[j]) new_dcs[j]->SetActiveAll(-1); // mark all new decay channels as inactive
      }
    }
  }
  dt->UpdateWidth();
  if (m_set_widths) dt->Flav().SetWidth(dt->TotalWidth()); // set width of particle to summed decay width of all channels
}


void Hard_Decay_Handler::SetHOSMWidths(ATOOLS::Scoped_Settings& s)
{
  // Higgs WG BRs 2022
  s["Channels"]["25,5,-5"]   ["Width"].SetDefault(2.382E-03);
  s["Channels"]["25,15,-15"] ["Width"].SetDefault(2.565E-04);
  s["Channels"]["25,13,-13"] ["Width"].SetDefault(8.901E-07);
  s["Channels"]["25,4,-4"]   ["Width"].SetDefault(1.182E-04);
  s["Channels"]["25,3,-3"]   ["Width"].SetDefault(1E-06);
  s["Channels"]["25,21,21"]  ["Width"].SetDefault(3.354E-04);
  s["Channels"]["25,22,22"]  ["Width"].SetDefault(9.307E-06);
  s["Channels"]["25,23,22"]  ["Width"].SetDefault(6.318E-06);
  s["Channels"]["24,2,-1"]   ["Width"].SetDefault(0.7041);
  s["Channels"]["24,4,-3"]   ["Width"].SetDefault(0.7041);
  s["Channels"]["24,12,-11"] ["Width"].SetDefault(0.2256);
  s["Channels"]["24,14,-13"] ["Width"].SetDefault(0.2256);
  s["Channels"]["24,16,-15"] ["Width"].SetDefault(0.2256);
  s["Channels"]["-24,-2,1"]  ["Width"].SetDefault(0.7041);
  s["Channels"]["-24,-4,3"]  ["Width"].SetDefault(0.7041);
  s["Channels"]["-24,-12,11"]["Width"].SetDefault(0.2256);
  s["Channels"]["-24,-14,13"]["Width"].SetDefault(0.2256);
  s["Channels"]["-24,-16,15"]["Width"].SetDefault(0.2256);
  s["Channels"]["23,1,-1"]   ["Width"].SetDefault(0.3828);
  s["Channels"]["23,2,-2"]   ["Width"].SetDefault(0.2980);
  s["Channels"]["23,3,-3"]   ["Width"].SetDefault(0.3828);
  s["Channels"]["23,4,-4"]   ["Width"].SetDefault(0.2980);
  s["Channels"]["23,5,-5"]   ["Width"].SetDefault(0.3828);
  s["Channels"]["23,11,-11"] ["Width"].SetDefault(0.0840);
  s["Channels"]["23,12,-12"] ["Width"].SetDefault(0.1663);
  s["Channels"]["23,13,-13"] ["Width"].SetDefault(0.0840);
  s["Channels"]["23,14,-14"] ["Width"].SetDefault(0.1663);
  s["Channels"]["23,15,-15"] ["Width"].SetDefault(0.0840);
  s["Channels"]["23,16,-16"] ["Width"].SetDefault(0.1663);
  s["Channels"]["6,24,5"]    ["Width"].SetDefault(1.32);
  s["Channels"]["-6,-24,-5"] ["Width"].SetDefault(1.32);
}

bool Hard_Decay_Handler::TriggerOffshell(Decay_Channel* dc, vector<Decay_Channel*> new_dcs) 
/* This function determines whether additional resolved offshell decay configurations should be applied to a given decay channel
*/
{
  DEBUG_FUNC(dc->Name()<<"... "<<new_dcs.size());

  if (m_offshell=="Threshold") { // sum the masses of the daughter particles and compare that to the mass of the decaying particle
    double outmass=0.0;
    for (size_t j=1; j<dc->Flavs().size(); ++j)
      outmass+=this->Mass(dc->Flavs()[j]);
    DEBUG_INFO(this->Mass(dc->Flavs()[0])<<" vs "<<outmass);
    return (this->Mass(dc->Flavs()[0])<outmass);
  }
  // sums the widths from all resolved decay channels (new_dcs) and returns true if their total exceeds the original channel’s width
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
/*The ResolveDecay method takes an existing two-body decay channel (dc1) and attempts to “resolve” it into additional 
offshell (or three-body) decay configurations. 
*/
{
  DEBUG_FUNC(dc1->Name());
  vector<Decay_Channel*> new_dcs;
  const std::vector<ATOOLS::Flavour> flavs1(dc1->Flavs());

  bool bbbar_channel = false; // flag to check if the decay channel is Higgs to b bbar

  // filter out the b bbar channel to manually add h0 to bbarg later
  if (flavs1[0].IDName() == "h0" && flavs1[1].IDName() == "b" && flavs1[2].IDName() == "bb"){
    bbbar_channel = true;
  }
  std::cout << "candidate: " << flavs1[0].IDName() << "  to  " << flavs1[1].IDName() << flavs1[2].IDName() << std::endl;
  if (flavs1[0].IDName() == "h0" && flavs1[1].IDName() == "Z" && flavs1[2].IDName() == "Z"){
    bbbar_channel = false;
  }


  for (size_t j=1;j<flavs1.size();++j) { // iterate over each daughter flavor (starting at index 1)
    bool ignore=false;
    if (flavs1[j].Width()<m_min_prop_width && !(bbbar_channel)) continue; // skip if width is below threshold
    for (size_t k=1; k<j; ++k) { // skip duplicates
      // TODO Do we really have to avoid double counting e.g. in h -> Z Z?
      // Further iterations: W+ -> b t -> b W b -> b b .. ?
      if (flavs1[j]==flavs1[k]) ignore=true;
    }
    if (ignore) continue;
    Vertex_Table::const_iterator it=s_model->VertexTable()->find(flavs1[j]); // create iterator which points at to the map entry of corresponding flavor
    const Vertex_List& vlist(it->second);  // retrieve vertex list for corresponding flavor
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
      if (!ProperVertex(sv) && !(bbbar_channel && (sv->in[2].IDName() == "G"))) continue; // let Higgs to bbbar + gluon pass
      // TODO so far special case 1->3 only
      Decay_Channel* dc=new Decay_Channel(flavs1[0], this);
      size_t nonprop(0), propi(0), propj(0); 
      dc->AddDecayProduct(flavs1[3-j]); // particle that does not decay   
      dc->AddDecayProduct(sv->in[1]);   // decay products of the decaying particle 
      dc->AddDecayProduct(sv->in[2]);
      DEBUG_FUNC("trying "<<dc->Name()); 
      // TODO what about W+ -> b t -> b W b -> b b ... two diagrams, factor 2, ...?
      // TODO what about identical particles like W' -> b t -> b W b ... two diagrams, factor 2, ...?
      // TODO what about W -> W gamma -> l v gamma?
      for (size_t l=1; l<4; ++l) { // iterate over the decay products 
        if (dc->Flavs()[l]==flavs1[3-j]) nonprop=l; // find the non-propagating particle
      }
      for (size_t l=1; l<4; ++l) { // first != nonprop leg assigned to propi, the other one to propj
        if (l!=nonprop && propi>0) propj=l;
        else if (l!=nonprop && propi==0) propi=l;
      }

      assert(dc1->GetDiagrams().size()==1); // assert that original two-body decay channel has only one diagram (to catch inconsistencies)
      DEBUG_VAR(dc->Flavs());
      DEBUG_VAR(flavs1[j]);

      Spin_Amplitudes* diagram = nullptr; // parent class for H_to_bbg_Real and Comix1to3
      if (bbbar_channel && (sv->in[2].IDName() == "G")) {
        diagram = new H_to_bbg_Real(dc->Flavs(),flavs1[1],nonprop, propi, propj);

        // here: second diagram needed for h0 -> b bbar g
        Spin_Amplitudes* diagram2 = nullptr;
        diagram2 = new H_to_bbg_Real(dc->Flavs(),flavs1[2],propj,propi,nonprop);
        dc->AddDiagram(diagram2);
      } else {
        diagram = new Comix1to3(dc->Flavs(),flavs1[j],
        nonprop, propi, propj);

        // test (delete this later):
        const std::vector<ATOOLS::Flavour> flavs1(dc->Flavs());
        std::cout << "flavs1[0].IDName(): " << flavs1[0].IDName() << "  to  " << flavs1[1].IDName() << flavs1[2].IDName() << flavs1[3].IDName() << std::endl;
        if (flavs1[0].IDName() == "h0" && flavs1[1].IDName() == "Z" && flavs1[2].IDName() == "b"){
          std::cout << "flavs1[0].IDName(): " << flavs1[0].IDName() << "  to  " << flavs1[1].IDName() << flavs1[2].IDName() << flavs1[3].IDName() << std::endl;
        }
      }

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
/* CalculateWidth first checks whether the decaying particle’s mass is kinematically allowed to decay into its outgoing 
particles (i.e. its mass exceeds the sum of the decay-product masses). If so, it either retrieves pre-computed width values 
from stored results or integrates to compute the decay width using specified accuracy and iteration parameters, and updates 
the decay channel's width values; if not, it disables the channel by setting its width values to zero.
*/
{
  // Integrate or use results read from decay table file
  double outmass=0.0;
  // This loop sums the masses of all decay products in the decay channel. It starts at index 1 (ignoring the initial 
  // decaying particle at index 0) and calls the handler’s Mass method for each daughter, storing the total in outmass. 
  // This sum is later used to check if the decaying particle is heavy enough to allow the decay kinematically.
  for (size_t l=1; l<dc->Flavs().size(); ++l)
    outmass+=this->Mass(dc->Flavs()[l]);
  if (this->Mass(dc->Flavs()[0])>outmass) { // check if decay is kinematically allowed
    DEBUG_FUNC("Starting calculation of "<<dc->Name()<<" width now.");
    if (m_store_results && m_read.find(dc->Flavs()[0])!=m_read.end()) { // check if results are already stored in m_read
      if (m_read[dc->Flavs()[0]].find(dc->IDCode())!=
          m_read[dc->Flavs()[0]].end()) { // check if current decay channel is in m_read
        const vector<double>& results(m_read[dc->Flavs()[0]][dc->IDCode()]); // get precomputed set of results
        dc->SetIWidth(results[0]); // set the width of the decay channel
        dc->SetIDeltaWidth(results[1]); // set the delta width of the decay channel
        dc->SetMax(results[2]);
      }
      else { // otherwise calculate the width
        msg_Tracking()<<"    Integrating "<<dc->Name()<<endl;
        dc->CalculateWidth(m_int_accuracy,
                           m_int_target_mode==0 ? dc->Flavs()[0].Width() : 0.0,
                           m_int_niter); // this method belongs to the Decay_Channel class
      }
    }
    else {
      msg_Tracking()<<"    Integrating "<<dc->Name()<<endl;
      dc->CalculateWidth(m_int_accuracy,
                         m_int_target_mode==0 ? dc->Flavs()[0].Width() : 0.0,
                         m_int_niter);
    }
  }
  else { // decay kinematically not allowed: disable decay channel
    dc->SetActiveAll(-1);
    dc->SetIWidth(0.0);
    dc->SetIDeltaWidth(0.0);
    dc->SetMax(0.0);
  }
  auto s = Settings::GetMainSettings()["HARD_DECAYS"]["Channels"][dc->IDCode()];
  dc->SetWidth(dc->IWidth());
  dc->SetDeltaWidth(dc->IDeltaWidth());
  return true; // width has been successfully calculated
}

bool Hard_Decay_Handler::ProperVertex(MODEL::Single_Vertex* sv)
/* The ProperVertex method checks whether a given vertex is valid for constructing a decay channel. It returns false if:
- the vertex has already been used (indicated by its dec flag)
- if any of its legs are dummy particles
- the vertex has exactly three legs and that none of the outgoing legs (index i) have the same KF code as the incoming leg (index 0) */
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
/* A decay blob is created for the incoming particle (inpart) if it has not already been assigned one.
*/
{
  DEBUG_FUNC(inpart->Flav());
  if(inpart->DecayBlob()) THROW(fatal_error,"Decay blob already exists.");
  if(!Decays(inpart->Flav())) return; // return if the particle does not decay
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
  blob->AddData("p_onshell",new Blob_Data<Vec4D>(inpart->Momentum())); // add on-shell momentum
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
/* BRFactor calculates an overall branching‐ratio factor for the blob’s decay chain. Here’s how it works:

- It initializes the factor (brfactor) to 1.0.
- Then, for each outgoing particle in the blob, it looks up a decay table using the particle’s reference flavor.
- If a decay table is found, it multiplies brfactor by the ratio of the particle’s active width (the width corresponding to the decay channels that are “active”) to the total width available for decays.
- If the particle itself decays (has a decay blob of type Hard_Decay), the function calls itself recursively on that blob and multiplies the result into brfactor.
- Finally, it returns the overall brfactor.
This factor is later used to correct event weights by accounting for the probabilities of the various decay channels.
*/
{
  double brfactor=1.0;
  for (size_t i=0; i<blob->NOutP(); ++i) {
    Particle* part=blob->OutParticle(i);
    Decay_Table* dt=p_decaymap->FindDecay(part->RefFlav());
    if (dt) {
      brfactor*=dt->ActiveWidth(0)/dt->TotalWidth();
      // If outgoing particle itself has a further decay (with type “Hard_Decay”), the factor from its decay blob is recursively multiplied in.
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
  /* Call the base‐class version of TreatInitialBlob (via Decay_Handler_Base::TreatInitialBlob). In that function, the following is done:
• The blob (representing a decaying system) is logged and a flag is reset (m_decaychainend).
• Any decay matrices from a previous event are cleared (if spin‐correlations are enabled).
• The list of outgoing (“daughter”) particles is obtained from the blob.
• A vector of indices is created and then randomly shuffled. This randomization helps to avoid bias when applying spin‐correlation effects.
• For each daughter, a check is performed to ensure that its momentum is consistent with its on‑shell mass; if not, the event is retried.
  */
  Decay_Handler_Base::TreatInitialBlob(blob, amps, origparts);

  double brfactor=m_br_weights ? BRFactor(blob) : 1.0;
  DEBUG_VAR(brfactor);
  Blob_Data_Base * bdbmeweight((*blob)["MEWeight"]); // retrieve the ME weight from the blob data
  if (bdbmeweight) {
    // msg_Out()<<METHOD<<"(ME = "<<bdbmeweight->Get<double>()<<", "
    // 	     <<"BR = "<<brfactor<<") for\n"<<"   "
    // 	     <<blob->InParticle(0)->Flav()<<" "
    // 	     <<blob->InParticle(1)->Flav()<<" -->";
    // for (size_t i=0;i<blob->NOutP();i++) 
    //   msg_Out()<<" "<<blob->OutParticle(i)->Flav();
    // msg_Out()<<".\n";
    bdbmeweight->Set<double>(brfactor*bdbmeweight->Get<double>()); // if data exists, multiply it by the branching ratio factor so that the 
                                                                  // weight reflects the decay’s branching ratio.
  }
  // update also MEWeightInfo and WeightsMap with the branching ratio factor
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

    if (p_newsublist) { // check, if the new sublist already exists
      // if so, clear the old one
      for (size_t i=0; i<p_newsublist->size(); ++i) delete (*p_newsublist)[i];
      p_newsublist->clear();
    }
    else p_newsublist=new NLO_subevtlist();
    for (size_t i=0; i<sublist->size(); ++i) {
      // iterate over sub events and replace decayed particles
      NLO_subevt* sub((*sublist)[i]);
      DEBUG_VAR(*(*sublist)[i]);

      if (sub->IsReal()) newn+=1; // "real": extra parton (vs. virtual)
      Flavour* newfls = new Flavour[newn];
      Vec4D* newmoms = new Vec4D[newn];
      size_t* newid = new size_t[newn];
      for (size_t n=0; n<newn; ++n) newid[n]=0;
      NLO_subevt* newsub=new NLO_subevt(*sub); // created new subevent as a copy of the current subevent
      newsub->m_n=newn; // number of legs
      newsub->p_id=newid;
      newsub->p_fl=newfls;
      newsub->p_mom=newmoms;
      newsub->m_delete=true;
      p_newsublist->push_back(newsub);

      int nin=blob->NInP();
      for (size_t j=0, jnew=0; j<nin+blob->NOutP()-1; ++j, ++jnew) { // iterate over all in- and out-particles
        if (j<2 || decayprods[j-nin].size()==1) { // If particle did not decay (or when the corresponding decay products list has only one particle):
                                                  // the code simply copies the original flavor and momentum from the subevent into the new arrays.
          newfls[jnew]=sub->p_fl[j];
          newmoms[jnew]=sub->p_mom[j];
          continue;
        }
        if (sub->p_fl[j]!=blob->OutParticle(j-nin)->Flav()) { // check that the flavor in the subevent agrees with the corresponding blob’s outgoing particle.
          THROW(fatal_error, "Internal Error 1");
        }

        Vec4D oldmom=blob->OutParticle(j-nin)->Momentum(); // get the momentum of the outgoing particle
        Vec4D newmom=sub->p_mom[j];                        // get the momentum of the corresponding subevent particle
        Poincare cms(oldmom);
        Poincare newframe(newmom);
        newframe.Invert();

        list<Particle*>::const_iterator it;
        for (it=decayprods[j-nin].begin(); it!=decayprods[j-nin].end(); ++it) { // iterate over a list of particles (the decay products for one blob particle)
          newfls[jnew]=(*it)->Flav(); // copy the decay product’s flavour into the new flavour array (newfls) at the index jnew
          newmoms[jnew]=Vec4D(newframe*(cms*(*it)->Momentum())); // calculate a new momentum for the decay product by applying a Lorentz transformation
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
  Particle_Vector daughters, splitphotonproducts;
  ParticlePair_Vector photons;
  ParticlePairPair_Vector splitphotons;
  for (size_t i(0);i<blob->GetOutParticles().size();++i) {
    Particle * p(blob->OutParticle(i));
    if      (p->Info()=='S') photons.push_back(make_pair(p,p));
    else if (p->Info()=='s') splitphotonproducts.push_back(p);
    else                     daughters.push_back(p);
  }
  msg_Debugging()<<"daughters: ";
  for (size_t i(0);i<daughters.size();++i)
    msg_Debugging()<<daughters[i]->Flav().IDName()<<" ";
  msg_Debugging()<<" +  "<<photons.size()<<" soft photon(s)"
                 <<" +  "<<splitphotonproducts.size()<<" photon splitting products"
                 <<std::endl;
  UnsplitPhotons(splitphotonproducts,splitphotons);
  std::sort(photons.begin(),photons.end(),ParticlePairFirstEnergySort());
  std::sort(splitphotons.begin(),splitphotons.end(),ParticlePairPairFirstEnergySort());
  AssignSplitPhotons(daughters,splitphotons);
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
    d1->SetMom(RecombinedMomentum(daughters[0],photons,splitphotons,stat1));
    d1->SetStat(stat1);
    d1->SetFlav(daughters[0]->Flav());
    d1->SetFromDec(true);
    copy->CreateLeg(RecombinedMomentum(daughters[1],photons,splitphotons,stat2),
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
    // always combine radiated (split) photons with
    // identified primary charged decay particle
    while (splitphotons.size())
      AddSplitPhotonsClustering(copy, daughters, splitphotons, imax, ids);
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
    d1->SetMom(RecombinedMomentum(daughters[0],photons,splitphotons,stat1));
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
    Vec4D momd2=RecombinedMomentum(daughters[1],photons,splitphotons,stat2);
    Vec4D momd3=RecombinedMomentum(daughters[2],photons,splitphotons,stat3);
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
    // always combine radiated (split) photons with
    // identified primary charged decay particle
    while (splitphotons.size())
      AddSplitPhotonsClustering(step2,daughters,splitphotons,imax,ids);
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

void Hard_Decay_Handler::AddSplitPhotonsClustering(Cluster_Amplitude*& ampl,
                                                   const Particle_Vector daughters,
                                                   ParticlePairPair_Vector& splitphotons,
                                                   size_t& imax,
                                                   const std::vector<size_t>& ids)
{
  DEBUG_FUNC(splitphotons.size()<<" split photons to be clustered");
  // will need to construct two cluster steps
  // 1) the (offshell) photon is radiated
  // 2) the (offshell) photon splits into the identified pair
  Particle * splitphoton1(splitphotons.back().first.first);
  Particle * splitphoton2(splitphotons.back().first.second);
  Particle * daughter(splitphotons.back().second);
  splitphotons.pop_back();
  size_t idmother(0),idphoton(0);
  if      (daughter==daughters[0]) idmother=ids[0];
  else if (daughter==daughters[1]) idmother=ids[1];
  else if (daughter==daughters[2]) idmother=ids[2];
  else THROW(fatal_error,"Did not find id for "+daughter->Flav().IDName());
  // construct recombined (off-shell) photon momentum
  Vec4D pmom(splitphoton1->Momentum()+splitphoton2->Momentum());
  msg_Debugging()<<"Cluster recombined photon with "<<pmom
                <<" with "<<daughter->Flav()<<" "<<ID(idmother)<<std::endl;

  // construct photon splitting cluster step
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
  for (size_t i=0; i<copy->Legs().size(); ++i)
    copy->Leg(i)->SetK(0);
  if (lij->Col().m_i!=0 || lij->Col().m_j!=0)
    THROW(fatal_error,"Adding QED to coloured particle.");
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
  d1->SetMom(RecombinedMomentum(daughter,ParticlePair_Vector(),splitphotons,stat1));
  d1->SetStat(stat1);
  d1->SetFlav(daughter->Flav());
  copy->CreateLeg(pmom,Flavour(kf_photon));
  idphoton=1<<(++imax);
  copy->Legs().back()->SetId(idphoton);
  copy->Legs().back()->SetStat(stat2);
  Cluster_Amplitude::SetColours(ampl->IdLeg(idmother),
                                copy->IdLeg(idmother),
                                copy->Legs().back());
  copy->SetIdNew(idphoton);
  DEBUG_VAR(*copy);
  // update IDs
  Cluster_Amplitude* tmp=copy;
  while (tmp->Next()) {
    tmp=tmp->Next();
    if (tmp->IdNew()&idmother) tmp->SetIdNew(tmp->IdNew()|idphoton);
    for (size_t i=0; i<tmp->Legs().size(); ++i) {
      if (tmp->Leg(i)->Id()&idmother) {
        tmp->Leg(i)->SetId(tmp->Leg(i)->Id()|idphoton);
      }
      if (tmp->Leg(i)->K()&idmother) {
        tmp->Leg(i)->SetK(tmp->Leg(i)->K()|idphoton);
      }
    }
    DEBUG_VAR(*tmp);
  }
  ampl=copy;

  // construct the photon splitting cluster step
  copy=ampl->InitPrev();
  copy->CopyFrom(ampl);
  copy->SetNLO(0);
  copy->SetFlag(1);
  copy->SetMS(ampl->MS());
  lij=ampl->IdLeg(idphoton);
  for (size_t i=0; i<ampl->Legs().size(); ++i)
    ampl->Leg(i)->SetStat(ampl->Leg(i)->Stat()|1);
  lij->SetStat(1|2|4);
  idk=0;
  for (size_t i=0; i<copy->Legs().size(); ++i)
    copy->Leg(i)->SetK(0);
  if (lij->Col().m_i!=0 || lij->Col().m_j!=0)
    THROW(fatal_error,"Adding QED to coloured particle.");
  // Ad hoc QED partner, must not be another soft photon
  ampl_nout=ampl->Legs().size()-ampl->NIn();
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
  d1=copy->IdLeg(idphoton);
  stat1=0; stat2=0;
  // we know the splitting products are not evolved further yet
  d1->SetMom(splitphoton1->Momentum());
  d1->SetStat(stat1);
  d1->SetFlav(splitphoton1->Flav());
  copy->CreateLeg(splitphoton2->Momentum(),splitphoton2->Flav());
  size_t idnew=1<<(++imax);
  copy->Legs().back()->SetId(idnew);
  copy->Legs().back()->SetStat(stat2);
  Cluster_Amplitude::SetColours(ampl->IdLeg(idmother),
                                copy->IdLeg(idmother),
                                copy->Legs().back());
  copy->SetIdNew(idnew);
  DEBUG_VAR(*copy);
  // update IDs
  tmp=copy;
  while (tmp->Next()) {
    tmp=tmp->Next();
    if (tmp->IdNew()&idphoton) tmp->SetIdNew(tmp->IdNew()|idnew);
    for (size_t i=0; i<tmp->Legs().size(); ++i) {
      if (tmp->Leg(i)->Id()&idphoton) {
        tmp->Leg(i)->SetId(tmp->Leg(i)->Id()|idnew);
      }
      if (tmp->Leg(i)->K()&idphoton) {
        tmp->Leg(i)->SetK(tmp->Leg(i)->K()|idnew);
      }
    }
    DEBUG_VAR(*tmp);
  }
  ampl=copy;
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
  d1->SetMom(RecombinedMomentum(daughter,photons,ParticlePairPair_Vector(),stat1));
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
  // update IDs
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

void Hard_Decay_Handler::UnsplitPhotons(const ATOOLS::Particle_Vector& spp,
                                        ParticlePairPair_Vector& splitphotons)
{
  // undo the splitting of a photon, recombine particle-anti-particle pair
  // into (off-shell) photon
  if (spp.size()) {
    if (spp.size()%2!=0)
      THROW(fatal_error,"Uneven number of photon-splitting products, cannot proceed.");
    if (spp.size()==2) {
      splitphotons.push_back(make_pair(make_pair(spp[0],spp[1]),spp[0]));
    }
    else {
      // find (particle,anti-particle)-pairs in spp,
      // order by invariant mass, starting with the smallest
      std::map<double,std::pair<size_t,size_t> > pairtable;
      for (size_t i(0);i<spp.size();++i) {
        for (size_t j(i+1);j<spp.size();++j) {
          if (spp[i]->Flav()==spp[j]->Flav().Bar()) {
            double m2((spp[i]->Momentum()+spp[j]->Momentum()).Mass());
            pairtable[m2]=make_pair(i,j);
          }
        }
      }
      if (pairtable.empty()) {
        msg_Debugging()<<"no pairs found"<<std::endl;
      }
      else if (pairtable.size()==1) {
        splitphotons.push_back(make_pair(make_pair(spp[0],spp[1]),spp[0]));
      }
      else {
        if (msg_LevelIsDebugging()) {
          msg_Debugging()<<"pairs found:\n";
          for (std::map<double,std::pair<size_t,size_t> >::const_iterator
               it=pairtable.begin();it!=pairtable.end();++it)
            msg_Debugging()<<it->second.first<<" "<<it->second.second<<": "
                <<", m2="<<it->first<<std::endl;
        }
        std::vector<size_t> usedindices;
        for (std::map<double,std::pair<size_t,size_t> >::const_iterator it=pairtable.begin();
             it!=pairtable.end();++it) {
          bool valid(true);
          for (size_t i(0);i<usedindices.size();++i)
            if (it->second.first==usedindices[i] ||
                it->second.second==usedindices[i]) { valid=false; break; }
          if (!valid) continue;
          usedindices.push_back(it->second.first);
          usedindices.push_back(it->second.second);
          msg_Debugging()<<"constructing split pair: P -> "
                        <<spp[it->second.first]->Flav()<<" "
                        <<spp[it->second.second]->Flav()<<" ,  m2 = "
                        <<it->first<<std::endl;
          splitphotons.push_back(make_pair(make_pair(spp[it->second.first],
                                                     spp[it->second.second]),
                                           spp[it->second.first]));
        }
        if (2*splitphotons.size()!=spp.size()) {
          msg_Error()<<METHOD<<"(): Found "<<splitphotons.size()<<" pairs in "
                     <<spp.size()<<" particles."<<std::endl;
          THROW(fatal_error,"Wrong number of pairs found.");
        }
      }
    }
  }
}

void Hard_Decay_Handler::AssignPhotons(const Particle_Vector& daughters,
                                       ParticlePair_Vector& photons)
{
  // for every photon, find charged particle that's closest (using dR)
  // ignore radiation off charged resonance for now
  if (photons.size()) {
    // first
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

void Hard_Decay_Handler::AssignSplitPhotons(const Particle_Vector& daughters,
                                            ParticlePairPair_Vector& splitphotons)
{
  // for every split photon, find charged particle that's closest (using dR)
  // ignore radiation off charged resonance for now
  if (splitphotons.size()) {
    // first
    Particle_Vector cdaughters;
    for (size_t i(0);i<daughters.size();++i)
      if (daughters[i]->Flav().Charge()) cdaughters.push_back(daughters[i]);
    if (cdaughters.size()==1) {
      for (size_t i(0);i<splitphotons.size();++i)
        splitphotons[i].second=cdaughters[0];
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
      for (size_t i(0);i<splitphotons.size();++i){
        Vec4D pmom(splitphotons[i].first.first->Momentum()
                   +splitphotons[i].first.second->Momentum());
        ccms.Boost(pmom);
        size_t id(0);
        double dR(pmom.DR(cmoms[0]));
        for (size_t j(1);j<cmoms.size();++j) {
          double dRj(pmom.DR(cmoms[j]));
          if (dRj<dR) { id=j; dR=dRj; }
        }
        splitphotons[i].second=cdaughters[id];
      }
    }
    for (size_t i(0);i<splitphotons.size();++i) {
      if (splitphotons[i].first.first==splitphotons[i].second ||
          splitphotons[i].first.second==splitphotons[i].second)
        THROW(fatal_error,"Split photon has not been assigned.");
      msg_Debugging()<<splitphotons[i].first.first->Flav()<<" "
                     <<splitphotons[i].first.first->Momentum()<<" and "
                     <<splitphotons[i].first.second->Flav()<<" "
                     <<splitphotons[i].first.second->Momentum()
                     <<" assigned to "<<splitphotons[i].second->Flav()<<std::endl;
    }
  }
}

Vec4D Hard_Decay_Handler::RecombinedMomentum(const Particle * daughter,
                                             const ParticlePair_Vector& photons,
                                             const ParticlePairPair_Vector& splitphotons,
                                             size_t& stat)
{
  Vec4D mom(0.,0.,0.,0.);
  for (size_t i(0);i<splitphotons.size();++i) {
    if (splitphotons[i].second==daughter) {
      mom+=splitphotons[i].first.first->Momentum()
           +splitphotons[i].first.second->Momentum();
      stat|=2|4;
    }
  }
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
/* Set Up Reader:
A Data_Reader object is created and configured with custom delimiters ("|", ";", "!"), comment markers (lines beginning with "#" or "//"), 
and word separators ("\t"). The reader’s input path is set to m_resultdir, and the input file is chosen based on decayer.ShellName().

Read and Process File Content:
The reader attempts to load the file contents into a matrix (a vector of string vectors). For each line that has exactly 4 tokens, 
the first token is treated as the decay channel identifier, and the next three tokens are converted to double values (using ToType<double>), 
resulting in a vector of three numbers (e.g., widths or related results). These pairs (decay channel and corresponding results) 
are then stored in the m_read map using the decayer as a key.
*/
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
        m_read[decayer].insert(make_pair(decaychannel, results));     // stores decay information
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
    for (dtit=dmit->second->begin(); dtit!=dmit->second->end(); ++dtit) {
      ostr<<setw(25)<<left<<(*dtit)->IDCode()<<"\t"
          <<setw(12)<<left<<(*dtit)->IWidth()<<"\t"
          <<setw(12)<<left<<(*dtit)->IDeltaWidth()<<"\t"
          <<setw(12)<<left<<(*dtit)->Max()<<endl;
    }
    ostr.close();
  }
}

bool Hard_Decay_Handler::Decays(const ATOOLS::Flavour& flav)
/* This method determines whether a given particle flavour should be allowed to decay. It checks three conditions:

If the flavour represents a hadron, it returns false.
If the particle is a tau (identified by its KF code) and decays for tau are disabled (m_decay_tau is false), it returns false.
If the particle is either marked as off or is stable, it also returns false.
If none of these conditions hold, the method returns true, meaning the particle is eligible for decay.
*/
{
  if (flav.IsHadron()) return false;
  if (flav.Kfcode()==kf_tau && !m_decay_tau) return false;
  if (!flav.IsOn() || flav.IsStable()) return false;
  return true;
}

double Hard_Decay_Handler::Mass(const ATOOLS::Flavour &fl) const
/* m_usemass is a flag (an instance variable) that determines whether the hard decay handler should 
consider a particle’s mass when computing decay kinematics and widths. m_decmass is a set (specifically a Flavour_Set) 
containing the particle flavours that are designated as “massive” for the purposes of decay; if a flavour is present in m_decmass, 
then its “massive” value (obtained by calling Mass(true)) is used, otherwise the default mass (or zero) is returned.
*/
{
  if (m_usemass==0) return fl.Mass();
  if (m_decmass.find(fl)!=m_decmass.end()) return fl.Mass(true);
  return fl.Mass();
}

