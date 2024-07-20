#include "ATOOLS/Phys/Fragmentation_Base.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Message.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Phys/KF_Table.H"
#include<set>
#include "Pythia8/Pythia.h"

using namespace ATOOLS;
using namespace std;

namespace SHERPA {

class Pythia8_Hadronisation : public Fragmentation_Base {

public:
  Pythia8_Hadronisation(const string& shower)
  {
    PRINT_INFO("Initialising Pythia8 hadronisation interface");
    rpa->gen.AddCitation(
        1, "The Pythia8 library is described in \\cite{Bierlich:2022pfr}.");

    // Initialise Pythia object
    m_pythia.readString("ProcessLevel:all = off");

    // Optionally switch off resonance decays, or only showers in them.
    m_pythia.readString("ProcessLevel:resonanceDecays = off");
    m_pythia.readString("PartonLevel:FSRinResonances = off");

    // Switch off automatic event listing in favour of manual.
    m_pythia.readString("Next:numberCount = 0");
    m_pythia.readString("Next:numberShowLHA = 0");
    m_pythia.readString("Next:numberShowInfo = 0");
    m_pythia.readString("Next:numberShowProcess = 0");
    m_pythia.readString("Next:numberShowEvent = 0");
    m_pythia.readString("Main:timesAllowErrors = 500");
    m_pythia.readString("Check:mTolErr = 1e-1");
    const double maxproperlifetime{m_settings["HADRON_DECAYS"]["Max_Proper_Lifetime"].Get<double>()};
    if (maxproperlifetime > 0.0) {
      m_pythia.readString("ParticleDecays:limitTau0 = on");
      m_pythia.readString("ParticleDecays:tau0Max = "+std::to_string(maxproperlifetime));
    }

    AssignDecays();
    ApplyPythiaSettings();
    HarmonizeMasses();

    m_pythia.init();
  }

  ~Pythia8_Hadronisation()
  {
    m_pythia.stat();
  }

  Return_Value::code Hadronize(Blob_List * blobs)
  {
    Pythia8::Event& event      = m_pythia.event;
    for (Blob_List::iterator blit=blobs->begin();blit!=blobs->end();++blit) {
      if ((*blit)->Has(blob_status::needs_hadronization) &&
	  (*blit)->Type()==btp::Fragmentation) {
	Blob * blob = (*blit);
	blob->SetTypeSpec("Pythia8");
	Sherpa2Pythia(blob, event);
	///Hadronization step
	if (!m_pythia.next()) {
          Blob * showerblob(blob->InParticle(0)->ProductionBlob());
          Blob * decblob(showerblob->InParticle(0)->ProductionBlob());
          if (decblob->Type()!=btp::Hadron_Decay) {
            event.list();
            msg_Error()<<"Pythia8 hadronisation failed.\n"<<endl;
            return Return_Value::Error;
          }
          else {
            msg_Tracking()<<"Error in "<<METHOD<<"."<<endl
                          <<"   Hadronization of partonic decay failed. Retry the event."<<endl;
            return Return_Value::Retry_Event;
          }
	}
	FillFragmentationBlob(blobs, blob, event);
        if (m_pythiadecays) {
          break;
        }
      }
    }
    if (m_shrink) Shrink(blobs);
    return Return_Value::Success;
  }

private:
  void Sherpa2Pythia(Blob * blob, Pythia8::Event& pevt)
  {
    /*
    pevt.append( id, status, col, acol, p, m)
    pevt.append( id, status, col, acol, px, py, pz, e, m)
    pevt.append( id, status, mother1, mother2, daughter1, daughter2, col, acol, p, m)
    pevt.append( id, status, mother1, mother2, daughter1, daughter2, col, acol, px, py, pz, e, m)
    The PDG particle code id and the Les Houches Accord colour col and anticolour acol tags must be set correctly.
    The four-momentum and mass have to be provided in units of GeV; if you omit the mass it defaults to 0.
    Outgoing particles that should hadronize should be given status code 23
    For normal hadronization/decays in pythia.next() the history encoded in the mother and daughter indices is not used.
    Therefore the first two append methods, which set all these indices vanishing, should suffice.
    The subsequent hadronization/decays will still be properly documented.
    The exception is when you want to include junctions in your string topology, i.e. have three
    string pieces meet.  Then you must insert in your event record the (decayed) particle that is
    the reason for the presence of a junction, e.g. a baryon beam remnant from which several valence
    quarks have been kicked out, or a neutralino that underwent a baryon-number-violating decay.
    This particle must have as daughters the three partons that together carry the baryon number.
    */
    // Reset event record to allow for new event.
    pevt.reset();
    int id(0);
    int status(23);
    int col(101);
    int acol(0);
    double px(0);
    double py(0);
    double pz(0);
    double e(0);
    double m(0);
    for (int i(0);i<blob->NInP();++i) {
      Particle * part = blob->InParticle(i);
      id = int(part->Flav());
      status = 23;
      col = part->GetFlow(1);
      acol = part->GetFlow(2);
      px = part->Momentum()[1];
      py = part->Momentum()[2];
      pz = part->Momentum()[3];
      e = part->Momentum()[0];
      m = part->FinalMass();
      pevt.append( IdToPythia(id), status, col, acol, px, py, pz, e, m);
      pevt[i].vProd(part->XProd()[1],part->XProd()[2],part->XProd()[3],part->XProd()[0]);
    }
  }

  void FillFragmentationBlob(Blob_List * bloblist, Blob * blob, Pythia8::Event& pevt)
  {
    /*
      Go through now hadronized Pythia event and fill sherpa fragmentation blob with particles
      resulting from hadronization.  If Pythia already did hadron decays the proper sherpa decay
      blobs should also be created.
    */
    m_processed.clear();

    for (int i = 1; i < pevt.size(); ++i) {
      // Status -23 are incoming partons.
      // Loop over them and iteratively check their daughters.
      // Only handle particles that have not been processed already.
      auto find = m_processed.find(i);
      if (pevt[i].status() == -23 && find == m_processed.end()) {
        m_processed.insert(i);
        HandleDaughters(bloblist, blob, pevt, i);
      }
    }
    if (!m_pythiadecays){
      blob->SetStatus(blob_status::needs_hadrondecays);
    }
  }

  void HandleDaughters(Blob_List * bloblist, Blob * decayblob, Pythia8::Event& pevt, int i){
    /*
      Loop over the daughters of the particle at position 'i' in the Pythia event.
      Either add them to outgoing particles directly or continue with their daughters.
      If the particle is unstable and Pythia handled the decays also initialize that particles's decay,
    */
    int begin;
    int end;

    // No daughters (unexpected)
    if (pevt[i].daughter1() == 0 && pevt[i].daughter2() == 0){
      THROW(fatal_error,"Particle does not have any daughters to handle.");
    }
    // One carbon copy as sole daughter. Example are recoil effects in the shower or oscillations.
    else if (pevt[i].daughter1() == pevt[i].daughter2() && pevt[i].daughter1() > 0){
      begin = pevt[i].daughter1();
      end = pevt[i].daughter1()+1;
    }
    // Only one (non-copy) daughter exists. Examples are beams, 2->1 hard interactions and clustering of nearby partons.
    else if (pevt[i].daughter1() > 0 && pevt[i].daughter2()== 0){
      begin = pevt[i].daughter1();
      end = pevt[i].daughter1()+1;
    }
    // Normal decay with daughters ranging from daughter1 to daughter2. (Daughter1 at index 5 and Daughter2 at index 8 means the loop works through 5,6,7,8)
    else if (pevt[i].daughter1() < pevt[i].daughter2() && pevt[i].daughter1() > 0){
      begin = pevt[i].daughter1();
      end = pevt[i].daughter2()+1;
    }
    // Two separately stored decay products (e.g. in backwards evolution of initial-state showers).
    else if (pevt[i].daughter1() > pevt[i].daughter2() && pevt[i].daughter2() > 0){
      THROW(fatal_error,"Two separetely stored decay products can not be handled at the moment.");
    }
    else {
      THROW(fatal_error,"Unknown configuration of daughgters.");
    }
    for (int d = begin; d < end; ++d) {
      auto find = m_processed.find(d);
      if (d > 0 && find == m_processed.end()){
        // Avoid handling this particle again. Would happen otherwise as in most cases each particle originiates from multiple initial partons at once.
        m_processed.insert(d);
        int Id = pevt[d].id();
        kf_code kfc = (kf_code) abs(IdToSherpa(Id));
        Flavour flav = Flavour(kfc, Id<0);
        const auto it = s_kftable.find(kfc);
        // Replace particles unknown to sherpa directly with their daughters.
        if (it == s_kftable.end()) {
          msg_Error() << "Sherpa does not know particle " << m_pythia.particleData.name(abs(Id)) << " with id " << abs(IdToSherpa(Id)) << " ! Replacing it with its decay products! If Sherpa knows this particle under a different ID it should be added to the conversion functions!" << std::endl;
          HandleDaughters(bloblist, decayblob, pevt, d);
        }
        // Replace partons with hadronization products.
        // Happens in hadronic decays or when initial partons are combined or experience recoil effect.
        else if (abs(pevt[d].status()) < 80) {
          msg_Tracking() << "Particle " <<  m_pythia.particleData.name(abs(Id)) << " with id " << abs(Id) << " should be a parton in preparation for hadronization (status in 70s) or from final state shower of partonic decay (status in 50s)!" << std::endl;
          msg_Tracking() << "Continuing with its daughters."  << std::endl;
          HandleDaughters(bloblist, decayblob, pevt, d);
        }
        // Check for particles that should be stable but are only intermediary
        else if (flav.IsStable() && pevt[d].status() < 0) {
          // Is fine for gluons or quarks as these are from partonic decays and thus hadronizing and not unexpectedly decaying.
          if (Id == 21 || Id < 9) {
            msg_Tracking() << "Particle " <<  m_pythia.particleData.name(abs(Id)) << " with id " << abs(Id) << " is quark or gluon from partonic decay. Continuing with daughters." << std::endl;
          }
          // This should not happen.(And has not in testing)
          else {
            msg_Error() << "Particle " <<  m_pythia.particleData.name(abs(Id)) << " with id " << abs(Id) << " was supposed to be stable but is only intermediary." << std::endl;
            THROW(fatal_error,"Particle is supposed to be stable.");
          }
          HandleDaughters(bloblist, decayblob, pevt, d);
        }
        // If none of previous exceptions occur the particle is added as outgoing to the decay or fragmentation blob.
        else {
          Vec4D momentum = Vec4D(pevt[d].e(),pevt[d].px(),pevt[d].py(),pevt[d].pz());
          Vec4D position = Vec4D(pevt[d].tProd(),pevt[d].xProd(),pevt[d].yProd(),pevt[d].zProd());
          Particle *daughter = new Particle(-1,flav,momentum);
          if (decayblob->Type()==btp::Fragmentation){
            daughter->SetInfo('P');
          }
          else {
            daughter->SetInfo('D');
          }
          daughter->SetNumber(0);
          daughter->SetFinalMass(pevt[d].mCalc());
          decayblob->SetPosition(position);
          decayblob->AddToOutParticles(daughter);
          // Initialize the particles decay blob if appropriate.
          if (m_pythiadecays && (pevt[d].status() < 0)){
            HandleDecays(bloblist, pevt, daughter, d);
          }
          else {
            daughter->SetStatus(part_status::active);
          }
        }
      }
    }
  }

  void HandleDecays(Blob_List * bloblist, Pythia8::Event& pevt, Particle* inpart, int i)
  /*
    Create decay blob and fill outgoing particles based on daughters in Pythia event.
  */
  {
    if(inpart->DecayBlob()){
      msg_Error() << inpart->DecayBlob() << std::endl;
      THROW(fatal_error,"Decay blob already exists.");
    }
    if(inpart->Flav().IsStable()){
      msg_Error() << inpart->Flav() << std::endl;
      pevt.list();
      THROW(fatal_error,"Particle is supposed to be stable.");
    }
    if(inpart->Time()==0.0) inpart->SetTime();
    inpart->SetStatus(part_status::decayed);
    Blob* blob = bloblist->AddBlob(btp::Hadron_Decay);
    blob->AddToInParticles(inpart);
    blob->SetTypeSpec("Pythia8");
    DEBUG_VAR(inpart->Momentum());
    HandleDaughters(bloblist, blob, pevt, i);
  }

  void AssignDecays() {
    /*
      Set variable for what should handle decays.
      If Sherpa does them they need to be turned off for Pythia and BreitWigner smearing needs to be turned off.
    */
    m_pythiadecays = m_settings["PYTHIA8"]["DECAYS"].SetDefault(true).Get<bool>();
    if  (!m_pythiadecays) {
      m_pythia.readString("HadronLevel:Decay = off");
      PRINT_INFO("Setting particles on-shell to allow sherpa decays.");
      m_pythia.readString("ParticleData:modeBreitWigner = 0");
    }
  }

  void ApplyPythiaSettings() {
    /*
      Apply all Pythia settings that are adjusted in the yaml via the readString mechanism.
     */
    PRINT_INFO("Applying Pythia8 settings");
    m_settings["PYTHIA8"]["PARAMETERS"].SetDefault("");
    for (auto& proc : m_settings["PYTHIA8"]["PARAMETERS"].GetItems()) {
      auto keys = proc.GetKeys();
      if (keys.size() != 1) {
    	if (!msg_LevelIsTracking()) msg_Info()<<"\n";
    	THROW(invalid_input, std::string{"Invalid Pythia8 setting.\n\n"});
      }
      auto pythiasetting = proc[keys[0]];
      std::string value = pythiasetting.SetDefault("").GetScalar<std::string>();
      std::string name = keys[0];
      m_pythia.readString(name+" = "+value);
    }
  }

  void HarmonizeMasses() {
    /*
      Harmonize particle settings between Sherpa and Pythia.
      Settings for one are adjusted to match those of the other depending on which option is chosen.
      Also possible to only adjust the settings for some particles.
     */
    bool SherpaValues;
    SherpaValues = m_settings["PYTHIA8"]["SHERPA_MASSES"].SetDefault(true).Get<bool>();
    bool MatchQuarks = m_settings["PYTHIA8"]["MATCH_QUARKS"].SetDefault(true).Get<bool>();
    bool MatchDiQuarks = m_settings["PYTHIA8"]["MATCH_DIQUARKS"].SetDefault(true).Get<bool>();
    bool MatchHadrons = m_settings["PYTHIA8"]["MATCH_HADRONS"].SetDefault(true).Get<bool>();
    bool MatchOther;     // Leptons + Bosons
    if (SherpaValues) {
      MatchOther = m_settings["PYTHIA8"]["MATCH_OTHER"].SetDefault(true).Get<bool>();
    }
    else {
      MatchOther = m_settings["PYTHIA8"]["MATCH_OTHER"].SetDefault(false).Get<bool>();
    }
    bool MatchOnlyUnstable = m_settings["PYTHIA8"]["MATCH_ONLY_UNSTABLE"].SetDefault(false).Get<bool>();

    PRINT_INFO("Harmonizing particle masses and widths!");
    if (SherpaValues){
      ReadMassParameters();
      ModifyPythiaValues(MatchQuarks,MatchDiQuarks,MatchHadrons,MatchOther,MatchOnlyUnstable);
    }
    else {
      ModifySherpaValues(MatchQuarks,MatchDiQuarks,MatchHadrons,MatchOther,MatchOnlyUnstable);
    }
  }

  void ModifyPythiaValues(bool MatchQuarks,bool MatchDiQuarks,bool MatchHadrons,bool MatchOther,bool MatchOnlyUnstable) {
    PRINT_INFO("Changing Pythia Values");
    for(KFCode_ParticleInfo_Map::const_iterator kfit(s_kftable.begin());
        kfit!=s_kftable.end();++kfit) {
      Flavour flav(kfit->first);
      int PythiaID = IdToPythia(flav.Kfcode());
      bool MatchParticleType = (((flav.IsHadron() && MatchHadrons) || (flav.IsQuark() && MatchQuarks) ||
                                 ((flav.IsLepton() || flav.IsBoson()) && MatchOther) && flav.IsOn()) ||
                                (flav.IsDiQuark() && MatchDiQuarks));
      bool MatchParticleConditions = (!flav.IsDummy() && flav.Size()==1 && flav.Kfcode()!=0);
      bool MatchParticleStability = (!MatchOnlyUnstable || (MatchOnlyUnstable && !flav.IsStable()));
      if (MatchParticleType && MatchParticleConditions && MatchParticleStability) {
        if (m_pythia.particleData.m0(PythiaID) && !(abs(flav.HadMass()-m_pythia.particleData.m0(PythiaID))/m_pythia.particleData.m0(PythiaID) < 1.e-2) ) {
          msg_Tracking()<<METHOD<<" Adjusted mass of "<<flav<<" ("<<flav.Kfcode()
                        <<") from "<<m_pythia.particleData.m0(PythiaID)<<" to "<<flav.HadMass()<<"."<<endl;
        }
        m_pythia.particleData.m0(PythiaID, flav.HadMass());
        if (m_pythia.particleData.mWidth(PythiaID)){
          m_pythia.particleData.mWidth(PythiaID, flav.Width());
        }
        m_pythia.particleData.mayDecay(PythiaID, !(flav.IsStable()));
      }
    }
  }

  void ModifySherpaValues(bool MatchQuarks,bool MatchDiQuarks,bool MatchHadrons,bool MatchOther,bool MatchOnlyUnstable) {
    PRINT_INFO("Changing Sherpa Values");
    for(KFCode_ParticleInfo_Map::const_iterator kfit(s_kftable.begin());
        kfit!=s_kftable.end();++kfit) {
      Flavour flav(kfit->first);
      int PythiaID = IdToPythia(flav.Kfcode());
      bool MatchParticleType = (((flav.IsHadron() && MatchHadrons) || (flav.IsQuark() && MatchQuarks) ||
                                 ((flav.IsLepton() || flav.IsBoson()) && MatchOther) && flav.IsOn()) ||
                                (flav.IsDiQuark() && MatchDiQuarks));
      bool MatchParticleConditions = (!flav.IsDummy() && flav.Size()==1 && flav.Kfcode()!=0);
      bool MatchParticleStability = (!MatchOnlyUnstable || (MatchOnlyUnstable && !flav.IsStable()));
      if (MatchParticleType && MatchParticleConditions && MatchParticleStability && m_pythia.particleData.isParticle(flav.Kfcode())) {
        if (flav.HadMass() && !(abs(flav.HadMass()-m_pythia.particleData.m0(PythiaID))/flav.HadMass() < 1.e-2) ) {
          msg_Tracking()<<METHOD<<" Adjusted mass of "<<flav<<" ("<<flav.Kfcode()
                        <<") from "<<flav.HadMass()<<" to "<<m_pythia.particleData.m0(PythiaID)<<"."<<endl;
        }
        flav.SetHadMass(m_pythia.particleData.m0(PythiaID));
        flav.SetMass(m_pythia.particleData.m0(PythiaID));
        if (m_pythia.particleData.mWidth(PythiaID)){
          flav.SetWidth(m_pythia.particleData.mWidth(PythiaID));
        }
        if (m_pythia.settings.flag("ParticleDecays:limitTau0") && m_pythia.particleData.tau0(PythiaID)>m_pythia.settings.parm("ParticleDecays:tau0Max")) {
          m_pythia.particleData.mayDecay(PythiaID, false);
        }
        flav.SetStable(!(m_pythia.particleData.mayDecay(PythiaID) && m_pythia.particleData.canDecay(PythiaID)));
      }
    }
  }

  int IdToPythia(const int& ID) {
    switch (ID) {
    case 102142: return 14122;
    case -102142: return -14122;
    case 102144: return 4124;
    case -102144: return -4124;
    default: return ID;
    }
  }

  int IdToSherpa(const int& ID) {
    switch (ID) {
    case 14122: return 102142;
    case -14122: return -102142;
    case 4124: return 102144;
    case -4124: return -102144;
    default: return ID;
    }
  }

  Pythia8::Pythia m_pythia;
  bool m_pythiadecays;
  Settings& m_settings = Settings::GetMainSettings();
  std::set<int> m_processed;
};

}

DEFINE_FRAGMENTATION_GETTER(SHERPA::Pythia8_Hadronisation, "Pythia8");
