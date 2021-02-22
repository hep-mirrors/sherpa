#include "ATOOLS/Phys/Fragmentation_Base.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Phys/KF_Table.H"

#include "Pythia8/Pythia.h"


using namespace ATOOLS;
using namespace std;

namespace SHERPA {

class Pythia8_Hadronisation : public Fragmentation_Base {

public:
  Pythia8_Hadronisation(const string& shower)
  {
    PRINT_INFO("Initialising Pythia8 hadronisation interface");

    // Initialise Pythia object
    m_pythia.readString("ProcessLevel:all = off");

    // Optionally switch off resonance decays, or only showers in them.
    m_pythia.readString("ProcessLevel:resonanceDecays = off");
    m_pythia.readString("PartonLevel:FSRinResonances = off");

    // Switch off automatic event listing in favour of manual.
    m_pythia.readString("Next:numberShowLHA = 0");
    m_pythia.readString("Next:numberShowInfo = 0");
    m_pythia.readString("Next:numberShowProcess = 0");
    m_pythia.readString("Next:numberShowEvent = 0");
    // m_pythia.readString("Next:numberCount = 1");

    // Optionally switch off ordinary decays.
    AssignDecays();

    ApplyPythiaSettings();

    // Settings for compressing and flagging partonic decays
    m_shrink = m_settings["COMPRESS_PARTONIC_DECAYS"].SetDefault(true).Get<bool>();
    m_flagpartonics = m_settings["FLAG_PARTONIC_DECAYS"].SetDefault(true).Get<bool>();

    HarmonizeMasses();

    // m_pythia.particleData.listAll();

    // ATOOLS::OutputParticles(msg->Info());
    // ATOOLS::OutputHadrons(msg->Info());

    m_pythia.init();

  }

  ~Pythia8_Hadronisation()
  {
    m_pythia.stat();
  }

  Return_Value::code Hadronize(Blob_List * blobs)
  {
    // PRINT_INFO("HADRONIZE WAS CALLED");
    Pythia8::Event& event      = m_pythia.event;
    for (Blob_List::iterator blit=blobs->begin();blit!=blobs->end();++blit) {
      // std::cout << "Blob status: " << (*blit)->Status() << std::endl;
      if ((*blit)->Has(blob_status::needs_hadronization) &&
	  (*blit)->Type()==btp::Fragmentation) {
	Blob * blob = (*blit);
	blob->SetTypeSpec("Pythia8");
        // std::cout << "Blob: " << std::endl;
        // std::cout << *(*blit) << std::endl;
	Sherpa2Pythia(blob, event);
	///Hadronization step
	if (!m_pythia.next()) {
          Blob * showerblob(blob->InParticle(0)->ProductionBlob());
          Blob * decblob(showerblob->InParticle(0)->ProductionBlob());
          if (decblob->Type()!=btp::Hadron_Decay) {
            // std::cout << "Bloblist:" << std::endl;
            // std::cout << *blobs << std::endl;
            // event.list();
            msg_Error()<<"Pythia8 hadronisation failed.\n"<<endl;
            return Return_Value::Error;
          }
          else {
            // std::cout << "Bloblist:" << std::endl;
            // std::cout << *blobs << std::endl;
            // event.list();
            msg_Tracking()<<"Error in "<<METHOD<<"."<<endl
                          <<"   Hadronization of partonic decay failed. Retry the event."<<endl;
            return Return_Value::Retry_Event;
          }
	}
	// event.list();
	FillFragmentationBlob(blob, event);
      }
    }
    // m_pythia.stat();
    // std::cout << "Bloblist:" << std::endl;
    // std::cout << *blobs << std::endl;
    if (m_shrink) Shrink(blobs);
    return Return_Value::Success;
  }

private:
  void Sherpa2Pythia(Blob * blob, Pythia8::Event& pevt)
  {
    // Reset event record to allow for new event.
    pevt.reset();
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
    The exception is when you want to include junctions in your string topology, i.e. have three string pieces meet.
    Then you must insert in your event record the (decayed) particle that is the reason for the presence of a junction,
    e.g. a baryon beam remnant from which several valence quarks have been kicked out, or a neutralino that underwent a baryon-number-violating decay.
    This particle must have as daughters the three partons that together carry the baryon number.
    */
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
      m = part->Momentum().Mass();
      pevt.append( id, status, col, acol, px, py, pz, e, m);
	}
    // pevt.list();
  }

  void FillFragmentationBlob(Blob * blob, Pythia8::Event& pevt)
  {
    Particle *particle;
    Flavour flav;
    Vec4D momentum, position;
    for (int i = 1; i < pevt.size(); ++i) {
      if (pevt[i].status() > 0) {
	kf_code kfc = (kf_code) abs(pevt[i].id());
	flav = Flavour(kfc, pevt[i].id()<0);
	momentum = Vec4D(pevt[i].e(),pevt[i].px(),pevt[i].py(),pevt[i].pz());
	position = Vec4D(pevt[i].tProd(),pevt[i].xProd(),pevt[i].yProd(),pevt[i].zProd());
	particle = new Particle(-1,flav,momentum);
	particle->SetNumber(0);
	particle->SetStatus(part_status::active);
	particle->SetInfo('P');
	particle->SetFinalMass(pevt[i].mCalc());
	blob->SetPosition(position);
	blob->AddToOutParticles(particle);
      }
    }
    blob->SetStatus(blob_status::needs_hadrondecays);
  }
  void AssignDecays() {
    bool pythiadecays = m_settings["PYTHIA8_DECAYS"].SetDefault(false).Get<bool>();
    if  (!pythiadecays) {
    m_pythia.readString("HadronLevel:Decay = off");
    PRINT_INFO("Setting particles on-shell to allow sherpa decays.");
    m_pythia.readString("ParticleData:modeBreitWigner = 0");
      }
  }

  void ApplyPythiaSettings() {
    PRINT_INFO("Applying Pythia8 settings");
    m_settings["PYTHIA8_PARAMETERS"].SetDefault("");
    for (auto& proc : m_settings["PYTHIA8_PARAMETERS"].GetItems()) {
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
    bool SherpaValues = m_settings["SHERPA_VALUES"].SetDefault(true).Get<bool>();
    bool MatchQuarks = m_settings["MATCH_QUARKS"].SetDefault(true).Get<bool>();
    bool MatchDiQuarks = m_settings["MATCH_DIQUARKS"].SetDefault(true).Get<bool>();
    bool MatchHadrons = m_settings["MATCH_HADRONS"].SetDefault(true).Get<bool>();
    bool MatchOther;     // Leptons + Bosons
    if (SherpaValues) {
     MatchOther = m_settings["MATCH_OTHER"].SetDefault(true).Get<bool>();
    }
    else {
     MatchOther = m_settings["MATCH_OTHER"].SetDefault(false).Get<bool>();
    }
    bool MatchOnlyUnstable = m_settings["MATCH_ONLY_UNSTALBE"].SetDefault(false).Get<bool>();

    // m_pythia.particleData.listAll();

    // ATOOLS::OutputParticles(msg->Info());
    // ATOOLS::OutputHadrons(msg->Info());

    PRINT_INFO("Harmonizing particle masses and widths!");
    if (SherpaValues){
      InitializeQuarkDiQuarkMasses();
      ModifyPythiaValues(MatchQuarks,MatchDiQuarks,MatchHadrons,MatchOther,MatchOnlyUnstable);
    }

    else {
      ModifySherpaValues(MatchQuarks,MatchDiQuarks,MatchHadrons,MatchOther,MatchOnlyUnstable);
    }
  }

  void InitializeQuarkDiQuarkMasses() {
      double mglue=  m_settings["HADRONIZATION_MASSES"]["M_GLUE"].SetDefault(0.00).Get<double>();
      double mud =   m_settings["HADRONIZATION_MASSES"]["M_UP_DOWN"].SetDefault(0.30).Get<double>();
      double ms =    m_settings["HADRONIZATION_MASSES"]["M_STRANGE"].SetDefault(0.40).Get<double>();
      double mc =    m_settings["HADRONIZATION_MASSES"]["M_CHARM"].SetDefault(1.80).Get<double>();
      double mb =    m_settings["HADRONIZATION_MASSES"]["M_BOTTOM"].SetDefault(5.10).Get<double>();
      double mdiq =  m_settings["HADRONIZATION_MASSES"]["M_DIQUARK_OFFSET"].SetDefault(0.30).Get<double>();
      double bind0 = m_settings["HADRONIZATION_MASSES"]["M_BIND_0"].SetDefault(0.12).Get<double>();
      double bind1 = m_settings["HADRONIZATION_MASSES"]["M_BIND_1"].SetDefault(0.50).Get<double>();
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


  void ModifyPythiaValues(bool MatchQuarks,bool MatchDiQuarks,bool MatchHadrons,bool MatchOther,bool MatchOnlyUnstable) {
    PRINT_INFO("Changing Pythia Values");
    for(KFCode_ParticleInfo_Map::const_iterator kfit(s_kftable.begin());
        kfit!=s_kftable.end();++kfit) {
      Flavour flav(kfit->first);
      bool MatchParticleType = (((flav.IsHadron() && MatchHadrons) || (flav.IsQuark() && MatchQuarks) ||
                                 ((flav.IsLepton() || flav.IsBoson()) && MatchOther) && flav.IsOn()) ||
                                (flav.IsDiQuark() && MatchDiQuarks));
      bool MatchParticleConditions = (!flav.IsDummy() && flav.Size()==1 && flav.Kfcode()!=0);
      bool MatchParticleStability = (!MatchOnlyUnstable || (MatchOnlyUnstable && !flav.IsStable()));
      if (MatchParticleType && MatchParticleConditions && MatchParticleStability) {
        if( !(abs(flav.HadMass()-m_pythia.particleData.m0(flav.Kfcode()))/m_pythia.particleData.m0(flav.Kfcode()) < 1.e-2) ) {
          msg_Tracking()<<METHOD<<" Adjusted mass of "<<flav<<" ("<<flav.Kfcode()
                        <<") from "<<m_pythia.particleData.m0(flav.Kfcode())<<" to "<<flav.HadMass()<<"."<<endl;
        }
        m_pythia.particleData.m0(flav.Kfcode(), flav.HadMass());
        m_pythia.particleData.mWidth(flav.Kfcode(), flav.Width());
      }
    }
  }

  void ModifySherpaValues(bool MatchQuarks,bool MatchDiQuarks,bool MatchHadrons,bool MatchOther,bool MatchOnlyUnstable) {
    PRINT_INFO("Changing Sherpa Values");
    for(KFCode_ParticleInfo_Map::const_iterator kfit(s_kftable.begin());
        kfit!=s_kftable.end();++kfit) {
      Flavour flav(kfit->first);
      bool MatchParticleType = (((flav.IsHadron() && MatchHadrons) || (flav.IsQuark() && MatchQuarks) ||
                                 ((flav.IsLepton() || flav.IsBoson()) && MatchOther) && flav.IsOn()) ||
                                (flav.IsDiQuark() && MatchDiQuarks));
      bool MatchParticleConditions = (!flav.IsDummy() && flav.Size()==1 && flav.Kfcode()!=0);
      bool MatchParticleStability = (!MatchOnlyUnstable || (MatchOnlyUnstable && !flav.IsStable()));
      if (MatchParticleType && MatchParticleConditions && MatchParticleStability && m_pythia.particleData.isParticle(flav.Kfcode())) {
        flav.SetHadMass(m_pythia.particleData.m0(flav.Kfcode()));
        flav.SetMass(m_pythia.particleData.m0(flav.Kfcode()));
        if( !(abs(flav.HadMass()-m_pythia.particleData.m0(flav.Kfcode()))/flav.HadMass() < 1.e-2) ) {
          msg_Tracking()<<METHOD<<" Adjusted mass of "<<flav<<" ("<<flav.Kfcode()
                        <<") from "<<flav.HadMass()<<" to "<<m_pythia.particleData.m0(flav.Kfcode())<<"."<<endl;
        }
        if (m_pythia.particleData.mWidth(flav.Kfcode())){
          flav.SetWidth(m_pythia.particleData.mWidth(flav.Kfcode()));
        }
      }
    }
  }

  // Shrink partonic blobs
  // Deletes the partonic decay products from the decay blob and replaces
  // with the results of their shower+fragmentation
  // Also removes the fragmentation and shower blob
  void Shrink(Blob_List * bloblist) {
    list<Blob *> deleteblobs;
    // Particle_Vector * parts;
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


  Pythia8::Pythia m_pythia;
  bool m_shrink;
  bool m_flagpartonics;
  Settings& m_settings = Settings::GetMainSettings();
};

}

DEFINE_FRAGMENTATION_GETTER(SHERPA::Pythia8_Hadronisation, "Pythia8");
