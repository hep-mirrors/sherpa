#include "ATOOLS/Phys/Fragmentation_Base.H"
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

    // Optionally switch off ordinary decays.
    m_pythia.readString("HadronLevel:Decay = off");

    // Switch off automatic event listing in favour of manual.
    m_pythia.readString("Next:numberShowInfo = 0");
    m_pythia.readString("Next:numberShowProcess = 0");
    m_pythia.readString("Next:numberShowEvent = 0");
    // m_pythia.readString("Next:numberCount = 1");

    Settings& s = Settings::GetMainSettings();

    PRINT_INFO("Applying Pythia8 settings");
    s["PYTHIA8_PARAMETERS"].SetDefault("");
    for (auto& proc : s["PYTHIA8_PARAMETERS"].GetItems()) {
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

    // m_pythia.particleData.listAll();

    // ATOOLS::OutputParticles(msg->Info());
    // ATOOLS::OutputHadrons(msg->Info());

    PRINT_INFO("Harmonizing particle masses and widths!");
    for(KFCode_ParticleInfo_Map::const_iterator kfit(s_kftable.begin());
	kfit!=s_kftable.end();++kfit) {
      Flavour flav(kfit->first);
      if ( (flav.IsHadron() || flav.IsDiQuark())  && (flav.Size()==1 && flav.Kfcode()!=0  && flav.IsOn() && !flav.Stable()) ) {
	// if( !(abs(flav.HadMass()-m_pythia.particleData.m0(flav.Kfcode()))/flav.HadMass() < 1.e-2) ) {
	//   msg_Info()<<METHOD<<" Adjusted mass of "<<flav<<" ("<<flav.Kfcode()
	// 	    <<") from "<<m_pythia.particleData.m0(flav.Kfcode())<<" to "<<flav.HadMass()<<"."<<endl;
	// }
	m_pythia.particleData.m0(flav.Kfcode(), flav.HadMass());
	m_pythia.particleData.mWidth(flav.Kfcode(), flav.Width());
	// m_pythia.particleData.mMin(flav.Kfcode(), flav.HadMass());
	// m_pythia.particleData.mMax(flav.Kfcode(), flav.HadMass());
      }
    }
    /// TODO
    PRINT_INFO("Setting particles on-shell to allow sherpa decays.");
    m_pythia.readString("ParticleData:modeBreitWigner = 0");

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
	  event.list();
	  msg_Error()<<"Pythia8 hadronisation failed.\n"<<endl;
	  return Return_Value::Error;
	}
	// event.list();
	FillFragmentationBlob(blob, event);
      }
    }
    // m_pythia.stat();

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

  Pythia8::Pythia m_pythia;
};

}

DEFINE_FRAGMENTATION_GETTER(SHERPA::Pythia8_Hadronisation, "Pythia8");
