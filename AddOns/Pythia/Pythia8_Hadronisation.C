#include "ATOOLS/Phys/Fragmentation_Base.H"

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

    /// TODO

    m_pythia.init();
  }

  ~Pythia8_Hadronisation()
  {
    m_pythia.stat();
  }

  Return_Value::code Hadronize(Blob_List * blobs)
  {
    Pythia8::Event& event      = m_pythia.event;
    Sherpa2Pythia(blobs, event);
    /// TODO

    if (!m_pythia.next()) {
      msg_Error()<<"Pythia8 hadronisation failed."<<endl;
      return Return_Value::Error;
    }

    return Return_Value::Success;
  }

private:
  void Sherpa2Pythia(Blob_List* blobs, Pythia8::Event& pevt)
  {
    /// TODO
  }

  Pythia8::Pythia m_pythia;
};

}

DEFINE_FRAGMENTATION_GETTER(SHERPA::Pythia8_Hadronisation, "Pythia8");
