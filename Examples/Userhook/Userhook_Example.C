#include "SHERPA/Tools/Userhook_Base.H"
#include "ATOOLS/Org/Message.H"
#include "SHERPA/Main/Sherpa.H"
#include "SHERPA/Tools/HepMC3_Interface.H"

#ifdef USING__HEPMC3
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"
#endif

using namespace ATOOLS;
using namespace SHERPA;
using namespace std;

class Userhook_Example : public Userhook_Base {

  Sherpa* p_sherpa;
  size_t  m_nevents, m_nvertices, m_nparticles;

public:

  Userhook_Example(const Userhook_Arguments args) :
    Userhook_Base("Example"), p_sherpa(args.p_sherpa),
    m_nevents(0), m_nvertices(0), m_nparticles(0)
  {
    PRINT_INFO("We are using a user hook within Sherpa and are using PDF "<<p_sherpa->PDFInfo());
  }

  ~Userhook_Example() {}

  ATOOLS::Return_Value::code Run(ATOOLS::Blob_List* blobs) {
    DEBUG_INFO("Let's do something with the bloblist for each event:");

    ++m_nevents;
    m_nvertices += blobs->size();
    for (auto blob : *blobs) {
      m_nparticles += blob->OutParticles()->size();
    }

#ifdef USING__HEPMC3
    // and now convert into a HepMC event
    SHERPA::HepMC3_Interface hepmc_converter;
    HepMC3::GenEvent hepmc_event;
    hepmc_converter.Sherpa2HepMC(blobs, hepmc_event);
    // ... do something with HepMC event here ...
    DEBUG_VAR(hepmc_event.particles().size());
#endif

    if(blobs->FourMomentumConservation()) {
      return Return_Value::Nothing;
    }
    else {
      return Return_Value::Error;
    }
  }

  void Finish() {
    PRINT_INFO("End of the run... "  << endl <<
               "  Number of events:  " << m_nevents << endl <<
               "  Average number of vertices per event: " << double(m_nvertices)/double(m_nevents) << endl <<
               "  Average number of particles per event: " << double(m_nparticles)/double(m_nevents) << endl);
  }

};

DECLARE_USERHOOK_GETTER(Userhook_Example, "Example")
