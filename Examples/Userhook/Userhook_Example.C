#include "SHERPA/Tools/Userhook_Base.H"
#include "ATOOLS/Org/Message.H"
#include "SHERPA/Main/Sherpa.H"

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

DECLARE_GETTER(Userhook_Example,"Example",
	       Userhook_Base,Userhook_Arguments);

Userhook_Base *ATOOLS::Getter<Userhook_Base,Userhook_Arguments,Userhook_Example>::
operator()(const Userhook_Arguments &args) const
{
  return new Userhook_Example(args);
}

void ATOOLS::Getter<Userhook_Base,Userhook_Arguments,Userhook_Example>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Example userhook";
}
