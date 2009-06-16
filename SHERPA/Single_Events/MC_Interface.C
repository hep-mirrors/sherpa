#include "SHERPA/Single_Events/MC_Interface.H"
#include "ATOOLS/Org/Message.H"

#ifdef PROFILE__Jet_Evolution
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif


using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

MC_Interface::MC_Interface(Lund_Interface * lund) :
  p_lund(lund), m_mode(1)
{
  m_name      = std::string("MC_Interface:Pythia");
  m_type      = eph::External_MC;
}

    
MC_Interface::~MC_Interface() { }


Return_Value::code MC_Interface::Treat(Blob_List * blobs, double & weight)
{
  PROFILE_LOCAL("MC_Interface::Treat");
  if (blobs->size()>0) return Return_Value::Nothing;
  switch (m_mode) {
    case 1: return p_lund->OneEvent(blobs,weight);
  }
  return Return_Value::Nothing;
}


void MC_Interface::CleanUp() { }

void MC_Interface::Finish(const std::string &) {}
