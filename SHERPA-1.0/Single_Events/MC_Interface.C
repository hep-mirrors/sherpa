#include "MC_Interface.H"

#ifdef PROFILE__Jet_Evolution
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif


using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

MC_Interface::MC_Interface(Lund_Interface * _pythia) :
  p_pythia(_pythia)
{
  m_name      = std::string("MC_Interface:Pythia");
  m_type      = eph::External_MC;
}

MC_Interface::~MC_Interface() { }


bool MC_Interface::Treat(Blob_List * blobs, double & weight)
{
  PROFILE_LOCAL("MC_Interface::Treat");
  if (blobs->size()>0) return false;
  return p_pythia->OneEvent(blobs,weight);
}


void MC_Interface::CleanUp() { }

void MC_Interface::Finish(const std::string &) {}
