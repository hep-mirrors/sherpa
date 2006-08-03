#include "MC_Interface.H"
#include "Message.H"

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

#ifdef USING__MCatNLO
MC_Interface::MC_Interface(Herwig_Interface * _herwig) :
  p_pythia(NULL), p_herwig(_herwig), p_mcatnlo(NULL), m_mode(2)
{
  m_name      = std::string("MC_Interface:Herwig");
  m_type      = eph::External_MC;
}

MC_Interface::MC_Interface(MCatNLO_Interface * _mcatnlo) :
  p_pythia(NULL), p_herwig(NULL), p_mcatnlo(_mcatnlo), m_mode(3)
{
  m_name      = std::string("MC_Interface:MC@NLO");
  m_type      = eph::External_MC;
}
#endif
    
MC_Interface::~MC_Interface() { }


Return_Value::code MC_Interface::Treat(Blob_List * blobs, double & weight)
{
  PROFILE_LOCAL("MC_Interface::Treat");
  if (blobs->size()>0) return Return_Value::Nothing;
  switch (m_mode) {
    case 1: return p_lund->OneEvent(blobs,weight);
#ifdef USING__MCatNLO
    case 2: return p_herwig->OneEvent(blobs,weight);
    case 3: return p_mcatnlo->OneEvent(blobs,weight);
#endif
  }
  return Return_Value::Nothing;
}


void MC_Interface::CleanUp() { }

void MC_Interface::Finish(const std::string &) 
{
  switch (m_mode) {
#ifdef USING__MCatNLO
  case 3: return p_mcatnlo->Terminate();
#endif
  }
}
