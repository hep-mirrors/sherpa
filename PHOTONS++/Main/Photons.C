#include "PHOTONS++/Main/Photons.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"


using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

// define statics
int    PHOTONS::Photons::s_mode     = 0;
bool   PHOTONS::Photons::s_useme    = 0;
double PHOTONS::Photons::s_ircutoff = 0;
double PHOTONS::Photons::s_accu     = 1E-6;

// member functions of class Photons

Photons::Photons(Data_Reader* reader, bool ana) :
  m_analyse(ana)
{
  rpa.gen.AddCitation
    (1,"Photons is published under \\cite{Schonherr:2008av}.");
  s_mode          = reader->GetValue<int>("YFS_MODE",2);
  s_useme         = (bool)reader->GetValue<int>("YFS_USE_ME",1);
  s_ircutoff      = reader->GetValue<double>("YFS_IR_CUTOFF",1E-3);
  s_accu          = sqrt(rpa.gen.Accu());
  m_success       = false;
  m_photonsadded  = false;
}

Photons::Photons(bool ana) :
  m_analyse(ana)
{
  rpa.gen.AddCitation
    (1,"Photons is published under \\cite{Schonherr:2008av}.");
  s_mode          = 2;
  s_useme         = true;
  s_ircutoff      = 1E-1;
  s_accu          = sqrt(rpa.gen.Accu());
  m_success       = false;
  m_photonsadded  = false;
}

Photons::~Photons() {
}

bool Photons::AddRadiation(Blob * blob) {

  if ((blob->Status() == blob_status::needs_extraQED) && (m_analyse == true)) {

    Define_Dipole dress(blob);
    dress.AddRadiation();
    m_photonsadded = dress.AddedAnything();
    m_success = dress.DoneSuccessfully();

#ifdef PHOTONS_DEBUG
    if (m_success == true)
      msg_Info()<<*blob<<endl;
    else
      msg_Info()<<"generation of YFS-radiation failed"<<endl;
#endif

  }
  blob->SetStatus(blob_status::needs_hadrondecays);
  return m_photonsadded;
}

