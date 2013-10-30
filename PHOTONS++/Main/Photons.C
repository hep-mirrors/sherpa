#include "PHOTONS++/Main/Photons.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Phys/Blob.H"
#include "PHOTONS++/Main/Define_Dipole.H"


using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

// define statics
int    PHOTONS::Photons::s_mode          = 2;
bool   PHOTONS::Photons::s_useme         = true;
double PHOTONS::Photons::s_ircutoff      = 1E-3;
int    PHOTONS::Photons::s_ircutoffframe = 0;
double PHOTONS::Photons::s_accu          = 1E-6;
int    PHOTONS::Photons::s_nmax          = std::numeric_limits<double>::max();
int    PHOTONS::Photons::s_nmin          = 0;
double PHOTONS::Photons::s_drcut         = 1000.;
bool   PHOTONS::Photons::s_strict        = false;
double PHOTONS::Photons::s_reducemax     = 1.;
bool   PHOTONS::Photons::s_checkfirst    = false;

double PHOTONS::Photons::s_alpha                = 0.;
bool   PHOTONS::Photons::s_userunningparameters = false;

// member functions of class Photons

Photons::Photons(Data_Reader* reader) :
  m_name("Photons")
{
  rpa->gen.AddCitation
    (1,"Photons is published under \\cite{Schonherr:2008av}.");
  s_mode          = reader->GetValue<int>("YFS_MODE",2);
  if (s_mode>2) s_mode=2;
  s_useme         = (bool)reader->GetValue<int>("YFS_USE_ME",1);
  s_ircutoff      = reader->GetValue<double>("YFS_IR_CUTOFF",1E-3);
  s_userunningparameters = (bool)reader->GetValue<int>("YFS_USE_RUNNING_PARAMETERS",0);
  std::string irframe
            = reader->GetValue<std::string>("YFS_IR_CUTOFF_FRAME","Multipole_CMS");
  if      (irframe == "Multipole_CMS")      s_ircutoffframe = 0;
  else if (irframe == "Lab")                s_ircutoffframe = 1;
  else if (irframe == "Decayer_Rest_Frame") s_ircutoffframe = 2;
  else {
    s_ircutoffframe = 0;
    msg_Info()<<"value '"<<irframe<<"' for the frame for applying the\n"
              <<"IR cut-off for soft photon radiation unkown ...\n"
              <<"setting it to 'Multipole_CMS' ...\n";
  }
  s_nmax          = reader->GetValue<int>("YFS_MAXEM",-1);
  s_nmin          = reader->GetValue<int>("YFS_MINEM",0);
  if (s_nmax<0) s_nmax = std::numeric_limits<double>::max();
  s_nmin          = reader->GetValue<int>("YFS_MINEM",0);
  s_drcut         = reader->GetValue<double>("YFS_DRCUT",1000.);
  s_strict        = reader->GetValue<int>("YFS_STRICTNESS",0);
  s_reducemax     = reader->GetValue<double>("YFS_REDUCE_MAXIMUM",1.);
  s_checkfirst    = (bool)reader->GetValue<double>("YFS_CHECK_FIRST",0);
  s_accu          = sqrt(rpa->gen.Accu());
  m_success       = true;
  m_photonsadded  = false;
  msg_Debugging()<<METHOD<<"(){\n"
		 <<"  Mode: "<<s_mode
		 <<" ,  MEs: "<<(s_mode>1?s_useme:0)
		 <<" ,  nmax: "<<s_nmax
		 <<" ,  nmin: "<<s_nmin
		 <<" ,  strict: "<<s_strict
		 <<" ,  dRcut: "<<s_drcut
		 <<" ,  reducemax: "<<s_reducemax
		 <<" ,  IR cut-off: "<<(s_mode>0?s_ircutoff:0)
		 <<" in frame "<<irframe<<" ("<<s_ircutoffframe<<")"
		 <<" ,  use running parameters "<<s_userunningparameters
		 <<"\n}"<<std::endl;
}

Photons::Photons() :
  m_name("Photons")
{
  PRINT_INFO("TODO: check whether running of alphaQED is MSbar");
  PRINT_INFO("TODO: evolve all particle masses in MSbar");
  s_mode          = 2;
  s_useme         = true;
  s_ircutoff      = 1E-3;
  s_ircutoffframe = 0;
  s_nmax          = std::numeric_limits<double>::max();
  s_nmin          = 0;
  s_drcut         = 1000.;
  s_strict        = false;
  s_reducemax     = 1.;
  s_checkfirst    = false;
  s_accu          = sqrt(rpa->gen.Accu());
  s_userunningparameters = false;
  m_success       = true;
  m_photonsadded  = false;
  msg_Debugging()<<METHOD<<"(){\n"
		 <<"  Mode: "<<s_mode
		 <<" ,  MEs: "<<(s_mode>1?s_useme:0)
		 <<" ,  nmax: "<<s_nmax
		 <<" ,  nmin: "<<s_nmin
		 <<" ,  strict: "<<s_strict
		 <<" ,  dRcut: "<<s_drcut
		 <<" ,  reducemax: "<<s_reducemax
		 <<" ,  IR cut-off: "<<(s_mode>0?s_ircutoff:0)
		 <<" in frame "<<s_ircutoffframe
		 <<" ,  use running parameters "<<s_userunningparameters
		 <<"\n}"<<std::endl;
}

bool Photons::AddRadiation(Blob * blob)
{
  if (!CheckStateBeforeTreatment(blob)) {
    m_photonsadded=false;
    return m_success=false;
  }
  ResetAlphaQED();
  Define_Dipole dress(blob);
  if (!dress.CheckMasses()) {
    msg_Error()<<METHOD<<"(): Found massless charged particles. Cannot cope."
               <<std::endl;
    m_photonsadded=false;
    return m_success=false;
  }
  dress.AddRadiation();
  m_photonsadded = dress.AddedAnything();
  m_success = dress.DoneSuccessfully();
  if (!blob->MomentumConserved()) {
    msg_Error()<<METHOD<<"(): Momentum not conserved after treatment: "
               <<blob->CheckMomentumConservation()<<std::endl
               <<*blob<<std::endl;
    return m_success=false;
  }
  return m_success;
}

bool Photons::CheckStateBeforeTreatment(Blob * blob)
{
  if (!s_checkfirst) return true;
  DEBUG_FUNC(blob->ShortProcessName());
  if (!blob->MomentumConserved()) {
    msg_Error()<<METHOD<<"(): Momentum not conserved before treatment: "
               <<blob->CheckMomentumConservation()<<std::endl
               <<*blob<<std::endl;
    return false;
  }
  bool rightmasses(true);
  for (size_t i(0);i<blob->NOutP();++i) {
    if (blob->OutParticle(i)->FinalMass()==0. &&
        blob->OutParticle(i)->Momentum().Mass()>1e-3) {
      msg_Error()<<METHOD<<"(): "<<blob->OutParticle(i)->Flav().IDName()
                         <<" not onshell: "
                         <<blob->OutParticle(i)->Momentum().Mass()
                         <<" vs "<<blob->OutParticle(i)->FinalMass()<<std::endl;
      rightmasses=false;
    }
    else if (blob->OutParticle(i)->FinalMass()!=0. &&
             !IsEqual(blob->OutParticle(i)->Momentum().Mass(),
                      blob->OutParticle(i)->FinalMass(),1e-3)) {
      msg_Error()<<METHOD<<"(): "<<blob->OutParticle(i)->Flav().IDName()
                         <<" not onshell: "
                         <<blob->OutParticle(i)->Momentum().Mass()
                         <<" vs "<<blob->OutParticle(i)->FinalMass()<<std::endl;
      rightmasses=false;
    }
  }
  if (!rightmasses) return false;
  return true;
}

