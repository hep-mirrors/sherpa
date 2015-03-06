#include "PHOTONS++/Main/Photons.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Phys/Blob.H"
#include "PHOTONS++/Main/Define_Dipole.H"

#ifdef PHOTONS_DEBUG
#include "ATOOLS/Math/Histogram_2D.H"
#include "ATOOLS/Org/Shell_Tools.H"
#endif

using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

// define statics
int    PHOTONS::Photons::s_mode          = 2;
bool   PHOTONS::Photons::s_useme         = true;
double PHOTONS::Photons::s_ircutoff      = 1E-3;
double PHOTONS::Photons::s_uvcutoff      = std::numeric_limits<double>::max();
int    PHOTONS::Photons::s_ircutoffframe = 0;
double PHOTONS::Photons::s_accu          = 1E-6;
int    PHOTONS::Photons::s_nmax          = std::numeric_limits<int>::max();
int    PHOTONS::Photons::s_nmin          = 0;
double PHOTONS::Photons::s_drcut         = 1000.;
bool   PHOTONS::Photons::s_strict        = false;
double PHOTONS::Photons::s_reducemaxenergy   = 1.;
double PHOTONS::Photons::s_increasemaxweight = 1.;
bool   PHOTONS::Photons::s_checkfirst    = false;
int    PHOTONS::Photons::s_ffrecscheme   = 0;
int    PHOTONS::Photons::s_firecscheme   = 0;

double PHOTONS::Photons::s_alpha                = 0.;
bool   PHOTONS::Photons::s_userunningparameters = false;

#ifdef PHOTONS_DEBUG
std::string  PHOTONS::Photons::s_histo_base_name("weights");
Histogram_2D PHOTONS::Photons::s_histo_dipole
                = Histogram_2D(101,1e-6,1e4,100,-0.5,10.5,11);
Histogram_2D PHOTONS::Photons::s_histo_jacobianM
                = Histogram_2D(101,1e-6,1e4,100,-0.5,10.5,11);
Histogram_2D PHOTONS::Photons::s_histo_jacobianL
                = Histogram_2D(101,1e-6,1e4,100,-0.5,10.5,11);
Histogram_2D PHOTONS::Photons::s_histo_higher
                = Histogram_2D(101,1e-6,1e4,100,-0.5,10.5,11);
Histogram_2D PHOTONS::Photons::s_histo_yfs
                = Histogram_2D(101,1e-6,1e4,100,-0.5,10.5,11);
Histogram_2D PHOTONS::Photons::s_histo_total
                = Histogram_2D(101,1e-6,1e4,100,-0.5,10.5,11);
Histogram_2D PHOTONS::Photons::s_histo_t_dipole
                = Histogram_2D(111,1e-6,1e4,100,1e-6,1e2,20);
Histogram_2D PHOTONS::Photons::s_histo_t_jacobianM
                = Histogram_2D(111,1e-6,1e4,100,1e-6,1e2,20);
Histogram_2D PHOTONS::Photons::s_histo_t_jacobianL
                = Histogram_2D(111,1e-6,1e4,100,1e-6,1e2,20);
Histogram_2D PHOTONS::Photons::s_histo_t_higher
                = Histogram_2D(111,1e-6,1e4,100,1e-6,1e2,20);
Histogram_2D PHOTONS::Photons::s_histo_t_yfs
                = Histogram_2D(111,1e-6,1e4,100,1e-6,1e2,20);
Histogram_2D PHOTONS::Photons::s_histo_t_total
                = Histogram_2D(111,1e-6,1e4,100,1e-6,1e2,20);
#endif


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
  s_uvcutoff      = reader->GetValue<double>("YFS_UV_CUTOFF",-1.);
  if (s_uvcutoff<0.) s_uvcutoff = std::numeric_limits<double>::max();
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
  if (s_nmax<0) s_nmax = std::numeric_limits<int>::max();
  s_nmin          = reader->GetValue<int>("YFS_MINEM",0);
  s_drcut         = reader->GetValue<double>("YFS_DRCUT",1000.);
  s_strict        = reader->GetValue<int>("YFS_STRICTNESS",0);
  s_reducemaxenergy = reader->GetValue<double>("YFS_REDUCE_MAXIMUM_ENERGY",1.);
  s_increasemaxweight = reader->GetValue<double>("YFS_INCREASE_MAXIMUM_WEIGHT",1.);
  s_checkfirst    = (bool)reader->GetValue<double>("YFS_CHECK_FIRST",0);
  s_ffrecscheme   = reader->GetValue<int>("YFS_FF_RECOIL_SCHEME",2);
  s_firecscheme   = reader->GetValue<int>("YFS_FI_RECOIL_SCHEME",2);
  s_accu          = sqrt(rpa->gen.Accu());
#ifdef PHOTONS_DEBUG
  s_histo_base_name = reader->GetValue<std::string>("YFS_HISTO_BASE_NAME",
                                                    "weights");
#endif
  m_success       = true;
  m_photonsadded  = false;
  msg_Debugging()<<METHOD<<"(){\n"
		 <<"  Mode: "<<s_mode
		 <<" ,  MEs: "<<(s_mode>1?s_useme:0)
		 <<" ,  nmax: "<<s_nmax
		 <<" ,  nmin: "<<s_nmin
		 <<" ,  strict: "<<s_strict
		 <<" ,  dRcut: "<<s_drcut
		 <<" ,  reducemaxenergy: "<<s_reducemaxenergy
		 <<" ,  increasemaxweight: "<<s_increasemaxweight
		 <<" ,  IR cut-off: "<<(s_mode>0?s_ircutoff:0)
		 <<" in frame "<<irframe<<" ("<<s_ircutoffframe<<")"
		 <<" ,  UV cut-off: "<<s_uvcutoff
		 <<" ,  use running parameters "<<s_userunningparameters
		 <<" ,  FF recoil scheme: "<<s_ffrecscheme
		 <<" ,  FI recoil scheme: "<<s_firecscheme
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
  s_uvcutoff      = std::numeric_limits<double>::max();
  s_ircutoffframe = 0;
  s_nmax          = std::numeric_limits<int>::max();
  s_nmin          = 0;
  s_drcut         = 1000.;
  s_strict        = false;
  s_reducemaxenergy = 1.;
  s_increasemaxweight = 1.;
  s_checkfirst    = false;
  s_ffrecscheme   = 0;
  s_firecscheme   = 0;
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
		 <<" ,  reducemaxenergy: "<<s_reducemaxenergy
		 <<" ,  increasemaxweight: "<<s_increasemaxweight
		 <<" ,  IR cut-off: "<<(s_mode>0?s_ircutoff:0)
		 <<" in frame "<<s_ircutoffframe
		 <<" ,  UV cut-off: "<<s_uvcutoff
		 <<" ,  use running parameters "<<s_userunningparameters
		 <<" ,  FF recoil scheme: "<<s_ffrecscheme
		 <<" ,  FI recoil scheme: "<<s_firecscheme
		 <<"\n}"<<std::endl;
}

Photons::~Photons()
{
#ifdef PHOTONS_DEBUG
  size_t pos(s_histo_base_name.find_last_of("/"));
  if (pos!=std::string::npos) MakeDir(s_histo_base_name.substr(0,pos));
  s_histo_dipole.Finalize();
  s_histo_dipole.Output(s_histo_base_name+"-dipole-n.dat");
  s_histo_jacobianM.Finalize();
  s_histo_jacobianM.Output(s_histo_base_name+"-jacobianM-n.dat");
  s_histo_jacobianL.Finalize();
  s_histo_jacobianL.Output(s_histo_base_name+"-jacobianL-n.dat");
  s_histo_higher.Finalize();
  s_histo_higher.Output(s_histo_base_name+"-higher-n.dat");
  s_histo_yfs.Finalize();
  s_histo_yfs.Output(s_histo_base_name+"-yfs-n.dat");
  s_histo_total.Finalize();
  s_histo_total.Output(s_histo_base_name+"-total-n.dat");
  s_histo_t_dipole.Finalize();
  s_histo_t_dipole.Output(s_histo_base_name+"-dipole-t.dat");
  s_histo_t_jacobianM.Finalize();
  s_histo_t_jacobianM.Output(s_histo_base_name+"-jacobianM-t.dat");
  s_histo_t_jacobianL.Finalize();
  s_histo_t_jacobianL.Output(s_histo_base_name+"-jacobianL-t.dat");
  s_histo_t_higher.Finalize();
  s_histo_t_higher.Output(s_histo_base_name+"-higher-t.dat");
  s_histo_t_yfs.Finalize();
  s_histo_t_yfs.Output(s_histo_base_name+"-yfs-t.dat");
  s_histo_t_total.Finalize();
  s_histo_t_total.Output(s_histo_base_name+"-total-t.dat");
#endif
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
               <<blob->CheckMomentumConservation()<<std::endl;
    msg_Debugging()<<*blob<<std::endl;
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
               <<blob->CheckMomentumConservation()<<std::endl;
    msg_Debugging()<<*blob<<std::endl;
    return false;
  }
  bool rightmasses(true);
  for (size_t i(0);i<blob->NOutP();++i) {
    if (blob->OutParticle(i)->FinalMass()==0. &&
        blob->OutParticle(i)->Momentum().Mass()>1e-3) {
      msg_Debugging()<<METHOD<<"(): "<<blob->OutParticle(i)->Flav().IDName()
                             <<" not onshell: "
                             <<blob->OutParticle(i)->Momentum().Mass()<<" vs "
                             <<blob->OutParticle(i)->FinalMass()<<std::endl;
      rightmasses=false;
    }
    else if (blob->OutParticle(i)->FinalMass()!=0. &&
             !IsEqual(blob->OutParticle(i)->Momentum().Mass(),
                      blob->OutParticle(i)->FinalMass(),1e-3)) {
      msg_Debugging()<<METHOD<<"(): "<<blob->OutParticle(i)->Flav().IDName()
                             <<" not onshell: "
                             <<blob->OutParticle(i)->Momentum().Mass()<<" vs "
                             <<blob->OutParticle(i)->FinalMass()<<std::endl;
      rightmasses=false;
    }
  }
  if (!rightmasses) {
    msg_Error()<<METHOD<<"(): Particle(s) not on their mass shell. Cannot cope."
               <<std::endl;
    return false;
  }
  return true;
}

