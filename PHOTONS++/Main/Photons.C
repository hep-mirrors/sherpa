#include "PHOTONS++/Main/Photons.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/Blob.H"
#include "PHOTONS++/Main/Define_Dipole.H"

#ifdef PHOTONS_DEBUG
#include "ATOOLS/Math/Histogram_2D.H"
#include "ATOOLS/Org/Shell_Tools.H"
#endif

using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

std::ostream &PHOTONS::operator<<(std::ostream &str,const yfsmode::code &ym)
{
  if      (ym==yfsmode::off)  return str<<"Off";
  else if (ym==yfsmode::soft) return str<<"Soft";
  else if (ym==yfsmode::full) return str<<"Full";
  return str<<"unknown";
}

std::istream &PHOTONS::operator>>(std::istream &str,yfsmode::code &ym)
{
  std::string tag;
  str>>tag;
  ym=yfsmode::full;
  if      (tag.find("None")!=std::string::npos) ym=yfsmode::off;
  else if (tag.find("Soft")!=std::string::npos) ym=yfsmode::soft;
  else if (tag.find("1")!=std::string::npos)    ym=yfsmode::soft;
  else if (tag.find("Full")!=std::string::npos) ym=yfsmode::full;
  else if (tag.find("2")!=std::string::npos)    ym=yfsmode::full;
  return str;
}

// define statics
yfsmode::code PHOTONS::Photons::s_mode   = yfsmode::full;
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
double PHOTONS::Photons::s_alpha_input          = 0.;
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

Photons::Photons() :
  m_name("Photons")
{
  RegisterDefaults();

  Scoped_Settings s{ Settings::GetMainSettings()["YFS"] };
  rpa->gen.AddCitation(1,
                       "Photons is published under \\cite{Schonherr:2008av}.");
  s_mode     = s["MODE"].Get<yfsmode::code>();
  s_useme    = (bool)s["USE_ME"].Get<int>();
  s_ircutoff = s["IR_CUTOFF"].Get<double>();
  s_uvcutoff = s["UV_CUTOFF"].Get<double>();
  s_alpha_input   = s["1/ALPHAQED"].Get<double>();
  s_alpha_input = (s_alpha_input?1./s_alpha_input:MODEL::aqed->AqedThomson());
  s_userunningparameters = (bool)s["USE_RUNNING_PARAMETERS"].Get<int>();
  std::string irframe = s["IR_CUTOFF_FRAME"].Get<std::string>();
  if      (irframe == "Multipole_CMS")      s_ircutoffframe = 0;
  else if (irframe == "Lab")                s_ircutoffframe = 1;
  else if (irframe == "Decayer_Rest_Frame") s_ircutoffframe = 2;
  else {
    s_ircutoffframe = 0;
    msg_Info()<<"value '"<<irframe<<"' for the frame for applying the\n"
              <<"IR cut-off for soft photon radiation unkown ...\n"
              <<"setting it to 'Multipole_CMS' ...\n";
  }
  s_nmax          = s["MAXEM"].Get<int>();
  s_nmin          = s["MINEM"].Get<int>();
  s_drcut         = s["DRCUT"].Get<double>();
  s_strict        = s["STRICTNESS"].Get<int>();
  s_reducemaxenergy = s["REDUCE_MAXIMUM_ENERGY"].Get<double>();
  s_increasemaxweight = s["INCREASE_MAXIMUM_WEIGHT"].Get<double>();
  s_checkfirst    = (bool)s["CHECK_FIRST"].Get<double>();
  s_ffrecscheme   = s["FF_RECOIL_SCHEME"].Get<int>();
  s_firecscheme   = s["FI_RECOIL_SCHEME"].Get<int>();
  s_accu          = sqrt(rpa->gen.Accu());
  m_splitphotons  = s["PHOTON_SPLITTER_MODE"].Get<int>();

  m_photonsplitter = Photon_Splitter(m_splitphotons);

#ifdef PHOTONS_DEBUG
  s_histo_base_name = s["HISTO_BASE_NAME"].Get<std::string>();
#endif
  m_success       = true;
  m_photonsadded  = false;
  msg_Debugging()<<METHOD<<"(){\n"
                 <<"  Mode: "<<s_mode;
  if ((int)s_mode>0) {
    msg_Debugging()<<" ,  MEs: "<<((int)s_mode>1?s_useme:0)
                   <<" ,  nmax: "<<s_nmax
                   <<" ,  nmin: "<<s_nmin
                   <<" ,  strict: "<<s_strict
                   <<" ,  dRcut: "<<s_drcut
                   <<" ,  reducemaxenergy: "<<s_reducemaxenergy
                   <<" ,  increasemaxweight: "<<s_increasemaxweight
                   <<" ,  IR cut-off: "<<((int)s_mode>0?s_ircutoff:0)
                   <<" in frame "<<irframe<<" ("<<s_ircutoffframe<<")"
                   <<" ,  UV cut-off: "<<s_uvcutoff
                   <<" ,  1/alpha: "<<1./s_alpha_input
                   <<" ,  use running parameters "<<s_userunningparameters
                   <<" ,  FF recoil scheme: "<<s_ffrecscheme
                   <<" ,  FI recoil scheme: "<<s_firecscheme;
  }
  msg_Debugging()<<"\n}"<<std::endl;
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

void Photons::RegisterDefaults()
{
  Scoped_Settings s{ Settings::GetMainSettings()["YFS"] };
  s["MODE"].SetDefault(yfsmode::full).UseNoneReplacements();
  s["PHOTON_SPLITTER_MODE"].SetDefault(15);
  s["USE_ME"].SetDefault(1);
  s["IR_CUTOFF"].SetDefault(1E-3);
  s["UV_CUTOFF"].SetDefault(std::numeric_limits<double>::max());
  s["1/ALPHAQED"].SetDefault(0.);
  s["USE_RUNNING_PARAMETERS"].SetDefault(0);
  s["IR_CUTOFF_FRAME"].SetDefault("Multipole_CMS");
  s["MAXEM"].SetDefault(std::numeric_limits<int>::max());
  s["MINEM"].SetDefault(0);
  s["DRCUT"].SetDefault(std::numeric_limits<double>::max());
  s["STRICTNESS"].SetDefault(0);
  s["REDUCE_MAXIMUM_ENERGY"].SetDefault(1.);
  s["INCREASE_MAXIMUM_WEIGHT"].SetDefault(1.);
  s["CHECK_FIRST"].SetDefault(0);
  s["FF_RECOIL_SCHEME"].SetDefault(2);
  s["FI_RECOIL_SCHEME"].SetDefault(2);
  s["HISTO_BASE_NAME"].SetDefault("weights");
}

bool Photons::AddRadiation(Blob * blob)
{
  if (s_mode==yfsmode::off) return m_success=true;
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
    msg_Error()<<METHOD<<"(): Momentum not conserved after photon radiation: "
               <<blob->CheckMomentumConservation()<<std::endl;
    msg_Debugging()<<*blob<<std::endl;
    return m_success=false;
  }
  if (m_success && m_photonsadded && m_splitphotons) {
    m_success = m_photonsplitter.SplitPhotons(blob);
  }
  if (!blob->MomentumConserved()) {
    msg_Error()<<METHOD<<"(): Momentum not conserved after photon splitting: "
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

