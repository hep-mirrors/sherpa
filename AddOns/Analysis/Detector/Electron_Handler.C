#include "AddOns/Analysis/Detector/Electron_Handler.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(Electron_Handler_Getter,"Electron_Handler",
	       Analysis_Object,Argument_Matrix);

Analysis_Object *
Electron_Handler_Getter::operator()(const Argument_Matrix &parameters) const
{			
  if (parameters.size()<1) return NULL;
  bool     mode(false);

  Electron_Handler * ehandler = new Electron_Handler(parameters());
  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
    if (cur[0]=="Energymode") 
      ehandler->SetEModes(ATOOLS::ToType<int>(cur[1]),ATOOLS::ToType<int>(cur[2]));
    else if (cur[0]=="Directionmode") 
      ehandler->SetDModes(ATOOLS::ToType<int>(cur[1]),ATOOLS::ToType<int>(cur[2]));
    else if (cur[0]=="Miss_prob")
      ehandler->SetMissProb(ATOOLS::ToType<double>(cur[1]));
    else if (cur[0]=="Threshold")  
      ehandler->SetThreshold(ATOOLS::ToType<double>(cur[1]));
    else if (cur[0]=="E_fraction") 
      ehandler->SetEfrac(ATOOLS::ToType<double>(cur[1]));
    else if (cur[0]=="E_visfrac") 
      ehandler->SetEvisfracParams(ATOOLS::ToType<double>(cur[1]),ATOOLS::ToType<double>(cur[2]));
    else if (cur[0]=="SegmentECal") {
      mode = (cur[1]=="Energy");
      ehandler->FillSegmentParameters(mode,true,cur,1);
    }
    else if (cur[0]=="SegmentHCal") {
      mode = (cur[1]=="Energy");
      ehandler->FillSegmentParameters(mode,false,cur,1);
    }
  }
  return ehandler;
}

void Electron_Handler_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"Specifications [keyword value(s)]:\n"
     <<std::setw(width+7)<<" "<<"- Energymode    (0=no smearing,1=3 params form,2 Gaussian)\n"
     <<std::setw(width+7)<<" "<<"- Directionmode (0=no smearing,1=3 params form,2 Gaussian)\n"
     <<std::setw(width+7)<<" "<<"- Miss_prob (probability to miss an electron)\n"
     <<std::setw(width+7)<<" "<<"- Threshold (threshold for the electron to depose any energy)\n"
     <<std::setw(width+7)<<" "<<"- Efrac     (energy fraction in ECal)\n"
     <<std::setw(width+7)<<" "<<"Further specifications [keyword parameters]:\n"
     <<std::setw(width+7)<<" "<<"- Segment mode(Energy or Direction) eta-interval smearing_parameters.\n"
     <<std::setw(width+4)<<" "<<"}\n";
}


Electron_Handler::Electron_Handler(Primitive_Analysis * ana) : 
  Particle_Smearer_Base(ana,"ElectronSmearer")
{ 
  p_qualifier = Particle_Qualifier_Getter::GetObject("91","electron");
}

Electron_Handler::~Electron_Handler() {}

void Electron_Handler::DetermineTracker() {
  m_track = true;
  m_eta_Track = p_part->Momentum().Eta();
  m_phi_Track = p_part->Momentum().Phi();
}

void Electron_Handler::DetermineECal() {
  Vec4D mom(p_part->Momentum());
  double E(mom[0]-m_E_deposed), eta(mom.Eta()), phi(mom.Phi());
  if (E>m_threshold) {
    Deflect(true,m_Dmode_ECal,E,eta,phi);
    SmearEnergy(true,m_Emode_ECal,m_efrac*m_evis*E,eta);
  }
  else m_efrac = 1.;
}

void Electron_Handler::DetermineHCal() {
  if (m_efrac>=0.999) return;
  Vec4D mom(p_part->Momentum());
  double E(mom[0]-m_E_deposed), eta(mom.Eta()), phi(mom.Phi());
  Deflect(false,m_Dmode_HCal,(1.-m_efrac)*E,eta,phi);
  SmearEnergy(false,m_Emode_HCal,(1.-m_efrac)*m_evis*E,eta);
}


