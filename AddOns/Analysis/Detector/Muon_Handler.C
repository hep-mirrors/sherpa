#include "AddOns/Analysis/Detector/Muon_Handler.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(Muon_Handler_Getter,"Muon_Handler",
	       Analysis_Object,Argument_Matrix);

Analysis_Object *
Muon_Handler_Getter::operator()(const Argument_Matrix &parameters) const
{			
  if (parameters.size()<1) return NULL;

  Muon_Handler * mhandler = new Muon_Handler(parameters());
  mhandler->SetEModes(-1,-1);
  mhandler->SetDModes(0,0);
  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
    if (cur[0]=="Miss_prob")
      mhandler->SetMissProb(ATOOLS::ToType<double>(cur[1]));
    else if (cur[0]=="Threshold")  
      mhandler->SetThreshold(ATOOLS::ToType<double>(cur[1]));
    else if (cur[0]=="E_fraction") 
      mhandler->SetEfrac(ATOOLS::ToType<double>(cur[1]));
    else if (cur[0]=="E_visfrac") 
      mhandler->SetEvisfracParams(ATOOLS::ToType<double>(cur[1]),ATOOLS::ToType<double>(cur[2]));
    else if (cur[0]=="DepositECal") {
      mhandler->SetECalParams(ATOOLS::ToType<double>(cur[1]),ATOOLS::ToType<double>(cur[2]));
    }
    else if (cur[0]=="DepositHCal") {
      mhandler->SetHCalParams(ATOOLS::ToType<double>(cur[1]),ATOOLS::ToType<double>(cur[2]));
    }
  }
  return mhandler;
}

void Muon_Handler_Getter::
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

Muon_Handler::Muon_Handler(Primitive_Analysis * ana) : 
  Particle_Smearer_Base(ana,"MuonSmearer")
{ 
  p_qualifier = Particle_Qualifier_Getter::GetObject("92","muon");
}

Muon_Handler::~Muon_Handler() {
  m_Emode_ECal = m_Emode_HCal = -1;
}

void Muon_Handler::SetECalParams(const double mean,const double width) {
  std::vector<double > * params(new std::vector<double >);
  params->push_back(mean);
  params->push_back(width);
  etainterval limits;
  limits.first  = -100.;
  limits.second = +100.;
  m_Esmearingparams_ECal[limits] = params;
}

void Muon_Handler::SetHCalParams(const double mean,const double width) {
  std::vector<double > * params(new std::vector<double >);
  params->push_back(mean);
  params->push_back(width);
  etainterval limits;
  limits.first  = -100.;
  limits.second = +100.;
  m_Esmearingparams_HCal[limits] = params;
}


void Muon_Handler::DetermineTracker() {
  m_track = true;
  m_eta_Track = p_part->Momentum().Eta();
  m_phi_Track = p_part->Momentum().Phi();
}

void Muon_Handler::DetermineECal() {
  if (m_Emode_HCal<0) {
    m_E_ECal = 0.;
    return;
  }
  Vec4D mom(p_part->Momentum());
  double E(mom[0]-m_E_deposed), eta(mom.Eta()), phi(mom.Phi());
  Deflect(true,m_Dmode_ECal,E,eta,phi);
  SmearEnergy(true,m_Emode_ECal,E,eta);
}

void Muon_Handler::DetermineHCal() {
  if (m_Emode_HCal<0) {
    m_E_HCal = 0.;
    return;
  }
  Vec4D mom(p_part->Momentum());
  double E(mom[0]-m_E_deposed), eta(mom.Eta()), phi(mom.Phi());
  Deflect(false,m_Dmode_HCal,E,eta,phi);
  SmearEnergy(false,m_Emode_HCal,E,eta);
}

void Muon_Handler::DetermineMuonChambers() {
  m_mc = true;
  m_eta_MChamb = p_part->Momentum().Eta();
  m_phi_MChamb = p_part->Momentum().Phi();
}

