#include "AddOns/Analysis/Detector/Hadron_Handler.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"

#ifdef USING__ROOT
#include "ATOOLS/Math/Scaling.H"
#include "TH1D.h"
#include "TH2D.h"
#endif 

using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(Hadron_Handler_Getter,"Hadron_Handler",
	       Analysis_Object,Argument_Matrix);

Analysis_Object *
Hadron_Handler_Getter::operator()(const Argument_Matrix &parameters) const
{			
  if (parameters.size()<1) return NULL;

  bool     mode(false);

  Hadron_Handler * hhandler = new Hadron_Handler(parameters());
  hhandler->SetEModes(1,1);
  hhandler->SetDModes(0,0);
  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
     if (cur[0]=="Miss_prob")
      hhandler->SetMissProb(ATOOLS::ToType<double>(cur[1]));
    else if (cur[0]=="ECal_Threshold")  
      hhandler->SetECalThreshold(ATOOLS::ToType<double>(cur[1]));
    else if (cur[0]=="HCal_Threshold")  
      hhandler->SetThreshold(ATOOLS::ToType<double>(cur[1]));
    else if (cur[0]=="E_fraction") 
      hhandler->SetEfrac(ATOOLS::ToType<double>(cur[1]));
    else if (cur[0]=="Ecrit") 
      hhandler->SetCriticals(ATOOLS::ToType<double>(cur[1]),ATOOLS::ToType<double>(cur[2]),
			     ATOOLS::ToType<double>(cur[3]));
    else if (cur[0]=="E_splitting") 
      hhandler->SetSplitting(ATOOLS::ToType<double>(cur[1]),ATOOLS::ToType<double>(cur[2]));
    else if (cur[0]=="E_visfrac") 
      hhandler->SetEvisfracParams(ATOOLS::ToType<double>(cur[1]),ATOOLS::ToType<double>(cur[2]));
    else if (cur[0]=="DepositECal") {
      hhandler->SetECalParams(ATOOLS::ToType<double>(cur[1]),ATOOLS::ToType<double>(cur[2]));
    }
    else if (cur[0]=="DepositHCal") {
      hhandler->SetHCalParams(ATOOLS::ToType<double>(cur[1]),ATOOLS::ToType<double>(cur[2]));
    }
    else if (cur[0]=="Energymode") 
      hhandler->SetEModes(ATOOLS::ToType<int>(cur[1]),ATOOLS::ToType<int>(cur[2]));
    else if (cur[0]=="Directionmode") 
      hhandler->SetDModes(ATOOLS::ToType<int>(cur[1]),ATOOLS::ToType<int>(cur[2]));
    else if (cur[0]=="SegmentECal") {
      mode = (cur[1]=="Energy");
      hhandler->FillSegmentParameters(mode,true,cur,1);
    }
    else if (cur[0]=="SegmentHCal") {
      mode = (cur[1]=="Energy");
      hhandler->FillSegmentParameters(mode,false,cur,1);
    }
  }
  return hhandler;
}

void Hadron_Handler_Getter::
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

Hadron_Handler::Hadron_Handler(Primitive_Analysis * ana) :
  Particle_Smearer_Base(ana,"HadronSmearer"),   
  m_Ecrit(15.),m_exponent(1.),m_threshold_ECal(.3),
  m_dep_mean_ECal(0.5),m_dep_width_ECal(0.1),m_dep_mean_HCal(0.5),m_dep_width_HCal(0.1)
{
  p_qualifier = Particle_Qualifier_Getter::GetObject("3","hadron");
#ifdef USING__ROOT
  std::string name = std::string("rEcal_had");
  (*MYROOT::myroot)(new TH2D(name.c_str(),name.c_str(),100,0.,200.,100,0.,2.0),name);
  name = std::string("rHcal_had");
  (*MYROOT::myroot)(new TH2D(name.c_str(),name.c_str(),100,0.,200.,100,0.,2.0),name);
#endif
}

Hadron_Handler::~Hadron_Handler() {}

void Hadron_Handler::SetECalParams(const double mean,const double width) {
  m_dep_mean_ECal  = mean; 
  m_dep_width_ECal = width;
}

void Hadron_Handler::SetCriticals(const double crit,const double exponent,const double ampl) {
  m_Ecrit    = crit;
  m_exponent = exponent;
  m_ampl     = ampl;
}

void Hadron_Handler::SetHCalParams(const double mean,const double width) {
  m_dep_mean_HCal  = mean; 
  m_dep_width_HCal = width;
}

void Hadron_Handler::SetSplitting(const double kappa, const double lambda) {
  m_kappa = kappa; 
  m_lambda = lambda;
}

void Hadron_Handler::DetermineTracker() {
  if (p_part->Flav().Charge()==0) return;
  m_track = true;
  m_eta_Track = p_part->Momentum().Eta();
  m_phi_Track = p_part->Momentum().Phi();
}

void Hadron_Handler::DetermineECal() {
  Vec4D mom(p_part->Momentum());
  double E(mom[0]-m_E_deposed), eta(mom.Eta()), phi(mom.Phi());
  CalculateDeposits(E,eta);
  E *= m_evis;
  int mode=(m_ECalfrac<=1.e-3)?2:m_Emode_ECal;

  double dep(E*m_Emode_ECal), norm(E/m_evis);
  if (mode==1) {
    std::vector<double> * params = GetSmearingParameters((&m_Esmearingparams_ECal),eta);
    if (params) {
      E    *= m_ECalfrac;
      norm *= m_ECalfrac;
      double rana,ranb,ranc,dummy,sigma;
      do {
	ran.Gaussian(rana,ranb);
	ran.Gaussian(ranc,dummy);
	sigma = (*params)[0]*rana + (*params)[1]/sqrt(E)*ranb + (*params)[0]/E*ranc; 
      } while (sigma<-1.);
      dep   = E*(1.+sigma);
    }
    else mode=2;
  }
  if (mode==2) {
    double rana,ranb;
    do {
      ran.Gaussian(rana,ranb);
      if (m_dep_mean_ECal>E) dep = E+m_dep_width_ECal*E/m_dep_mean_ECal*ranb;
      else dep = m_dep_mean_ECal + m_dep_width_ECal*rana;
    } while (dep<0.);
  }

  m_E_ECal = dep;
  m_eta_ECal = eta; 
  m_phi_ECal = phi;
  m_E_deposed += dep;
#ifdef USING__ROOT
  std::string name = std::string("rEcal_had");
  ((TH2D*)(*MYROOT::myroot)[name])->Fill(norm,m_E_ECal/norm,1.);
#endif
}

void Hadron_Handler::DetermineHCal() {
  Vec4D mom(p_part->Momentum());
  double E(m_evis*mom[0]-m_E_deposed), eta(mom.Eta()), phi(mom.Phi());
  if (E<0.) return; 
  if (E<m_threshold) { m_E_deposed += m_E_HCal = E*ran.Get(); }
  else {
    Deflect(false,m_Dmode_HCal,E,eta,phi);
    SmearEnergy(false,m_Emode_HCal,E,eta);
  }
  m_eta_HCal = eta; 
  m_phi_HCal = phi;
#ifdef USING__ROOT
  std::string name = std::string("rHcal_had");
  ((TH2D*)(*MYROOT::myroot)[name])->Fill(E/m_evis,m_E_HCal*m_evis/E,1.);
#endif
}

void Hadron_Handler::CalculateEvis(const double E) {
  m_evis = 1.-(1.-m_evis_mean)*exp(m_evis_slope/pow(E,0.13));
  double rana,ranb,sigma;
  do {
    ran.Gaussian(rana,ranb);
    sigma = 0.007*rana+0.15/sqrt(E)*ranb;
  } while (sigma<-1. || sigma>1);
  m_evis *= (1.+sigma);
}

/* 
 * Hadrons act at least as minimally ionizing particles (MIPs), deposing energy with
 * a fixed mean of dep_mean which is Gaussian smeared with a width of dep_width.
 * Typically these parameters are around 0.5 GeV and 0.1 GeV.
 * On the other hand, they may depose up to the maximal energy fraction (of the 
 * order of 0.9) in the ECal.
 */

void Hadron_Handler::CalculateDeposits(const double E,const double eta) {
  double crit(pow(dabs(m_Ecrit-E)/m_Ecrit,m_exponent)),test,rana,dummy;
  if (E<m_Ecrit) crit = -crit;
  test = m_ampl/2.-atan(-crit)/M_PI;
  if (test<m_threshold_ECal/E) m_ECalfrac = -1.;
  else {
    do ran.Gaussian(rana,dummy); while (rana<0.);
    if (test>ran.Get()) m_ECalfrac = 1.-exp(-test)*rana;
    else m_ECalfrac = 0.+exp(-dabs(crit))*rana;
  }
  return;
}

