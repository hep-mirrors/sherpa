#include "AddOns/Analysis/Detector/Particle_Smearer_Base.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "AddOns/Analysis/Detector/Detector.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/Random.H"

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__ROOT
#include "ATOOLS/Math/Scaling.H"
#include "TH2D.h"
#include "ATOOLS/Org/My_Root.H"
#endif 


using namespace ANALYSIS;
using namespace ATOOLS;

Particle_Smearer_Base::
Particle_Smearer_Base(Primitive_Analysis * ana,const std::string name) :
  m_Emode_ECal(0),m_Dmode_ECal(0),m_Emode_HCal(0),m_Dmode_HCal(0),
  m_missprob(0.),m_threshold(0.),m_evis_mean(1.),m_evis_slope(0.),m_efrac(1.),
  m_name(name),p_qualifier(NULL)
{
  p_ana = ana;
  Analysis_Object * det = p_ana->GetObject("Detector");
  if (det==NULL) {
    det = new Detector();
    p_ana->AddObject(det);
  }
  dynamic_cast<Detector *>(det)->AddParticleSmearer(this);
#ifdef USING__ROOT
  if (m_name=="HadronSmearer") {
    (*MYROOT::myroot)(new TH2D(ToString(this).c_str(),("smear_pion"),
			       50,0.,200.,100,-5.,5.),("smear_pion"));
  }
#endif
}

Particle_Smearer_Base::~Particle_Smearer_Base() {
  if (p_qualifier) { delete p_qualifier; p_qualifier=NULL; }
  for (m_spiter=m_Esmearingparams_ECal.begin();m_spiter!=m_Esmearingparams_ECal.end();m_spiter++) {
    if (m_spiter->second!=NULL) { delete m_spiter->second; m_spiter->second=NULL; }
  }
  m_Esmearingparams_ECal.clear();

  for (m_spiter=m_Dsmearingparams_ECal.begin();m_spiter!=m_Dsmearingparams_ECal.end();m_spiter++) {
    if (m_spiter->second!=NULL) { delete m_spiter->second; m_spiter->second=NULL; }
  }
  m_Dsmearingparams_ECal.clear();
  for (m_spiter=m_Esmearingparams_HCal.begin();m_spiter!=m_Esmearingparams_HCal.end();m_spiter++) {
    if (m_spiter->second!=NULL) { delete m_spiter->second; m_spiter->second=NULL; }
  }
  m_Esmearingparams_HCal.clear();

  for (m_spiter=m_Dsmearingparams_HCal.begin();m_spiter!=m_Dsmearingparams_HCal.end();m_spiter++) {
    if (m_spiter->second!=NULL) { delete m_spiter->second; m_spiter->second=NULL; }
  }
  m_Dsmearingparams_HCal.clear();
}

void Particle_Smearer_Base::
SetSmearingParams(etainterval interval,std::vector<double> * params,
		  std::map<etainterval, std::vector<double > *> * SPmap) 
{
  if (SPmap->find(interval)==SPmap->end()) (*SPmap)[interval] = params;
}

bool Particle_Smearer_Base::
FillSegmentParameters(const bool mode,const bool cal,
		      std::vector<std::string> cur,const size_t pos)
{
  std::pair<double,double> interval;
  std::vector<double>    * params = new std::vector<double>;
  size_t crit(0);

  int test = 0;
  if (mode) test = (cal)?m_Emode_ECal:m_Emode_HCal;
       else test = (cal)?m_Dmode_ECal:m_Dmode_HCal;
  if (test==1)      crit=3;
  else if (test==2) crit=2;

  interval.first  = ATOOLS::ToType<double>(cur[pos+1]);
  interval.second = ATOOLS::ToType<double>(cur[pos+2]);
  for (size_t i=3;i<crit+3;i++) params->push_back(ATOOLS::ToType<double>(cur[pos+i]));

  std::map<etainterval, std::vector<double > *> * SPmap;
  if (mode) SPmap = (cal)?(&m_Esmearingparams_ECal):(&m_Esmearingparams_HCal); 
       else SPmap = (cal)?(&m_Dsmearingparams_ECal):(&m_Dsmearingparams_HCal);  
  (*SPmap)[interval] = params;
  return true;
}

std::vector<double > * Particle_Smearer_Base::
GetSmearingParameters(std::map<etainterval, std::vector<double > *> * SPmap,const double eta) {
  std::vector<double > * params(NULL);
  for (std::map<etainterval, std::vector<double > *>::iterator piter=SPmap->begin();
       piter!=SPmap->end();piter++) {
    if (piter->first.first<=eta&&eta<=piter->first.second) {
      params=piter->second; 
      break;
    }
  }
  return params;
}

void Particle_Smearer_Base::CalculateEvis(const double E) {
  double rana,dummy;
  m_evis = -1.;
  do {
    ran.Gaussian(rana,dummy);
    m_evis = m_evis_mean - m_evis_slope*rana;
  } while (m_evis<0. || m_evis>1.);
}

bool Particle_Smearer_Base::TreatParticle(Particle * part) {
  if (!(*p_qualifier)(part)) return false;
  m_E_deposed = m_eta_ECal = m_phi_ECal = m_eta_HCal = m_phi_HCal = m_E_ECal = m_E_HCal = 0.;
  m_eta_Track = m_phi_Track = m_eta_MChamb = m_phi_MChamb = 0.;
  m_track = m_mc = false;
  if (m_missprob>ran.Get()) return true;
  p_part = part;
  CalculateEvis(part->Momentum()[0]);
  DetermineTracker();
  DetermineECal();
  DetermineHCal();
  DetermineMuonChambers();

  //std::cout<<METHOD<<" for "<<part->Flav()<<"  "<<part->Momentum()<<" : "<<std::endl
  //	   <<"        "
  //	   <<" ECal = "<<m_E_ECal<<" ("<<m_eta_ECal<<","<<m_phi_ECal<<" ),"
  //	   <<" HCal = "<<m_E_HCal<<" ("<<m_eta_HCal<<","<<m_phi_HCal<<" )."<<std::endl;
  return true;
}

bool Particle_Smearer_Base::GivesTrack(double & eta,double & phi) {
  if (!m_track) return false;
  eta = m_eta_Track; 
  phi = m_phi_Track; 
  return true;
}

double Particle_Smearer_Base::EnergyInECal(double & eta,double & phi) {
  eta = m_eta_ECal; 
  phi = m_phi_ECal; 
  return m_E_ECal;
}

double Particle_Smearer_Base::EnergyInHCal(double & eta,double & phi) {
  eta = m_eta_HCal; 
  phi = m_phi_HCal; 
  return m_E_HCal;
}

bool Particle_Smearer_Base::GivesMuon(double & eta,double & phi) {
  if (!m_mc) return false;
  eta = m_eta_MChamb; 
  phi = m_phi_MChamb; 
  return true;
}

void Particle_Smearer_Base::Deflect(const bool cal,const int mode,
				    const double E,const double eta,const double phi) 
{
  double sigma(0.),ct((exp(eta)-exp(-eta))/(exp(eta)+exp(-eta))),angle(0.);
  std::vector<double> * params = 
    GetSmearingParameters((cal)?(&m_Dsmearingparams_ECal):(&m_Dsmearingparams_HCal),eta);
  if (params) {
    double rana,ranb,ranc,dummy;
    angle = 2.*M_PI*ran.Get();
    switch (mode) {
    case 1:
      ran.Gaussian(rana,ranb);
      ran.Gaussian(ranc,dummy);
      sigma = (*params)[0]*rana + (*params)[1]/sqrt(E)*ranb + (*params)[0]/E*ranc; 
      break;
    case 0:
    default:
      break;
    }
    ct   = cos(acos(ct)+sigma*cos(angle));
  }
  if (cal) { 
    m_eta_ECal = 0.5*log((1+ct)/(1-ct)); 
    m_phi_ECal = phi+sigma*sin(angle); 
  }
  else { 
    m_eta_HCal = 0.5*log((1+ct)/(1-ct)); 
    m_phi_HCal = phi+sigma*sin(angle); 
  }
}

void Particle_Smearer_Base::SmearEnergy(const bool cal,const int mode,
					const double E,const double eta)
{
  double sigma(1.),dep(E);
  std::vector<double> * params = 
    GetSmearingParameters((cal)?(&m_Esmearingparams_ECal):(&m_Esmearingparams_HCal),eta);
  if (params) {
    double rana,ranb,ranc,dummy;
    switch (mode) {
    case 1:
      ran.Gaussian(rana,ranb);
      ran.Gaussian(ranc,dummy);
      sigma = (*params)[0]*rana + (*params)[1]/sqrt(E)*ranb + (*params)[2]/E*ranc; 
      dep   = E*(1.+sigma);
      break;
    case 2:
      ran.Gaussian(rana,dummy);
      dep = (*params)[0] + (*params)[1]*rana;
      break;
    case 0:
    default:
      break;
    }
  }
  if (dep<0.) dep=0.;

  m_E_deposed += dep;
  if (cal) m_E_ECal = dep; 
  else {
    m_E_HCal = dep;
#ifdef USING__ROOT
    double logdiff = log(dabs(m_E_HCal-E)/E);
    ((TH2D*)(*MYROOT::myroot)[("smear_pion")])->Fill(E,logdiff,1.);
#endif
  }
}
