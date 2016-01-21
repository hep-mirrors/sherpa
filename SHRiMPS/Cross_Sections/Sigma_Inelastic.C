#include "SHRiMPS/Cross_Sections/Sigma_Inelastic.H"
#include "SHRiMPS/Tools/Special_Functions.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace ATOOLS;

Sigma_Inelastic::~Sigma_Inelastic() {
  m_xsecs.clear();
  for (std::map<Omega_ik *,std::vector<double> * >::iterator
	 bg=m_Bgrids.begin();bg!=m_Bgrids.end();bg++) {
    delete bg->second;
  }
  m_Bgrids.clear();
}

void Sigma_Inelastic::FillDifferentialGrids() {
  m_sigma = 0.;
  for (std::list<Omega_ik *>::iterator eikonal=p_eikonals->begin();
       eikonal!=p_eikonals->end(); eikonal++) {
    m_sigma += m_xsecs[(*eikonal)] = FillBGrid((*eikonal));
  }
  msg_Info()<<METHOD<<" yields effective inelastic cross section "
	    <<"sigma = "<<m_sigma/1.e9<<" mbarn.\n";
}

double Sigma_Inelastic::XSec(Omega_ik * eik) {
  if (eik==NULL) return m_sigma;
  if (m_xsecs.find(eik)==m_xsecs.end()) abort();
  return m_xsecs[eik];
}

const double Sigma_Inelastic::FillBGrid(Omega_ik * eikonal) {
  p_eikonal = eikonal;
  std::vector<double> * grid = new std::vector<double>;
  double deltaB(eikonal->DeltaB()), B(0.), sigma(0.), val1(0.), val2(0.);
  grid->push_back(0.);
  do {
    B     += deltaB;
    val2   = 2.*M_PI*B*GetValue(B);
    sigma += deltaB*(val1+val2)/2.;
    grid->push_back(sigma);
    val1   = val2;
  } while (B<MBpars.GetEikonalParameters().bmax);
  m_Bgrids[eikonal] = grid;
  return sigma*rpa->Picobarn();
}  
  
Omega_ik * Sigma_Inelastic::SelectEikonal() {
  p_eikonal = 0;
  while (p_eikonal==NULL) {
    double disc = ran->Get()*m_sigma;
    for (std::map<Omega_ik *,double>::iterator eikiter=m_xsecs.begin();
	 eikiter!=m_xsecs.end();eikiter++) {
      disc-=eikiter->second;
      if (disc<=1.e-12) {
	p_eikonal = eikiter->first;
	break;
      }
    }
  }
  return p_eikonal;
}

double Sigma_Inelastic::SelectB() {
  if (p_eikonal==0) {
    msg_Error()<<"Error in "<<METHOD<<": no eikonal selected.\n";
    return -1.;
  }
  std::vector<double> * grid = m_Bgrids[p_eikonal];  
  double deltaB(p_eikonal->DeltaB()), B(-1.);
  do {
    double random = ran->Get()*(*grid)[grid->size()-1];
    size_t bin(0);
    while (bin<grid->size()-1 && (random-(*grid)[bin]>=0)) bin++;
    if (bin>=grid->size()) continue;
    double inthigh((*grid)[bin]), intlow((*grid)[bin-1]);
    double Bhigh(bin*deltaB), Blow((bin-1)*deltaB);
    B  = (Blow*(random-intlow)+Bhigh*(inthigh-random))/(inthigh-intlow);
  } while (B<0.);
  return B;
}


double Sigma_Inelastic::GetValue(const double & B) { 
  return p_eikonal->Prefactor()*(1.-exp(-(*p_eikonal)(B))); 
}

double Sigma_Inelastic::GetCombinedValue(const double & B) { 
  double value(0.);
  for (std::list<Omega_ik *>::iterator eikonal=p_eikonals->begin();
       eikonal!=p_eikonals->end(); eikonal++) {
    value += (*eikonal)->Prefactor()*(1.-exp(-(**eikonal)(B))); 
  }
  return value;
}

double Sigma_Inelastic::Test() {
  const double EulerGamma= 0.577215664901532860606512090082 ;
  double a(MBpars.GetFFParameters().Lambda2/
	   (8.*(1.+MBpars.GetFFParameters().kappa)));
  double c(MBpars.GetEikonalParameters().beta02*
	   MBpars.GetFFParameters().Lambda2*
	   (1.+MBpars.GetFFParameters().kappa)*
	   exp(2.*MBpars.GetEikonalParameters().Delta*
	       MBpars.GetEikonalParameters().Ymax)/(8.*M_PI));
  double alpha(2.*M_PI*MBpars.GetFFParameters().norm);
  return alpha*(EulerGamma+log(c))/(2.*a)*rpa->Picobarn();
}
