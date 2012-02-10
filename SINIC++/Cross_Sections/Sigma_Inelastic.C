#include "SINIC++/Cross_Sections/Sigma_Inelastic.H"
#include "SINIC++/Tools/Special_Functions.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace SINIC;
using namespace ATOOLS;

Sigma_Inelastic::~Sigma_Inelastic() {
  m_xsecs.clear();
  m_grids.clear();
}

double Sigma_Inelastic::FixEikonalAndImpact(Omega_ik *& eikonal) {
  double random(m_sigma*ATOOLS::ran->Get()*0.99999999999);
  for (std::map<Omega_ik *,double>::iterator xseciter=m_xsecs.begin();
       xseciter!=m_xsecs.end();xseciter++) {
    random -= xseciter->second;
    if (random<0.) {
      eikonal = xseciter->first;
      break;
    }
  }
  if (eikonal==NULL) {
    msg_Error()<<"Error in "<<METHOD<<": "<<std::endl
	       <<"   No eikonal selected, take the first one."<<std::endl;
    eikonal = m_xsecs.begin()->first;
  }
  if (m_grids.find(eikonal)==m_grids.end()) {
    msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
	       <<"   Did not find eikonal in grid-map, will exit."<<std::endl;
    exit(1);
  }
  std::vector<double> & grid = m_grids[eikonal];  
  unsigned int i;
  double B1, B2, B;
  do {
    random = ATOOLS::ran->Get()*0.99999999999;
    i = 0;
    while (i<grid.size()-1 && (random-grid[i]>=0)) i++;
    
    msg_Debugging()<<"In "<<METHOD<<"("<<random<<" --> "<<i<<")"<<std::endl;
    
    B1 = m_Bmin+(i-1)*m_deltaB;
    B2 = i==grid.size()-1?m_Bmax:B1+m_deltaB;
    B  = (B2*(random-grid[i-1])+B1*(grid[i]-random))/(grid[i]-grid[i-1]);
    
  } while (B>0.8*m_Bmax);
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

void Sigma_Inelastic::FillGrid(const double & Bmin,const double & Bmax,
			       const double & deltaB,const double & sigma) 
{
  m_Bmin   = Bmin;
  m_Bmax   = Bmax;
  m_deltaB = deltaB;
  
  double B,val1,val2,cumul(0.);

  std::vector<double> grid;
  grid.push_back(0.);
  B    = m_Bmin;
  val1 = 2.*M_PI*B*GetValue(B);
  while (B<m_Bmax) {
    B     += m_deltaB;
    val2   = 2.*M_PI*B*GetValue(B);
    cumul += m_deltaB*(val1+val2)/2.;
    grid.push_back(cumul);
    val1   = val2;
  }
  B = m_Bmin;
  for (size_t i=0;i<grid.size();i++) {
    grid[i] /= cumul;
    B       += m_deltaB;
  }

  m_grids[p_eikonal] = grid;
  m_xsecs[p_eikonal] = sigma;
}

void Sigma_Inelastic::SetSigma(const double & sigma) {
  if (sigma>=0) m_sigma = sigma;
  else {
    m_sigma = 0.;
    for (std::map<Omega_ik *,double>::iterator xseciter=m_xsecs.begin();
	 xseciter!=m_xsecs.end();xseciter++) {
      m_sigma += xseciter->second;
    }
  }
}


void Sigma_Inelastic::TestInelasticCrossSection(){
  const double EulerGamma= 0.577215664901532860606512090082 ;
  double m_a,m_c,m_alpha,m_res;
  double m_Delta,m_prefactor,m_Lambda2,m_beta0,m_kappa;
  m_Delta=(*p_eikonals).front()->Delta();
  m_prefactor=(*p_eikonals).front()->Prefactor();
  m_kappa=(*p_eikonals).front()->Kappa_i();
  m_Lambda2=(*p_eikonals).front()->Lambda2();
  m_beta0=(*p_eikonals).front()->FF1()->Beta0();
  m_a=m_Lambda2/(8.*(1.+m_kappa));
  m_c=ATOOLS::sqr(m_beta0)*m_Lambda2*(1.+m_kappa)*exp(2.*m_Delta*m_Y)/(8.*M_PI);
  m_alpha=2.*M_PI*m_prefactor;
  m_res=m_alpha*(EulerGamma+log(m_c))/(2.*m_a);
  msg_Out() << "In " << METHOD << " sigma_inelas = "<< m_res <<" 1/GeV^2 = "
	    <<m_res*rpa->Picobarn()/1.e9<<" mb ."<<std::endl;
}



