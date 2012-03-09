#include "SINIC++/Cross_Sections/Cross_Sections.H"
#include "SINIC++/Tools/MinBias_Parameters.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Random.H"

using namespace SINIC;
using namespace ATOOLS;


Cross_Sections::
Cross_Sections(std::list<Omega_ik *> * eikonals,const double & energy,
	       const int & test) :
  p_eikonals(eikonals), p_selected(NULL), 
  m_originalY(MBpars("originalY")), m_cutoffY(MBpars("deltaY")), 
  m_Y(m_originalY-m_cutoffY),
  m_B(0.), m_Bmin(MBpars("bmin")), m_Bmax(MBpars("bmax")),  
  m_sigma_inelastic(Sigma_Inelastic(p_eikonals)), 
  m_sigma_elastic(Sigma_Elastic(p_eikonals,energy,test)),
  m_sigma_SD(Sigma_SD(&m_sigma_elastic,test)),
  m_sigma_DD(Sigma_DD(&m_sigma_elastic,&m_sigma_SD,test)),
  m_stot(0.), m_sinel(0.), m_sel(0.), m_sSD(0.), 
  m_test(test)
{  }

Cross_Sections::~Cross_Sections()
{ }

void Cross_Sections::CalculateTotalCrossSections()
{
  Sigma_Tot sigma_tot(p_eikonals,m_originalY,m_cutoffY,1.e-3);

  msg_Info()
    <<"==========================================================="<<std::endl
    <<"In "<<METHOD<<"(Y = "<<m_Y<<" from E = "
    <<(ATOOLS::Flavour(kf_p_plus).HadMass()*exp(m_originalY))<<")."<<std::endl;

  m_stot  = sigma_tot.Calculate(0.,m_Bmax);
  m_sinel = m_sigma_inelastic.Calculate(0.,m_Bmax);
  m_sel   = m_sigma_elastic.Calculate(0.,m_Bmax);
  m_sSD   = m_sigma_SD.Calculate(0.,m_Bmax);
  m_sDD   = m_sigma_DD.Calculate(0.,m_Bmax);

  Elastic_Slope slope(p_eikonals,m_originalY,m_cutoffY,1.e-3,m_stot);
  m_sval  = slope.Calculate(0.,m_Bmax+2.);
  msg_Info()
    <<"   sigma_tot = "<<m_stot/1.e9<<" mbarn (B = "<<m_sval<<"),"<<std::endl
    <<"   "
    <<"sigma_in ("<<m_sinel/1.e9<<"), "
    <<"sigma_el ("<<m_sel/1.e9<<"), "
    <<"sigma_SD ("<<m_sSD/1.e9<<"), "
    <<"sigma_DD ("<<m_sDD/1.e9<<") mbarn."<<std::endl;
  if (1.-dabs((m_sinel+m_sel+m_sSD+m_sDD)/m_stot)>1.e-2) {
    msg_Info()<<"   Sum = "<<(m_sinel+m_sel+m_sSD+m_sDD)/1.e9
	      <<" vs. "<<m_stot/1.e9<<" mbarn, "
	      <<"should maybe adjust sigma_inel to fit sigma_tot.\n";
  }
  if(m_test==1) {
    sigma_tot.TestTotalCrossSection();
    m_sigma_elastic.TestElasticCrossSection();
    m_sigma_inelastic.TestInelasticCrossSection();
  }
  m_modemap[0]  = m_sel/m_stot;
  m_modemap[1]  = m_sSD/m_stot;
  m_modemap[2]  = m_sDD/m_stot;
  m_modemap[10] = m_sinel/m_stot;
}

int Cross_Sections::SelectCollisionMode() {
  double random(ran->Get());
  for (std::map<int,double>::iterator miter=m_modemap.begin();
       miter!=m_modemap.end();miter++) {
    random -= miter->second;
    if (random<=0) return miter->first;
  }
  return -1;
}


