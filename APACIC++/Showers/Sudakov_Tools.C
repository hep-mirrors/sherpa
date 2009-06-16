#include "APACIC++/Showers/Sudakov_Tools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"

#include "ATOOLS/Org/Message.H"

using namespace APACIC;
using namespace MODEL;
using namespace ATOOLS;

Sudakov_Tools::Sudakov_Tools(MODEL::Model_Base * model) :
  p_as(model->GetScalarFunction(std::string("alpha_S"))),
  p_aqed(model->GetScalarFunction(std::string("alpha_QED")))
{ 
}

void Sudakov_Tools::CalculateMaxCouplings(const int scheme,
					  const double tmin,const double tmax)
{
  if (scheme>0) {
    m_alphaQEDmax = (*p_aqed)(tmax);    
    double cutq2  = static_cast<Running_AlphaS*>(p_as)->CutQ2();
    if (tmin<cutq2) 
      m_alphaSmax=AlphaS(cutq2); 
    else m_alphaSmax  = AlphaS(tmin);    
    /*
      std::ofstream was;
      was.open("astest.dat");
      for (double Q(0.01);Q<100.;Q*=1.05) {
      PRINT_INFO(Q<<" "<<(*p_as)(sqr(Q)));
      was<<Q<<" "<<(*p_as)(sqr(Q))<<"\n";
      }
      was.close();
      abort();
    */
    FixLambda2(sqr((Flavour(kf_Z)).Mass(true))); 
  }
  else {
    m_alphaQEDmax       = 1./128.;     
    m_alphaSmax         = 0.2;         
    m_beta0 = m_lambda2 = 0.;          
  }
  Output();
}

void Sudakov_Tools::Output() {
  msg_Tracking()<<"Initialise Sudakov-Tools : "<<std::endl
		<<"   beta0        = "<<m_beta0<<std::endl
		<<"   lambda2      = "<<m_lambda2<<std::endl	
		<<"   alphaS(MZ)   = "<<AlphaS(sqr((Flavour(kf_Z)).Mass(true)))<<"  (exact)"<<std::endl
		<<"   alphaQED(MZ) = "<<Alpha(sqr((Flavour(kf_Z)).Mass(true)))<<"  (exact)"<<std::endl
		<<"   alphaQED_max = "<<m_alphaQEDmax<<"  (exact)"<<std::endl;
}

double Sudakov_Tools::AlphaS(double t) const 
{
  return (*p_as)(dabs(t));
}

double Sudakov_Tools::Alpha(double t) const 
{
  return (*p_aqed)(ATOOLS::dabs(t));
}

double Sudakov_Tools::Nf(double t) const 
{
  return ((Running_AlphaS*)p_as)->Nf(dabs(t));
}

void Sudakov_Tools::FixLambda2(double t) 
{ 
  t         = dabs(t);
  m_beta0   = ((Running_AlphaS*)p_as)->Beta0(t)/M_PI;  
  m_lambda2 = t*exp(-1./(m_beta0*AlphaS(t)));
}

