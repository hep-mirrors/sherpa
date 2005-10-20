#include "Sudakov_Tools.H"
#include "Run_Parameter.H"
#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"

#include "Message.H"

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
    m_alphaQEDmax = (*p_aqed)(tmax*rpa.gen.RenormalizationScaleFactor());    
    double cutq2  = static_cast<Running_AlphaS*>(p_as)->CutQ2();
    if (.25*tmin*rpa.gen.RenormalizationScaleFactor()<cutq2) 
      m_alphaSmax=AlphaS(cutq2); 
    else m_alphaSmax  = AlphaS(0.25*tmin);    
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
    FixLambda2(sqr((Flavour(kf::Z)).Mass())); 
    Setscalefac(tmin);   
  }
  else {
    m_alphaQEDmax       = 1./128.;     
    m_alphaSmax         = 0.2;         
    m_beta0 = m_lambda2 = 0.;          
    m_scalefac          = 1.;          
  }
  Output();
}

void Sudakov_Tools::Output() {
  msg.Out()<<"Initialise Sudakov-Tools : "<<std::endl
	   <<"   beta0      = "<<m_beta0<<std::endl
	   <<"   lambda2    = "<<m_lambda2<<std::endl	
	   <<"   alphaS(MZ) = "<<CrudeAlphaS(sqr((Flavour(kf::Z)).Mass()))<<"  (estimated)"<<std::endl
	   <<"   alphaS(MZ) = "<<AlphaS(sqr((Flavour(kf::Z)).Mass()))<<"  (exact)"<<std::endl
	   <<"   scalefac   = "<<m_scalefac<<"."<<std::endl;
}

double Sudakov_Tools::CrudeAlphaS(const double t) const 
{
  return m_scalefac/
    (m_beta0*log(dabs(t)*rpa.gen.RenormalizationScaleFactor()/m_lambda2));
}

double Sudakov_Tools::AlphaS(double t) const 
{
  double tmp((*p_as)(dabs(t)*rpa.gen.RenormalizationScaleFactor()));
  //std::cout<<"   AlphaS("<<t<<") = "<<tmp<<" vs. "
  //	   <<m_alphaSmax<<"("<<(static_cast<Running_AlphaS*>(p_as)->CutQ2())<<") -> "
  //	   <<(*p_as)(static_cast<Running_AlphaS*>(p_as)->CutQ2())<<std::endl;
  return tmp;
}

double Sudakov_Tools::Alpha(double t) const 
{
  return (*p_aqed)(ATOOLS::dabs(t)*
		   ATOOLS::rpa.gen.RenormalizationScaleFactor());
}

double Sudakov_Tools::Nf(double t) const 
{
  return ((Running_AlphaS*)p_as)->
    Nf(dabs(t)*rpa.gen.RenormalizationScaleFactor());
}

void Sudakov_Tools::Setscalefac(double t0) 
{
  m_scalefac = AlphaS(dabs(t0))/CrudeAlphaS(dabs(t0));
}

void Sudakov_Tools::FixLambda2(double t) 
{ 
  t         = dabs(t*rpa.gen.RenormalizationScaleFactor());
  m_beta0   = ((Running_AlphaS*)p_as)->Beta0(t)/M_PI;  
  m_lambda2 = t*exp(-1./(m_beta0*AlphaS(t)));
}

