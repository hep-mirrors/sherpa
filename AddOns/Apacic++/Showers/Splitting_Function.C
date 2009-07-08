#include "AddOns/Apacic++/Showers/Splitting_Function.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace APACIC;


const double Splitting_Function::s_Nc = 3.;
const double Splitting_Function::s_CF = (s_Nc*s_Nc-1.)/(2.*s_Nc);
const double Splitting_Function::s_CA = s_Nc;
const double Splitting_Function::s_TR = 1./2.;
const double Splitting_Function::s_Tf = s_TR*3.0;
const double Splitting_Function::s_kappa  = s_CA*(67.0/18.0-ATOOLS::sqr(M_PI)/6.0) -
                                            s_Tf*10.0/9.0;
const double Splitting_Function::s_kappaG = s_CA*(67.0/18.0-ATOOLS::sqr(M_PI)/6.0);
const double Splitting_Function::s_kappaQ = s_TR*10.0/9.0;

int Splitting_Function::s_kfactorscheme = 0;

// some virtual functions

Splitting_Function::Splitting_Function(ATOOLS::Mass_Selector *&ms):
  p_ms(ms)
{
}

Splitting_Function::~Splitting_Function()
{
}

double Splitting_Function::GetPhi(double z) 
{ 
  return 2.*M_PI*ATOOLS::ran.Get(); 
}
const ATOOLS::Simple_Polarisation_Info Splitting_Function::GetPolB(double z, double phi) 
{ 
  return ATOOLS::Simple_Polarisation_Info(); 
}

const ATOOLS::Simple_Polarisation_Info Splitting_Function::GetPolC(double z, double phi, double phi_b) 
{ 
  return ATOOLS::Simple_Polarisation_Info();
}


double Splitting_Function::CrudeInt(double z0)  
{
  return CrudeInt(z0,1.-z0);
}

double Splitting_Function::Integral(double zmin,double zmax)  
{ 
  return 0.;
}

double Splitting_Function::GetLastInt() 
{
  return m_lastint;
}

void Splitting_Function::Add(Splitting_Function *) {
  msg_Error()<<" Error in Spliting_Function: something nasty is going on"<<std::endl;
}

void Splitting_Function::SelectOne() {}

const ATOOLS::Flavour &Splitting_Function::GetFlA() const 
{ 
  return m_flavs[0];
}

const ATOOLS::Flavour &Splitting_Function::GetFlB() const 
{ 
  return m_flavs[1];
}

const ATOOLS::Flavour &Splitting_Function::GetFlC() const
{ 
  return m_flavs[2];
}

void Splitting_Function::PrintStat(int mode) {
  if (!msg_LevelIsDebugging()) return;
  if (mode>0) for(int i=0;i<mode;++i) msg_Out()<<' ';
  msg_Out()<<"Splitting Function: "
		   <<GetFlA()<<" -> "<<GetFlB()<<" + "<<GetFlC()<<std::endl;
}
