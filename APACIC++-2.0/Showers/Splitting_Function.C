#include "Splitting_Function.H"
#include "Message.H"

using namespace APACIC;


const double Splitting_Function::s_Nc = 3.;
const double Splitting_Function::s_CF = (s_Nc*s_Nc-1.)/(2.*s_Nc);
const double Splitting_Function::s_CA = s_Nc;
const double Splitting_Function::s_TR = 1./2.;

// some virtual functions

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
  ATOOLS::msg.Error()<<" Error in Spliting_Function: something nasty is going on"<<std::endl;
}

void Splitting_Function::SelectOne() {}

ATOOLS::Flavour & Splitting_Function::GetFlA() 
{ 
  return m_flavs[0];
}

ATOOLS::Flavour & Splitting_Function::GetFlB() 
{ 
  return m_flavs[1];
}

ATOOLS::Flavour & Splitting_Function::GetFlC() 
{ 
  return m_flavs[2];
}

void Splitting_Function::PrintStat(int mode) {
  if (!ATOOLS::msg.LevelIsDebugging()) return;
  if (mode>0) for(int i=0;i<mode;++i) ATOOLS::msg.Out()<<' ';
  ATOOLS::msg.Out()<<"Splitting Function: "
		   <<GetFlA()<<" -> "<<GetFlB()<<" + "<<GetFlC()<<std::endl;
}
