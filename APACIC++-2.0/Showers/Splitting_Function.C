#include "Splitting_Function.H"
//#include "Sudakov_Tools.H"

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

double Splitting_Function::GetLastInt() 
{
  return m_lastint;
}

void Splitting_Function::Add(Splitting_Function *) {
  std::cerr<<" Error in Spliting_Function: something nasty is going on"<<std::endl;
  std::cout<<" Error in Spliting_Function: something nasty is going on"<<std::endl;
}

void Splitting_Function::SelectOne() {}

APHYTOOLS::Flavour & Splitting_Function::GetFlA() 
{ 
  return m_flavs[0];
}

APHYTOOLS::Flavour & Splitting_Function::GetFlB() 
{ 
  return m_flavs[1];
}

APHYTOOLS::Flavour & Splitting_Function::GetFlC() 
{ 
  return m_flavs[2];
}

void Splitting_Function::PrintStat(int mode) {
  if (mode>0) for(int i=0;i<mode;++i) std::cout<<' ';
  std::cout<<"Splitting Function: "
	   <<GetFlA()<<" -> "<<GetFlB()<<" + "<<GetFlC()<<std::endl;
}
