#include "Splitting_Function.H"

using namespace APACIC;


const double Splitting_Function::Nc = 3.;
const double Splitting_Function::CF = (Nc*Nc-1.)/(2.*Nc);
const double Splitting_Function::CA = Nc;
const double Splitting_Function::TR = 1./2.;

// some virtual functions

double Splitting_Function::CrudeInt(double z0)  {return CrudeInt(z0,1.-z0);};
double Splitting_Function::GetLastInt() {return lastint;};

void Splitting_Function::Add(Splitting_Function *) {
  std::cerr<<" Error in Spliting_Function: something nasty is going on"<<std::endl;
  std::cout<<" Error in Spliting_Function: something nasty is going on"<<std::endl;
}

void Splitting_Function::SelectOne() {}

APHYTOOLS::Flavour & Splitting_Function::GetFlA() { return flavs[0];}
APHYTOOLS::Flavour & Splitting_Function::GetFlB() { return flavs[1];}
APHYTOOLS::Flavour & Splitting_Function::GetFlC() { return flavs[2];}

void Splitting_Function::PrintStat(int mode=0) {
  if (mode>0) for(int i=0;i<mode;++i) std::cout<<' ';
  std::cout<<"Splitting Function: "
	   <<GetFlA()<<" -> "<<GetFlB()<<" + "<<GetFlC()<<std::endl;
}
