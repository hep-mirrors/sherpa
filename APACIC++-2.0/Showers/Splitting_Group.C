#include "Splitting_Group.H"

using namespace APACIC;

Splitting_Group::Splitting_Group(Splitting_Function * spl):partsums(0) 
{
  if (spl) group.Append(spl);
  selected = spl;
}

Splitting_Group::~Splitting_Group() 
{
  if (partsums) delete [] partsums;
}

double Splitting_Group::CrudeInt(double _zmin, double _zmax) 
{
  if (!partsums) partsums = new double[group.GetLength()];
  lastint = 0;
  int i   = 0;
  for (SplFunIter iter(group);iter();++iter,++i)
    partsums[i] = lastint += iter()->CrudeInt(_zmin,_zmax);
  return lastint;
}        


void Splitting_Group::SelectOne() {
  double rr = lastint*AMATOOLS::ran.Get();
  int i;
  for (i=0;partsums[i]<rr;++i) {  
//     cout<<lastint<<" : "<<partsums[i]<<" : "<<rr<<endl; 
//     cout<<group[i]->GetFlA()<<" -> "<<group[i]->GetFlB()<<" + "<<group[i]->GetFlC()<<endl;
  };
//   cout<<" i = "<<i<<endl;
  selected=group[i];
}


void Splitting_Group::PrintStat(int mode) {
  if (mode>0) for(int i=0;i<mode;++i) std::cout<<' ';
  std::cout<<"Splitting Group:"<<GetFlA()<<" -> "<<GetFlB()<<" + "<<GetFlC()<<std::endl<<std::endl;
  for (SplFunIter iter(group);iter();++iter) {
    iter()->PrintStat(mode+4);
  }
    
}


double Splitting_Group::operator()(double z)                        
{ 
  return (*selected)(z);
}

double Splitting_Group::GetZ()
{ 
  return selected->GetZ();
}
             
double Splitting_Group::GetCoupling() 
{ 
  return selected->GetCoupling();
}

double Splitting_Group::GetCoupling(double t) 
{ 
  return selected->GetCoupling(t);
}

double Splitting_Group::GetWeight(double z,double pt2,bool masses)  
{ 
  return selected->GetWeight(z,pt2,masses);
}

APHYTOOLS::Flavour & Splitting_Group::GetFlA()                      
{ 
  return selected->GetFlA();
}

APHYTOOLS::Flavour & Splitting_Group::GetFlB() 
{ 
  return selected->GetFlB();
}

APHYTOOLS::Flavour & Splitting_Group::GetFlC() 
{ 
  return selected->GetFlC();
}  

void Splitting_Group::Add(Splitting_Function * spl) {
  group.Append(spl);
  selected = spl;
}
