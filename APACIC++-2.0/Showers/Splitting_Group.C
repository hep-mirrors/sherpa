#include "Splitting_Group.H"
#include "Random.H"
#include "Message.H"

using namespace APACIC;

Splitting_Group::Splitting_Group(Splitting_Function * spl):p_partsums(0) 
{
  if (spl) m_group.Append(spl);
  p_selected = spl;
}

Splitting_Group::~Splitting_Group() 
{
  if (p_partsums) delete [] p_partsums;
}

double Splitting_Group::CrudeInt(double _zmin, double _zmax) 
{
  if (!p_partsums) p_partsums = new double[m_group.GetLength()];
  m_lastint = 0;
  int i   = 0;
  for (SplFunIter iter(m_group);iter();++iter,++i)
    p_partsums[i] = m_lastint += iter()->CrudeInt(_zmin,_zmax);
  return m_lastint;
}        


void Splitting_Group::SelectOne() {
  double rr = m_lastint*ATOOLS::ran.Get();
  int i;
  for (i=0;p_partsums[i]<rr;++i) {  
  }
  p_selected = m_group[i];
}


void Splitting_Group::PrintStat(int mode) {
  if (mode>0) for(int i=0;i<mode;++i) msg_Debugging()<<' ';
  msg_Debugging()<<"Splitting Group:"<<GetFlA()<<" -> "<<GetFlB()<<" + "<<GetFlC()<<std::endl<<std::endl;
  for (SplFunIter iter(m_group);iter();++iter) {
    iter()->PrintStat(mode+4);
  }
    
}


double Splitting_Group::operator()(double z)                        
{ 
  return (*p_selected)(z);
}

double Splitting_Group::GetZ()
{ 
  return p_selected->GetZ();
}
             
double Splitting_Group::GetCoupling() 
{ 
  return p_selected->GetCoupling();
}

double Splitting_Group::GetCoupling(double t) 
{ 
  return p_selected->GetCoupling(t);
}

double Splitting_Group::GetWeight(double z,double pt2,bool masses)  
{ 
  return p_selected->GetWeight(z,pt2,masses);
}

ATOOLS::Flavour & Splitting_Group::GetFlA()                      
{ 
  return p_selected->GetFlA();
}

ATOOLS::Flavour & Splitting_Group::GetFlB() 
{ 
  return p_selected->GetFlB();
}

ATOOLS::Flavour & Splitting_Group::GetFlC() 
{ 
  return p_selected->GetFlC();
}  

void Splitting_Group::Add(Splitting_Function * spl) 
{
  m_group.Append(spl);
  p_selected = spl;
}

double Splitting_Group::Integral(double zmin,double zmax)
{
  double value=0.; 
  for (SplFunIter iter(m_group);iter();++iter) {
    value+=iter()->Integral(zmin,zmax);
  }
  return value;
}
