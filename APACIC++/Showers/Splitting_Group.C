#include "APACIC++/Showers/Splitting_Group.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace APACIC;

Splitting_Group::Splitting_Group(ATOOLS::Mass_Selector *&ms,
				 Splitting_Function *const spl):
  Splitting_Function(ms), m_partsums(0) 
{
  if (spl!=NULL) m_splittings.push_back(spl);
  p_selected = spl;
}

Splitting_Group::~Splitting_Group() 
{
  while (m_splittings.size()) {
    delete m_splittings.back();
    m_splittings.pop_back();
  }
}

double Splitting_Group::CrudeInt(double zmin, double zmax) 
{
  if (m_partsums.empty()) m_partsums.resize(m_splittings.size());
  m_lastint=0.0;
  for (size_t size(m_splittings.size()), i(0);i<size;++i)
    m_partsums[i]=m_lastint+=m_splittings[i]->CrudeInt(zmin,zmax);
  return m_lastint;
}        

void Splitting_Group::SelectOne() 
{
  double rr(m_lastint*ATOOLS::ran.Get());
  size_t i(0);
  while (m_partsums[i++]<=rr);
  p_selected=m_splittings[--i];
}

double Splitting_Group::operator()(double z)                        
{ 
  return (*p_selected)(z);
}

double Splitting_Group::GetZ()
{ 
  return p_selected->GetZ();
}

double Splitting_Group::GetPhi(double z)
{
 return p_selected->GetPhi(z);
}

const ATOOLS::Simple_Polarisation_Info 
Splitting_Group::GetPolB(double z, double phi)
{
  return p_selected->GetPolB(z,phi);
}

const ATOOLS::Simple_Polarisation_Info 
Splitting_Group::GetPolC(double z, double phi, double phi_b)
{
  return p_selected->GetPolC(z,phi,phi_b);
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

const ATOOLS::Flavour & Splitting_Group::GetFlA() const
{ 
  return p_selected->GetFlA();
}

const ATOOLS::Flavour & Splitting_Group::GetFlB() const
{ 
  return p_selected->GetFlB();
}

const ATOOLS::Flavour & Splitting_Group::GetFlC() const
{ 
  return p_selected->GetFlC();
}  

void Splitting_Group::Add(Splitting_Function * spl) 
{
  m_splittings.push_back(spl);
  p_selected = spl;
}

double Splitting_Group::Integral(double zmin,double zmax)
{
  double value=0.; 
  for (Splitting_Vector::iterator sit(m_splittings.begin());
       sit!=m_splittings.end();++sit) 
    value+=(*sit)->Integral(zmin,zmax);
  return value;
}

void Splitting_Group::PrintStat(int mode) 
{
  if (mode>0) for(int i(0);i<mode;++i) msg_Debugging()<<' ';
  msg_Debugging()<<"Splitting Group: "<<GetFlA()<<" -> "
		 <<GetFlB()<<" + "<<GetFlC()<<std::endl<<std::endl;
  for (size_t i(0);i<m_splittings.size();++i) {
    m_splittings[i]->PrintStat(mode+4);
  }
}
