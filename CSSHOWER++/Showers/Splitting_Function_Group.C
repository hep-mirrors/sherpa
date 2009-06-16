#include "CSSHOWER++/Showers/Splitting_Function_Group.H"
#include "ATOOLS/Math/Random.H"

using namespace CSSHOWER;
using namespace ATOOLS;
using namespace std;

ostream& CSSHOWER::operator<<(std::ostream& str, Splitting_Function_Group &group) {
  str<<"Splitting_Function_Group : "<<group.m_lastint<<endl;
  for (std::list<Splitting_Function_Base *>::iterator splitter=group.m_splittings.begin();
       splitter!=group.m_splittings.end();splitter++) {
    str<<(**splitter);
  }
  str<<"-------------------------------------------------------------"<<endl;
  return str;
}

Splitting_Function_Group::~Splitting_Function_Group() {
  m_splitter=m_splittings.begin();
  do {
    if (*m_splitter) { delete (*m_splitter); (*m_splitter=NULL); }
    m_splitter = m_splittings.erase(m_splitter);
  } while (m_splitter!=m_splittings.end());
  m_splittings.clear();
}


void Splitting_Function_Group::Add(Splitting_Function_Base * split) {
  m_splittings.push_back(split);
}


void Splitting_Function_Group::SelectOne() {
  double disc = m_lastint;
  double random;
  bool valid=false;
  while (!valid) {
    random = ran.Get();
    if (random!=0.0) valid=true;
  }
  disc *= random;
  m_splitter  = m_splittings.begin();
  do {
    disc -= (*m_splitter)->Last();
    if (disc<=0.) break;
    m_splitter++;
  } while (m_splitter!=m_splittings.end());
  if (m_splitter!=m_splittings.end()) p_selected = (*m_splitter);
  else {
    msg_Error()<<"Error in Splitting_Function_Group::SelectOne() : "<<endl
	       <<"   Try to select a splitting with "<<m_lastint<<" "<<disc<<endl;
    abort();
  }
}

double Splitting_Function_Group::OverIntegrated(const double zmin,const double zmax,
						const double scale,const double xbj) {
  for (m_splitter=m_splittings.begin();m_splitter!=m_splittings.end();m_splitter++) 
    m_lastint += (*m_splitter)->OverIntegrated(zmin,zmax,scale,xbj); 
  return m_lastint;
}


double Splitting_Function_Group::operator() (const double z,const double y,
					     const double eta = 1.,const double scale=0.,
					     const double Q2=0.,int mode) { 
  return (*p_selected)(z,y,eta,scale,Q2,mode); 
}

double Splitting_Function_Group::Overestimated(const double z,const double y) { 
  return p_selected->Overestimated(z,y); 
}

double Splitting_Function_Group::RejectionWeight(const double z,const double y,
						 const double eta,const double scale,const double Q2) { 
  return p_selected->RejectionWeight(z,y,eta,scale,Q2); 
}

double Splitting_Function_Group::Z() { 
  return p_selected->Z(); 
}         

void Splitting_Function_Group::ClearSpecs()
{
  m_specs.clear();
  for (m_splitter=m_splittings.begin();m_splitter!=m_splittings.end();m_splitter++) 
    (*m_splitter)->ClearSpecs();
}

void Splitting_Function_Group::ResetLastInt()
{
  m_lastint=0.0;
  for (m_splitter=m_splittings.begin();m_splitter!=m_splittings.end();m_splitter++) 
    (*m_splitter)->ResetLastInt();
}
