#include "Splitting_Function_Group.H"
#include "Random.H"

using namespace CS_SHOWER;
using namespace ATOOLS;
using namespace std;

ostream& CS_SHOWER::operator<<(std::ostream& str, Splitting_Function_Group &group) {
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
}


void Splitting_Function_Group::Add(Splitting_Function_Base * split) {
  m_splittings.push_back(split);
}


void Splitting_Function_Group::SelectOne() {
  double disc = m_lastint * ran.Get();
  m_splitter  = m_splittings.begin();
  do {
    disc -= (*m_splitter)->Last();
    if (disc<=0.) break;
    m_splitter++;
  } while (m_splitter!=m_splittings.end());
  if (m_splitter!=m_splittings.end()) p_selected = (*m_splitter);
  else {
    msg.Error()<<"Error in Splitting_Function_Group::SelectOne() : "<<endl
	       <<"   Try to select a splitting with "<<m_lastint<<" "<<disc<<endl;
    abort();
  }
}

double Splitting_Function_Group::OverIntegrated(const double zmin,const double zmax,
						const double scale) {
  m_lastint = 0.;
  for (m_splitter=m_splittings.begin();m_splitter!=m_splittings.end();m_splitter++) {
    m_lastint += (*m_splitter)->OverIntegrated(zmin,zmax,scale);
  }
  return m_lastint;
}


double Splitting_Function_Group::operator() (const double z,const double y) { 
  return (*p_selected)(z,y); 
}

double Splitting_Function_Group::Overestimated(const double z,const double y) { 
  return p_selected->Overestimated(z,y); 
}

double Splitting_Function_Group::RejectionWeight(const double z,const double y) { 
  return p_selected->RejectionWeight(z,y); 
}

double Splitting_Function_Group::Z() { 
  return p_selected->Z(); 
}         
