#include "QCD_Remnant_Info.H"

#include <iomanip>

#ifdef PROFILE__all
#define PROFILE__QCD_Remnant_Info
#endif
#ifdef PROFILE__QCD_Remnant_Info
#include "prof.hh" 
#else
#define PROFILE_HERE
#endif

using namespace SHERPA;

std::set<QCD_Remnant_Info*> QCD_Remnant_Info::m_tested;

std::ostream &SHERPA::operator<<(std::ostream &ostr,QCD_Remnant_Info &info)
{
  ostr<<"m_this = ("<<info()<<") {\n";
  for (short unsigned int i=0;i<2;++i) ostr<<"   m_this["<<i<<"]      = "<<*info(i)<<"\n";
  for (short unsigned int i=0;i<2;++i) {
    if (info[i]!=NULL) ostr<<"   m_connected["<<i<<"] = "<<*(*info[i])(1-i)<<" <-> ("
			   <<std::setw(3)<<(info-i)<<" -> "<<std::setw(3)<<(info+i)<<")\n";
  }
  if (++info!=NULL) ostr<<"   m_finalstate   = "<<*++info<<"\n";
  return ostr<<"   m_anti         = "<<!info<<"\n}"<<std::flush;
}

QCD_Remnant_Info::QCD_Remnant_Info(ATOOLS::Particle *const particle): 
  p_finalstate(NULL),
  m_anti((unsigned int)(particle->Flav().IsAnti()^particle->Flav().IsDiQuark()))
{ 
  for (short unsigned int i=0;i<2;++i) { 
    p_connected[i]=NULL; 
    m_oldcolor[i]=0; 
  } 
  p_this[m_anti]=particle;
  SelectCompanion();
}

void QCD_Remnant_Info::SelectCompanion() 
{ 
  if (p_this[m_anti]->Flav().IsGluon()) {
    p_this[1-m_anti]=p_this[m_anti];
  }
  else {
    p_this[1-m_anti] = new ATOOLS::Particle(1,p_this[m_anti]->Flav().Bar());
    p_this[1-m_anti]->SetFlow(2-m_anti,ATOOLS::Flow::Counter());
  }
}

bool QCD_Remnant_Info::Find(const QCD_Remnant_Info *info,const bool forward) 
{
  if (this==info) return true;
  if (m_tested.find(this)!=m_tested.end()) return false;
  m_tested.insert(this);
  if ((*this)[forward]!=NULL) return (*this)[forward]->Find(info,forward);
  return false;
}

int QCD_Remnant_Info::Find(const QCD_Remnant_Info *info) 
{
  PROFILE_HERE;
  m_tested.clear();
  for (short unsigned int i=0;i<2;++i) if (Find(info,i)) return i;
  return -1;
}

bool QCD_Remnant_Info::Find(const ATOOLS::Particle *particle,const bool forward) 
{
  if (p_this[0]==particle || p_this[1]==particle) return true;
  if (m_tested.find(this)!=m_tested.end()) return false;
  m_tested.insert(this);
  //  msg_Debugging()<<"   "<<forward<<" "<<*p_this[m_anti]<<std::endl;
  if ((*this)[forward]!=NULL) return (*this)[forward]->Find(particle,forward);
  return false;
}

int QCD_Remnant_Info::Find(const ATOOLS::Particle *particle) 
{
  //  msg_Debugging()<<"Find "<<*particle<<" {"<<std::endl;
  PROFILE_HERE;
  for (short unsigned int i=0;i<2;++i) {
      m_tested.clear();
      if (Find(particle,i)) return i;
  }
  return -1;
}

