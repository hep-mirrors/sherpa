#include "ATOOLS/Phys/Parton_Finder.H"

#include "ATOOLS/Org/Exception.H"

size_t s_maxdepth=1000;

using namespace ATOOLS;

Parton_Tester::~Parton_Tester()
{
}

void Parton_Tester::Turn()
{
  THROW(fatal_error,"Virtual method called.");
}

bool Parton_Tester::Test(const Particle *parton) const
{
  THROW(fatal_error,"Virtual method called.");
  return false;
}

bool Parton_Finder::Test(const Particle *cur) 
{
  return p_criterion->Test(cur);
}

void Parton_Finder::Turn()
{
  p_criterion->Turn();
}

Parton_Finder::Parton_Finder(Parton_Tester &criterion):
  p_criterion(&criterion), m_forward(true) 
{
  m_excludeblobs.insert(btp::Signal_Process);
}

const Particle *Parton_Finder::
FindConstConnectedForward(const Particle *start)
{
  if (m_track.size()>s_maxdepth) {
    THROW(critical_error,"Caught in infinite loop.");
    return start;
  }
  if (!Test(start) ||
      m_excludeflavours.find(start->Flav().Kfcode())!=m_excludeflavours.end())
    return NULL;
  m_track.push_back((Particle*)start);
  Blob *decay=start->DecayBlob();
  if (decay==NULL) return m_end=start;
  if (m_excludeblobs.find(decay->Type())!=m_excludeblobs.end())
    return NULL;
  const Particle *stop=NULL;
  for (int i=decay->NOutP()-1;i>=0;--i) {
    const Particle *next=decay->ConstOutParticle(i);
    if (m_forward && next==m_track.front()) continue;
    if ((stop=FindConstConnectedForward(next))!=NULL) break;
  }
  if (stop==NULL) {
    Turn();
    for (int i=decay->NInP()-1;i>=0;--i) {
      const Particle *next=decay->ConstInParticle(i);
      if (next==start) continue;
      if (!m_forward && next==m_track.front()) continue;
      if ((stop=FindConstConnectedBackward(next))!=NULL) break;
    }
  }
  if (stop==NULL) stop=start;
  return m_end=stop;
}

const Particle *Parton_Finder::
FindConstConnectedBackward(const Particle *start)
{
  if (m_track.size()>s_maxdepth) {
    THROW(critical_error,"Caught in infinite loop.");
    return start;
  }
  if (!Test(start) ||
      m_excludeflavours.find(start->Flav().Kfcode())!=m_excludeflavours.end())
    return NULL;
  m_track.push_back((Particle*)start);
  Blob *production=start->ProductionBlob();
  if (production==NULL) return m_end=start;
  if (m_excludeblobs.find(production->Type())!=m_excludeblobs.end())
    return NULL;
  const Particle *stop=NULL;
  for (int i=production->NInP()-1;i>=0;--i) {
    const Particle *previous=production->ConstInParticle(i);
    if (!m_forward && previous==m_track.front()) continue;
    if ((stop=FindConstConnectedBackward(previous))!=NULL) break;
  }
  if (stop==NULL) {
    Turn();
    for (int i=production->NOutP()-1;i>=0;--i) {
      const Particle *previous=production->ConstOutParticle(i);
      if (previous==start) continue;
      if (m_forward && previous==m_track.front()) continue;
      if ((stop=FindConstConnectedForward(previous))!=NULL) break;
    }
  }
  if (stop==NULL) stop=start;
  return m_end=stop;
}

const Particle *Parton_Finder::
FindConstConnected(const Particle *start,bool forward)
{
  m_track.clear();
  m_forward=forward;
  for (short unsigned int i=0;i<2;++i) {
    if (forward) {
      if (FindConstConnectedForward(start)!=NULL) break;
    }
    else {
      if (FindConstConnectedBackward(start)!=NULL) break;
    }
    forward=!forward;
  }
  if (m_end==NULL) m_end=start;
  if (msg_LevelIsDebugging()) {
    msg_Out()<<"Parton_Finder::FindConstConnected(..): {\n"
	     <<"   "<<*start<<" -> ("<<forward<<")\n";
    for (size_t i=0;i<m_track.size();++i) msg_Out()<<"\n   "<<*m_track[i];
    msg_Out()<<"\n}"<<std::endl;
  }
  return m_end;
}

void Parton_Finder::Clear()
{
  m_excludeblobs.clear();
  m_excludeblobs.insert(btp::Signal_Process);
  m_excludeflavours.clear();
  m_start=NULL;
  m_end=NULL;
}

Particle *Parton_Finder::FindConnected(const Particle *start,
				       bool forward)
{ 
  return (Particle*)FindConstConnected(start,forward); 
}

const Particle *Parton_Finder::FindConstConnected() 
{ 
  return (Particle*)FindConstConnected(m_start); 
}

