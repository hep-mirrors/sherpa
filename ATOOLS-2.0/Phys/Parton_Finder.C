#include "Parton_Finder.H"

#include "Exception.H"

#ifdef PROFILE__all
#define PROFILE__Parton_Finder
#endif
#ifdef PROFILE__Parton_Finder
#include "prof.hh" 
#else
#define PROFILE_HERE
#endif

using namespace ATOOLS;

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
  p_criterion(&criterion) {}

const Particle *Parton_Finder::
FindConstConnectedForward(const Particle *start)
{
  if (!Test(start) ||
      m_excludeflavours.find(start->Flav().Kfcode())!=m_excludeflavours.end())
    return NULL;
  m_track.push_back(start);
  Blob *decay=start->DecayBlob();
  if (decay==NULL) return m_end=start;
  if (m_excludeblobs.find(decay->Type())!=m_excludeblobs.end())
    return start;
  const Particle *stop=NULL;
  for (size_t i=0;i<(size_t)decay->NOutP();++i) {
    const Particle *next=decay->ConstOutParticle(i);
    if ((stop=FindConstConnectedForward(next))!=NULL) break;
  }
  if (stop==NULL) {
    Turn();
    for (size_t i=0;i<(size_t)decay->NInP();++i) {
      const Particle *next=decay->ConstInParticle(i);
      if (next==start) continue;
      if ((stop=FindConstConnectedBackward(next))!=NULL) break;
    }
  }
  return m_end=stop;
}

const Particle *Parton_Finder::
FindConstConnectedBackward(const Particle *start)
{
  if (!Test(start) ||
      m_excludeflavours.find(start->Flav().Kfcode())!=m_excludeflavours.end())
    return NULL;
  m_track.push_back(start);
  Blob *production=start->ProductionBlob();
  if (production==NULL) return m_end=start;
  if (m_excludeblobs.find(production->Type())!=m_excludeblobs.end())
    return start;
  const Particle *stop=NULL;
  for (size_t i=0;i<(size_t)production->NInP();++i) {
    const Particle *previous=production->ConstInParticle(i);
    if ((stop=FindConstConnectedBackward(previous))!=NULL) break;
  }
  if (stop==NULL) {
    Turn();
    for (size_t i=0;i<(size_t)production->NOutP();++i) {
      const Particle *previous=production->ConstOutParticle(i);
      if (previous==start) continue;
      if ((stop=FindConstConnectedForward(previous))!=NULL) break;
    }
  }
  return m_end=stop;
}

const Particle *Parton_Finder::
FindConstConnected(const Particle *start,bool forward)
{
  PROFILE_HERE;
  m_track.clear();
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
  if (msg.LevelIsDebugging()) {
    msg.Out()<<"Parton_Finder::FindConstConnected(..): {\n"
	     <<"   "<<*start<<" -> ("<<forward<<")\n";
    for (size_t i=0;i<m_track.size();++i) msg.Out()<<"\n   "<<*m_track[i];
    msg.Out()<<"\n}"<<std::endl;
  }
  return m_end;
}

void Parton_Finder::Clear()
{
  m_excludeblobs.clear();
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

