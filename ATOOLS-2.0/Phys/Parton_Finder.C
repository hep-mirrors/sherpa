#include "Parton_Finder.H"

#ifdef PROFILE__all
#define PROFILE__Parton_Finder
#endif
#ifdef PROFILE__Parton_Finder
#include "prof.hh" 
#else
#define PROFILE_HERE
#endif

using namespace ATOOLS;

bool Parton_Finder::Test(const Particle *cur) 
{
  switch (m_criterion) {
  case pfc::color: if (cur->GetFlow(m_color.first)==m_color.second) return true;
  }
  return false;
}

Parton_Finder::Parton_Finder():
  m_criterion(pfc::color) {}

const Particle *Parton_Finder::FindConstConnectedForward(const Particle *start)
{
  if (m_tested.find(start)!=m_tested.end()) return NULL;
  m_tested.insert(start);
  if (!Test(start) ||
      m_excludeflavours.find(start->Flav().Kfcode())!=m_excludeflavours.end()) return NULL;
  Blob *decay=start->DecayBlob();
  if (decay!=NULL) {
    if (m_excludeblobs.find(decay->Type())!=m_excludeblobs.end()) return start;
    const Particle *stop=NULL;
    for (size_t i=0;i<(size_t)decay->NOutP();++i) {
      const Particle *next=decay->ConstOutParticle(i);
      if ((stop=FindConstConnectedForward(next))!=NULL) {
	m_end=stop;
	break;
      }
    }
    if (stop==NULL) {
      Turn();
      for (size_t i=0;i<(size_t)decay->NInP();++i) {
	const Particle *next=decay->ConstInParticle(i);
	if ((stop=FindConstConnectedBackward(next))!=NULL) {
	  m_end=stop;
	  break;
	}
      }
    }
  }
  else {
    m_end=start;
  }
  return m_end;
}

const Particle *Parton_Finder::FindConstConnectedBackward(const Particle *start)
{
  if (m_tested.find(start)!=m_tested.end()) return NULL;
  m_tested.insert(start);
  if (!Test(start) ||
      m_excludeflavours.find(start->Flav().Kfcode())!=m_excludeflavours.end()) return NULL;
  Blob *production=start->ProductionBlob();
  if (production!=NULL) {
    if (m_excludeblobs.find(production->Type())!=m_excludeblobs.end()) return start;
    const Particle *stop=NULL;
    for (size_t i=0;i<(size_t)production->NInP();++i) {
      const Particle *previous=production->ConstInParticle(i);
      if ((stop=FindConstConnectedBackward(previous))!=NULL) {
	m_end=stop;
	break;
      }
    }
    if (stop==NULL) {
      Turn();
      for (size_t i=0;i<(size_t)production->NOutP();++i) {
	const Particle *previous=production->ConstOutParticle(i);
	if ((stop=FindConstConnectedForward(previous))!=NULL) {
	  m_end=stop;
	  break;
	}
      }
    }
  }
  else {
    m_end=start;
  }
  return m_end;
}

const Particle *Parton_Finder::FindConstConnected(const Particle *start,bool forward)
{
  PROFILE_HERE;
  const Particle *stop=start;
  do {
    start=stop;
    if (forward) stop=FindConstConnectedForward(start);
    else stop=FindConstConnectedBackward(start);
    if (stop==NULL) break;
    forward=!forward;
  } while (stop!=NULL);
  m_end=start;
  return m_end;
}

void Parton_Finder::Clear()
{
  m_excludeblobs.clear();
  m_excludeflavours.clear();
  m_start=NULL;
  m_end=NULL;
}

void Parton_Finder::Turn()
{
  m_color.first=3-m_color.first;
}

Particle *Parton_Finder::FindConnected(const Particle *start,bool forward,unsigned int index)
{ 
  m_end=NULL; 
  m_tested.clear(); 
  if (index==0) index=1+(unsigned int)(start->Flav().IsAnti()^start->Flav().IsDiQuark());
  m_color.first=index;
  m_color.second=start->GetFlow(index);
  return (Particle*)FindConstConnected(start,forward); 
}

const Particle *Parton_Finder::FindConstConnected() 
{ 
  m_end=NULL; 
  m_tested.clear(); 
  return (Particle*)FindConstConnected(m_start); 
}

