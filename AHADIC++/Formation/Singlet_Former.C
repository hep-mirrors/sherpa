#include "AHADIC++/Formation/Singlet_Former.H"
#include "ATOOLS/Org/Message.H"
#include <cassert>

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Singlet_Former::Singlet_Former(list<Singlet *> * singlets) :
  p_singlets(singlets)
{}

Singlet_Former::~Singlet_Former() {}

void Singlet_Former::ExtractOutgoingCols(Blob * blob) {
  // Filter relevant - i.e. coloured - particles and put a copy of them
  // into the "active" list of coloured particles.  It must be a copy,
  // because they are purely internal, will decay, and their decay products
  // will also not show up in the event record.
  Particle * part(NULL);
  for (unsigned int i=0;i<blob->NInP();i++) {
    part = blob->InParticle(i); 
    if ((part->Status()!=part_status::active  &&
	 part->Status()!=part_status::fragmented) ||
	(part->GetFlow(1)==0 && part->GetFlow(2)==0)) continue;
    m_colparts.push_back(part);
  }
}

void Singlet_Former::FormSinglets() {
  // While there are particles in the active list, singlets must be
  // constructed.  This will be done recursively along the following lines:
  // - Take first particle of "active" list, make a new singlet list and
  //   put this particle into it.  At the same time, erase the first particle
  //   from the "active" list.
  // - Look for colour connected partners - col2 will link up with col1.
  //   When found, the colour-connected particle will be erased from
  //   "active" list and put into the singlet list
  // - The process stops when the ring is closed
  // I/we still have to think about possible errors and how to catch them ... . 
  while (!m_colparts.empty()) {
    p_singlets->push_back(MakeAnother());
  }
}

Singlet * Singlet_Former::MakeAnother() {
  Singlet * partlist = new Singlet();
  Particle * part = FindStart();
  partlist->push_back(new Proto_Particle(*part));
  if (part->Flav().IsQuark()) partlist->back()->SetLeading(true);
  unsigned int col1 = part->GetFlow(1);
  unsigned int col2 = part->GetFlow(2);
  while (col2!=col1) {
    for (list<Particle *>::iterator pliter=m_colparts.begin();
	 pliter!=m_colparts.end();pliter++) {
      part = *pliter;
      if (col1==part->GetFlow(2)) {
	m_colparts.erase(pliter);
	col1 = part->GetFlow(1);
	partlist->push_back(new Proto_Particle(*part));
	if (part->Flav().IsQuark()) partlist->back()->SetLeading(true);
	break;
      }
    }
  }
  return partlist;
}

Particle * Singlet_Former::FindStart() {
  Particle * part(NULL);
  // Look for quark
  for (list<Particle *>::iterator pliter=m_colparts.begin();
       pliter!=m_colparts.end();pliter++) {
    if ((*pliter)->GetFlow(1)!=0 && (*pliter)->GetFlow(2)==0) {
      part = (*pliter);
      m_colparts.erase(pliter);
      break;
    }
  }
  // no quark found - take gluon
  if (part==NULL) {
    part = m_colparts.front();
    m_colparts.pop_front();
  }
  return part;
}

bool Singlet_Former::Extract(Blob * blob) {
  // We construct colour singlets, i.e. closed colour structures, from all
  // coloured active particles entering the blob - these incoming partons
  // originate from either shower or hadron decays.  
  assert(m_colparts.empty());
  ExtractOutgoingCols(blob);
  FormSinglets();

  return true;
}
