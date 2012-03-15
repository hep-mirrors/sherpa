#include "SHRiMPS/Beam_Remnants/Singlet_Sorter.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;

Singlet_Sorter::Singlet_Sorter() {}

void Singlet_Sorter::Sort(PartList * inlist,PSetYSet * outlists) {
  m_inlist = inlist;
  unsigned int ref1, ref2, cont;
  Particle * part;
  PartYSet  * pset;
  while (!m_inlist->empty()) {
    part = FindNextStart();
    if (!part) {
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   No new particle found in non-empty list.\n";
      exit(1);
    } 
    pset = new PartYSet;
    pset->insert(part);
    ref1 = part->GetFlow(1);
    ref2 = part->GetFlow(2);
    //if (part) msg_Out()<<"found a start: ["<<ref1<<", "<<ref2<<"].\n";
    while (part) {
      part = FindNext(ref1,ref2);
      if (part) {
	ref1 = part->GetFlow(1);
	ref2 = part->GetFlow(2);
	//msg_Out()<<"   continue with: "<<ref1<<", "<<ref2<<".\n";
	pset->insert(part);
	if (part->GetFlow(1)==cont) break;
      }
    }
    outlists->insert(pset);
  }
}

void Singlet_Sorter::Sort(PartList * inlist,PCList * outlist) {
  m_inlist = inlist;
  unsigned int ref, ref2, cont;
  Particle * part;
  while (!m_inlist->empty()) {
    part = FindNextStart();
    if (!part) {
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   No new particle found in non-empty list.\n";
      exit(1);
    } 
    ref  = part->GetFlow(1);
    cont = part->GetFlow(2);
    outlist->push_back(make_pair(part,make_pair(ref,cont)));
    //if (part) msg_Out()<<"found a start: "<<ref<<".\n";
    while (part) {
      part = FindNext(ref);
      if (part) {
	ref  = part->GetFlow(1);
	ref2 = part->GetFlow(2);
	//msg_Out()<<"   continue with: "<<ref<<".\n";
	outlist->push_back(make_pair(part,make_pair(ref,ref2)));
	if (part->GetFlow(1)==cont) break;
      }
    }
  }
}

Particle * Singlet_Sorter::FindNextStart() {
  double maxFB(0.), y;
  PartList::iterator pit(m_inlist->end()),start(pit);
  for (pit=m_inlist->begin();pit!=m_inlist->end();pit++) {
    if ((*pit)->Flav().IsGluon() ||
	(*pit)->GetFlow(1)==0) continue;
    y = dabs((*pit)->Momentum().Y());
    if (y>maxFB) { 
      maxFB = y;
      start = pit;
    }
  }
  if (start==m_inlist->end()) {
    start = m_inlist->begin();
    for (pit=m_inlist->begin();pit!=m_inlist->end();pit++) {
      y = dabs((*pit)->Momentum().Y());
      if (y>maxFB) { 
	maxFB = y;
	start = pit;
      }
    }
  }
  if (start==m_inlist->end()) exit(1);
  m_inlist->erase(start);
  return (*start);
}

Particle * Singlet_Sorter::
FindNext(const unsigned int & ref1,const unsigned int & ref2) {
  Particle * part(NULL);
  for (PartList::iterator pit=m_inlist->begin();pit!=m_inlist->end();pit++) {
    part = (*pit);
    //msg_Out()<<"     test  ["<<part->GetFlow(1)<<", "<<part->GetFlow(2)<<"]"
    //	     <<" in list with "<<m_inlist->size()<<" members.\n";
    if (part->GetFlow(2)==ref1 || part->GetFlow(1)==ref2) {
      m_inlist->erase(pit);
      return part;
    }
  }
  return NULL;
}
