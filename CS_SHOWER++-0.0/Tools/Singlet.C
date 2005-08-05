#include "Singlet.H"
#include "Parton.H"
#include "Particle.H"
#include "Message.H"
#include <list>

using namespace CS_SHOWER;
using namespace ATOOLS;
using namespace std;

std::ostream& CS_SHOWER::operator<<(std::ostream& str, Singlet & singlet) {
  str<<"Singlet parton list from CS_Shower : "<<endl;
  Parton * part;
  for (PLiter plit=singlet.begin();plit!=singlet.end();plit++) {
    part = (*plit);
    str<<part<<" "<<(*part);
  }
  str<<"-------------------------------------------------------------------------"<<endl;
  return str;
}

std::ostream& CS_SHOWER::operator<<(std::ostream & str,All_Singlets & all) {
  str<<"Singlet list from CS_Shower : "<<endl;
  Singlet * sing;
  for (ASiter asit=all.begin();asit!=all.end();asit++) {
    sing = (*asit);
    str<<sing<<" "<<sing->size()<<" "<<(*sing);
  }
  str<<"-------------------------------------------------------------------------"<<endl;
  return str;
}


Singlet::~Singlet() {
  if (!empty()) {
    PLiter plit = begin();
    do {
      if ((*plit)) { delete (*plit); (*plit) = NULL; }
      plit = erase(plit);
    } while (plit!=end());
    clear();
  }
}

int Singlet::SplitParton(PLiter & plit, Parton * part1, Parton * part2) 
{
  Parton * mother = (*plit);
  Flavour flav    = mother->GetFlavour(), flav1 = part1->GetFlavour(), flav2 = part2->GetFlavour();

  PLiter pos1,pos2;
  plit = insert(plit,part1);
  pos1 = plit;
  plit++;
  plit = insert(plit,part2);
  pos2 = plit;

  if (mother->GetType()==pst::FS) {
    if (flav.IsGluon()) { 
      if (part1->GetType()==pst::FS && part2->GetType()==pst::FS) { 
	if (flav1.IsQuark() && flav2.IsQuark()) {
	  if (!flav1.IsAnti()) Swap(pos1,pos2);
	  plit++;
	  delete mother; 
	  plit = erase(plit);
	  return 1;
	}
      }
    }
    else if (flav.IsQuark()) {
      if (flav.IsAnti()) {
	if (flav1.IsQuark()) Swap(pos1,pos2);
      }
      else {
	if (flav2.IsQuark()) Swap(pos1,pos2);
      }
    }
  }
  plit++;
  delete mother; 
  plit = erase(plit);
  return 0;
}

Singlet * Singlet::SplitList(PLiter plit) 
{
  Singlet * singlet = new Singlet();
  PLiter plit1;
  for (plit1=begin();plit1!=end();plit1++) {
    if (plit1==plit) break;
  }
  if (plit1==end()) {
    msg.Error()<<"ERROR in Singlet::SplitList : "<<endl
	       <<"   iterator not found : "<<endl
	       <<(*plit)<<endl<<(*this)<<endl;
    abort();
  }
  do {
    singlet->push_back((*plit1));
    plit1 = erase(plit1);
  } while (plit1!=end());
  return singlet;
}

void Singlet::ReshuffleList(PLiter plit) 
{}

void Singlet::ExtractFSPartons(ATOOLS::Blob * blob) 
{
  Particle * part;
  int flow1, flow2;
  PLiter plit1;
  for (PLiter plit=begin();plit!=end();plit++) {
    if ((*plit)->GetType()==pst::FS) {
      part = new Particle(-1,(*plit)->GetFlavour(),(*plit)->Momentum(),'F');
      part->SetNumber(0);
      part->SetStatus(1);
      if (plit==begin()) {
	if ((*plit)->GetFlavour().IsQuark()) {
	  if ((*plit)->GetFlavour().IsAnti()) abort();
	  else {
	    flow1 = m_col;
	    flow2 = 0;
	    part->SetFlow(1,flow1);
	    flow1 = part->GetFlow(1);
	  }
	}
      }
      else {
	plit1 = plit;
	plit1++;
	if (plit1!=end()) {
	  if ((*plit)->GetFlavour().IsQuark()) abort();
	  else {
	    if (flow1!=0) {
	      part->SetFlow(2,flow1);
	      part->SetFlow(1,-1);
	      flow1 = part->GetFlow(1);
	    }
	    else abort();
	  }
	}
	else {
	  if (flow1!=0) {
	    if ((*plit)->GetFlavour().IsQuark() && (*plit)->GetFlavour().IsAnti()) {
	      part->SetFlow(2,flow1);
	      part->SetFlow(1,0);
	    }
	  }
	}
      }
      blob->AddToOutParticles(part);
    }
    else abort();
  }
}

void Singlet::Swap(PLiter pl1,PLiter pl2)
{
  Parton * tmp = (*pl2);
  (*pl2) = (*pl1);
  (*pl1) = tmp;
}
