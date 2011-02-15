#include "CSSHOWER++/Tools/Singlet.H"
#include "CSSHOWER++/Tools/Parton.H"
#include "CSSHOWER++/Showers/Sudakov.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
#include "ATOOLS/Math/ZAlign.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_Limits.H"
#include <list>

using namespace CSSHOWER;
using namespace ATOOLS;
using namespace std;

std::ostream& CSSHOWER::operator<<(std::ostream& str, Singlet & singlet) {
  str<<"Singlet parton list from CS_Shower : "<<&singlet<<endl;
  Parton * part;
  for (PLiter plit=singlet.begin();plit!=singlet.end();plit++) {
    part = (*plit);
    str<<(*part);
  }
  if (singlet.GetSplit() || singlet.GetLeft() || singlet.GetRight() || singlet.GetSpec()) {
    if (singlet.GetSplit()) str<<"Split:  "<<singlet.GetSplit()<<"  ";
    if (singlet.GetLeft()) str<<"Left:  "<<singlet.GetLeft()<<"  ";
    if (singlet.GetRight()) str<<"Right:  "<<singlet.GetRight()<<"  ";
    if (singlet.GetSpec()) str<<"Spec:  "<<singlet.GetSpec()<<"  ";
    str<<"\n";
  }
  str<<"-------------------------------------------------------------------------"<<endl;
  return str;
}

std::ostream& CSSHOWER::operator<<(std::ostream & str,All_Singlets & all) {
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
      if ((*plit)) { 
	delete (*plit); (*plit) = NULL; 
      }
       plit = erase(plit);
    } while (plit!=end());
    clear();
  }
  if (p_ref) delete p_ref;
}

Singlet *Singlet::RefCopy(All_Singlets *const all,std::map<Parton*,Parton*> &pmap)
{
  if (p_ref) delete p_ref;  
  p_ref = new Singlet();
  all->push_back(p_ref);
  p_ref->p_all=all;
  p_ref->p_ms=p_ms;
  p_ref->p_jf=p_jf;
  for (const_iterator it(begin());it!=end();++it) {
    Parton *c(new Parton(**it));
    p_ref->push_back(c);
    c->SetSing(p_ref);
    pmap[*it]=c;
  }
  for (const_iterator it(begin());it!=end();++it) {
    Parton *c(pmap[*it]);
    std::map<Parton*,Parton*>::iterator lit(pmap.find((*it)->GetLeft()));
    if (lit!=pmap.end()) c->SetLeft(lit->second);
    std::map<Parton*,Parton*>::iterator rit(pmap.find((*it)->GetRight()));
    if (rit!=pmap.end()) c->SetRight(rit->second);
  }
  return p_ref;
}

Parton *Singlet::IdParton(const size_t &id) const
{
  for (const_iterator it(begin());it!=end();++it)
    if ((*it)->Id()==id) return *it;
  return NULL;
}

bool Singlet::JetVeto(Sudakov *const sud) const
{
  DEBUG_FUNC("");
  msg_Debugging()<<*(Singlet*)this<<"\n";
  for (const_iterator iit(begin());iit!=end();++iit) {
    bool ii((*iit)->GetType()==pst::IS);
    Flavour fi((*iit)->GetFlavour());
    for (const_iterator jit(iit);jit!=end();++jit) {
      bool ji((*jit)->GetType()==pst::IS);
      Flavour fj((*jit)->GetFlavour());
      if (jit==iit || (ii&&ji)) continue;
      for (const_iterator kit(begin());kit!=end();++kit) {
	if (kit==iit || kit==jit) continue;
	bool ki((*kit)->GetType()==pst::IS);
	cstp::code et((ii||ji)?(ki?cstp::II:cstp::IF):(ki?cstp::FI:cstp::FF));
	if (sud->HasKernel(fi,fj,(*kit)->GetFlavour(),et)) {
	  double q2ijk(p_jf->Qij2(ii?-(*iit)->Momentum():(*iit)->Momentum(),
				  ji?-(*jit)->Momentum():(*jit)->Momentum(),
				  ki?-(*kit)->Momentum():(*kit)->Momentum(),
				  ii?fi.Bar():fi,ji?fj.Bar():fj,p_jf->DR()));
 	  msg_Debugging()<<"Q_{"<<(*iit)->Id()<<(*jit)->Id()
			 <<","<<(*kit)->Id()<<"} = "<<sqrt(q2ijk)<<"\n";
	  if (q2ijk<(*kit)->KtVeto()) return false;
	}
	else {
	  msg_Debugging()<<"No kernel for "<<fi<<" "<<fj<<" <-> "
			 <<(*kit)->GetFlavour()<<" ("<<et<<")\n";
	}
      }
    }
  }
  msg_Debugging()<<"--- Jet veto ---\n";
  return true;
}

void Singlet::SetColours(Singlet *const sing,Parton *const np,
			 const int oldc[2],const int newc[2])
{
//   msg_Debugging()<<METHOD<<"(): np = "<<np<<" {\n";
//   msg_Indent();
//   msg_Debugging()<<"Shift ("<<oldc[0]<<","<<oldc[1]
// 		 <<") -> ("<<newc[0]<<","<<newc[1]<<")\n";
//   msg_Debugging()<<*sing;
  for (PLiter pit(sing->begin());pit!=sing->end();++pit) {
    if (*pit==np) continue;
    for (int i(0);i<2;++i)
      if ((*pit)->GetFlow(i+1)==oldc[i]) 
	(*pit)->SetFlow(i+1,newc[i]);
  }
  for (PLiter lit(sing->begin());lit!=sing->end();++lit)
    if ((*lit)->GetFlavour().StrongCharge()==8 ||
	((*lit)->GetType()==pst::FS?
	 (*lit)->GetFlavour().StrongCharge()>0:
	 (*lit)->GetFlavour().StrongCharge()<0))
      for (PLiter rit(sing->begin());rit!=sing->end();++rit)
	if ((*lit)->GetFlow(1)==(*rit)->GetFlow(2)) {
	  (*lit)->SetLeft(*rit);
	  (*rit)->SetRight(*lit);
	  break;
	}
//   msg_Debugging()<<*sing;
  if (sing->GetLeft()) {
    SetColours(sing->GetLeft()->GetSing(),np?np->GetNext():NULL,oldc,newc);
    if (sing->GetRight()->GetSing()!=sing->GetLeft()->GetSing())
      SetColours(sing->GetRight()->GetSing(),np?np->GetNext():NULL,oldc,newc);
  }
//   msg_Debugging()<<"}\n";
}

void Singlet::SetColours(Parton *const pa,Parton *const p,
			 const int oldc[2],const int newc[2])
{
  msg_Debugging()<<METHOD<<"(): ("<<oldc[0]<<","<<oldc[1]
		 <<") -> ("<<newc[0]<<","<<newc[1]<<")\n";
  if (pa==NULL || pa->GetSing()->GetLeft()==NULL) return;
  SetColours(pa->GetSing()->GetLeft()->GetSing(),p?p->GetNext():NULL,oldc,newc);
  if (pa->GetSing()->GetLeft()->GetSing()!=
      pa->GetSing()->GetRight()->GetSing())
    SetColours(pa->GetSing()->GetRight()->GetSing(),p?p->GetNext():NULL,oldc,newc);
}

int Singlet::SplitParton(Parton * mother, Parton * part1, Parton * part2) 
{
  iterator plit(begin());
  for (;plit!=end();++plit) if (*plit==mother) break;
  if (plit==end()) THROW(fatal_error,"Internal error");

  Flavour flav    = mother->GetFlavour(), flav1 = part1->GetFlavour(), flav2 = part2->GetFlavour();

  PLiter pos1,pos2;
  plit = insert(plit,part1);
  pos1 = plit;
  plit++;
  plit = insert(plit,part2);
  pos2 = plit;

  part1->SetSing(this);
  part2->SetSing(this);

  if (part2->GetNext()) part2->GetNext()->GetSing()->AddParton(part2->GetNext());

  if (mother->GetType()==pst::IS) {
    if ((flav.IsGluon()  || flav.IsGluino())  && 
	(flav1.IsQuark() || flav1.IsSquark()) && flav1.IsAnti() && 
	(flav2.IsQuark() || flav2.IsSquark())) {  
      std::swap<PLiter>(pos1,pos2);
    }
    else if ((flav.IsQuark()  || flav.IsSquark())   && !flav.IsAnti()  && 
	     (flav1.IsQuark() || flav1.IsSquark()) && !flav1.IsAnti() && 
	     (flav2.IsGluon() || flav2.IsGluino())) {
      std::swap<PLiter>(pos1,pos2);
    }
    else if ((flav.IsQuark()  || flav.IsSquark())  && flav.IsAnti()   && 
	     (flav2.IsQuark() || flav2.IsSquark()) && !flav2.IsAnti() && 
	     (flav1.IsGluon() || flav1.IsGluino())) {
      std::swap<PLiter>(pos1,pos2);
    }
  }
  if (mother->GetType()==pst::FS) {
    if ((flav.IsQuark()  || flav.IsSquark()) && 
	(flav2.IsQuark() || flav2.IsSquark())) {
      if (!flav2.IsAnti()) std::swap<PLiter>(pos1,pos2);
    }
    else if ((flav.IsQuark()  || flav.IsSquark()) && 
	     (flav1.IsQuark() || flav1.IsSquark())) {
      if (flav1.IsAnti()) std::swap<PLiter>(pos1,pos2);
    }
    else if ((flav.IsGluon()  || flav.IsGluino()) && 
	     (flav1.IsQuark() || flav1.IsSquark())) {
      if (!flav1.IsAnti()) std::swap<PLiter>(pos1,pos2);
    }
  }
  
  plit++;
  delete mother; 
  plit = erase(plit);
  if ((flav.IsGluon()  || flav.IsGluino()) && 
      (flav1.IsQuark() || flav1.IsSquark()) && 
      (flav2.IsQuark() || flav2.IsSquark())) { return 1; }
  return 0;
}

void Singlet::ExtractPartons
(ATOOLS::Blob * blob,ATOOLS::Mass_Selector *const ms) 
{
  Particle * part;
  for (PLiter plit=begin();plit!=end();plit++) {
    if ((*plit)->Stat()==1) continue;
      part = new Particle(-1,(*plit)->GetFlavour(),(*plit)->Momentum(),'F');
      part->SetNumber(0);
      if ((*plit)->GetType()==pst::IS) {
	part->SetInfo('I');
	blob->AddToInParticles(part);
      } 
      else blob->AddToOutParticles(part);
      if ((*plit)->GetType()==pst::FS) {
	part->SetFlow(1,(*plit)->GetFlow(1));
	part->SetFlow(2,(*plit)->GetFlow(2));
      }
      else if ((*plit)->GetType()==pst::IS) {
	part->SetFlow(1,(*plit)->GetFlow(2));
	part->SetFlow(2,(*plit)->GetFlow(1));
      }
      part->SetFinalMass(ms->Mass((*plit)->GetFlavour()));
  }
}

void Singlet::RemoveParton(Parton *const p,const int mode)
{
  for (iterator pit(begin());pit!=end();++pit)
    if (*pit==p) {
      if (p->GetNext()) p->GetNext()->GetSing()->
	RemoveParton(p->GetNext(),mode);
      if (mode) {
	if (p->GetPrev()) p->GetPrev()->SetNext(NULL);
	delete p;
      }
      erase(pit);
      return;
    }
  THROW(fatal_error,"Parton not found");
}

void Singlet::AddParton(Parton *const p)
{
  push_back(p);
  p->SetSing(this);
  if (p_left) {
    Parton *np(p->GetNext());
    if (np==NULL) {
      np = new Parton(p->GetFlavour(),p->Momentum(),p->GetType());
      p->SetStat(1);
      p->SetNext(np);
      np->SetPrev(p);
      np->SetStart(p->KtStart());
      np->SetVeto(p->KtVeto());
      np->SetKtMax(p->KtMax());
    }
    p_left->GetSing()->AddParton(np);
  }
}

bool Singlet::RearrangeColours(Parton * mother, Parton * daughter1, Parton * daughter2)
{
  int oldc[2]={mother->GetFlow(1),mother->GetFlow(2)}, newc[2]={oldc[0],oldc[1]};
  if (mother->GetType()==pst::IS) std::swap<Parton*>(mother,daughter1);
  Flavour mo(mother->GetFlavour()), d1(daughter1->GetFlavour()), d2(daughter2->GetFlavour());
  if (mother->GetType()==pst::FS) {
    if (mo.StrongCharge()==-3) {
      if (d1.StrongCharge()==-3) {
	if (d2.StrongCharge()==8) {
	  newc[1]=daughter1->GetFlow(2);
	  SetColours(daughter1,NULL,newc,oldc);
	  daughter2->SetRightOf(mother);
	  return true;
	}
	else if (d2.StrongCharge()==0) {
	  daughter1->SetRightOf(mother);
	  return true;
	}
      }
      else if (d2.StrongCharge()==-3) {
	if (d1.StrongCharge()==8) {
	  newc[1]=daughter2->GetFlow(2);
	  SetColours(daughter1,NULL,newc,oldc);
	  daughter1->SetRightOf(mother);
	  return true;
	}
	else if (d1.StrongCharge()==0) {
	  daughter2->SetRightOf(mother);
	  return true;
	}
      }
    }
    else if (mo.StrongCharge()==3) {
      if (d1.StrongCharge()==3) {
	if (d2.StrongCharge()==8) {
	  newc[0]=daughter1->GetFlow(1);
	  SetColours(daughter1,NULL,newc,oldc);
	  daughter2->SetLeftOf(mother);
	  return true;
	}
	else if (d2.StrongCharge()==0) {
	  daughter1->SetLeftOf(mother);
	  return true;
	}
      }
      else if (d2.StrongCharge()==3) {
	if (d1.StrongCharge()==8) {
	  newc[0]=daughter2->GetFlow(1);
	  SetColours(daughter1,NULL,newc,oldc);
	  daughter1->SetLeftOf(mother);
	  return true;
	}
	else if (d1.StrongCharge()==0) {
	  daughter2->SetLeftOf(mother);
	  return true;
	}
      }
    }
    else if (mo.StrongCharge()==8) {
      if (d1.StrongCharge()==3 && 
	  d2.StrongCharge()==-3) {  
	daughter1->SetLeftOf(mother);
	daughter2->SetRightOf(mother);
	return true;
      }
      else if (d1.StrongCharge()==-3 && 
	       d2.StrongCharge()==3) {  
	daughter2->SetLeftOf(mother);
	daughter1->SetRightOf(mother);
	return true;
      }
      else if (d1.StrongCharge()==8 && 
	       d2.StrongCharge()==8) {
 	newc[0]=daughter1->GetFlow(1);
	SetColours(daughter1,NULL,newc,oldc);
	daughter2->SetLeftOf(mother);
	daughter1->SetRightOf(mother);
	return true;
      }
    }
    else if (mo.StrongCharge()==0) {
      if (abs(d1.StrongCharge())==3 && 
	  abs(d2.StrongCharge())==3) {  
	return true;
      }
      else if (d1.StrongCharge()==0 && 
	       d2.StrongCharge()==0) {
	return true;
      }
    }
  }
  else if (daughter1->GetType()==pst::IS) {
    if (d1.StrongCharge()==-3) {
      if (mo.StrongCharge()==-3) {
	if (d2.StrongCharge()==8) {
	  newc[0]=mother->GetFlow(1);
	  SetColours(mother,NULL,newc,oldc);
	  daughter2->SetLeftOf(daughter1);
	  return true;
	}
	else if (d2.StrongCharge()==0) {
	  mother->SetLeftOf(daughter1);
	  return true;
	}
      }
      else if (d2.StrongCharge()==3) {
	if (mo.StrongCharge()==8) {
	  newc[1]=mother->GetFlow(2);
	  SetColours(mother,NULL,newc,oldc);
	  mother->SetLeftOf(daughter1);
	  return true;
	}
	if (mo.StrongCharge()==0) {
	  daughter2->SetLeftOf(daughter1);
	  return true;
	}
      }
    }
    else if (d1.StrongCharge()==3) {
      if (mo.StrongCharge()==3) {
	if (d2.StrongCharge()==8) {
	  newc[1]=mother->GetFlow(2);
	  SetColours(mother,NULL,newc,oldc);
	  daughter2->SetRightOf(daughter1);
	  return true;
	}
	else if (d2.StrongCharge()==0) {
	  mother->SetRightOf(daughter1);
	  return true;
	}
      }
      else if (d2.StrongCharge()==-3) {
	if (mo.StrongCharge()==8) {
	  newc[0]=mother->GetFlow(1);
	  SetColours(mother,NULL,newc,oldc);
	  mother->SetRightOf(daughter1);
	  return true;
	}
	else if (mo.StrongCharge()==0) {
	  daughter2->SetRightOf(daughter1);
	  return true;
	}
      }
    }
    else if (d1.StrongCharge()==8) {
      if (abs(mo.StrongCharge())==3) {
	if (d2.StrongCharge()==-3) {
	  mother->SetLeftOf(daughter1);
	  daughter2->SetRightOf(daughter1);
	  return true;
	}
	else if (d2.StrongCharge()==3) {
	  mother->SetRightOf(daughter1);
	  daughter2->SetLeftOf(daughter1);
	  return true;
	}
      }
      else if (mo.StrongCharge()==8 && 
	       d2.StrongCharge()==8) {
	newc[0]=mother->GetFlow(1);
	SetColours(mother,NULL,newc,oldc);
	mother->SetRightOf(daughter1);
	daughter2->SetLeftOf(daughter1);
	return true;
      }
    }
    else if (d1.StrongCharge()==0) {
      if (abs(mo.StrongCharge())==3) {
	if (d2.StrongCharge()==-3) {
	  return true;
	}
	else if (d2.StrongCharge()==3) {
	  return true;
	}
      }
    }
  }
  return false;
}

bool Singlet::ArrangeColours(Parton * mother, Parton * daughter1, Parton * daughter2)
{
  daughter1->SetSing(this);
  daughter2->SetSing(this);
  int oldc[2]={mother->GetFlow(1),mother->GetFlow(2)}, newc[2]={oldc[0],oldc[1]};
  if (mother->GetType()==pst::IS) std::swap<Parton*>(mother,daughter1);
  Flavour mo(mother->GetFlavour()), d1(daughter1->GetFlavour()), d2(daughter2->GetFlavour());
  if (mother->GetType()==pst::FS) {
    if (mo.StrongCharge()==-3) {
      if (d1.StrongCharge()==-3) {
	if (d2.StrongCharge()==8) {
	  daughter2->SetFlow(2,mother->GetFlow(2));
	  daughter2->SetMEFlow(2,mother->GetMEFlow(2));
	  daughter2->SetFlow(1,-1);
	  daughter1->SetFlow(2,daughter2->GetFlow(1));
	  newc[1]=daughter1->GetFlow(2);
	  daughter2->UpdateColours();
	  SetColours(daughter1,daughter2,oldc,newc);
	  mother->SetRightOf(daughter2);
	  daughter2->SetLeft(daughter1);
	  daughter1->SetRight(daughter2);
	  return true;
	}
	else if (d2.StrongCharge()==0) {
	  daughter1->SetFlow(2,mother->GetFlow(2));
	  daughter1->SetMEFlow(2,mother->GetMEFlow(2));
	  daughter2->SetFlow(1,0);
	  daughter2->SetFlow(2,0);
	  mother->SetRightOf(daughter1);
	  return true;
	}
      }
      else if (d2.StrongCharge()==-3) {
	if (d1.StrongCharge()==8) {
	  daughter1->SetFlow(2,mother->GetFlow(2));
	  daughter1->SetMEFlow(2,mother->GetMEFlow(2));
	  daughter1->SetFlow(1,-1);
	  daughter2->SetFlow(2,daughter1->GetFlow(1));
	  newc[1]=daughter2->GetFlow(2);
	  daughter2->UpdateColours();
	  SetColours(daughter1,daughter2,oldc,newc);
	  mother->SetRightOf(daughter1);
	  daughter1->SetLeft(daughter2);
	  daughter2->SetRight(daughter1);
	  return true;
	}
	else if (d1.StrongCharge()==0) {
	  daughter2->SetFlow(2,mother->GetFlow(2));
	  daughter2->SetMEFlow(2,mother->GetMEFlow(2));
	  daughter1->SetFlow(1,0);
	  daughter1->SetFlow(2,0);
	  mother->SetRightOf(daughter2);
	  return true;
	}
      }
    }
    else if (mo.StrongCharge()==3) {
      if (d1.StrongCharge()==3) {
	if (d2.StrongCharge()==8) {
	  daughter2->SetFlow(1,mother->GetFlow(1));
	  daughter2->SetMEFlow(1,mother->GetMEFlow(1));
	  daughter2->SetFlow(2,-1);
	  daughter1->SetFlow(1,daughter2->GetFlow(2));
	  newc[0]=daughter1->GetFlow(1);
	  daughter2->UpdateColours();
	  SetColours(daughter1,daughter2,oldc,newc);
	  mother->SetLeftOf(daughter2);
	  daughter2->SetRight(daughter1);
	  daughter1->SetLeft(daughter2);
	  return true;
	}
	else if (d2.StrongCharge()==0) {
	  daughter1->SetFlow(1,mother->GetFlow(1));
	  daughter1->SetMEFlow(1,mother->GetMEFlow(1));
	  daughter2->SetFlow(1,0);
	  daughter2->SetFlow(2,0);
	  mother->SetLeftOf(daughter1);
	  return true;
	}
      }
      else if (d2.StrongCharge()==3) {
	if (d1.StrongCharge()==8) {
	  daughter1->SetFlow(1,mother->GetFlow(1));
	  daughter1->SetMEFlow(1,mother->GetMEFlow(1));
	  daughter1->SetFlow(2,-1);
	  daughter2->SetFlow(1,daughter1->GetFlow(2));
	  newc[0]=daughter2->GetFlow(1);
	  daughter2->UpdateColours();
	  SetColours(daughter1,daughter2,oldc,newc);
	  mother->SetLeftOf(daughter1);
	  daughter1->SetRight(daughter2);
	  daughter2->SetLeft(daughter1);
	  return true;
	}
	else if (d1.StrongCharge()==0) {
	  daughter2->SetFlow(1,mother->GetFlow(1));
	  daughter2->SetMEFlow(1,mother->GetMEFlow(1));
	  daughter1->SetFlow(1,0);
	  daughter1->SetFlow(2,0);
	  mother->SetLeftOf(daughter2);
	  return true;
	}
      }
    }
    else if (mo.StrongCharge()==8) {
      if (d1.StrongCharge()==3 && 
	  d2.StrongCharge()==-3) {  
	daughter1->SetFlow(1,mother->GetFlow(1));
	daughter1->SetMEFlow(1,mother->GetMEFlow(1));
	daughter1->SetFlow(2,0);
	daughter2->SetFlow(1,0);
	daughter2->SetFlow(2,mother->GetFlow(2));
	daughter2->SetMEFlow(2,mother->GetMEFlow(2));
	mother->SetLeftOf(daughter1);
	mother->SetRightOf(daughter2);
	return true;
      }
      else if (d1.StrongCharge()==-3 && 
	       d2.StrongCharge()==3) {  
	daughter2->SetFlow(1,mother->GetFlow(1));
	daughter2->SetMEFlow(1,mother->GetMEFlow(1));
	daughter2->SetFlow(2,0);
	daughter1->SetFlow(1,0);
	daughter1->SetFlow(2,mother->GetFlow(2));
	daughter1->SetMEFlow(2,mother->GetMEFlow(2));
	mother->SetLeftOf(daughter2);
	mother->SetRightOf(daughter1);
	return true;
      }
      else if (d1.StrongCharge()==8 && 
	       d2.StrongCharge()==8) {
	daughter1->SetFlow(2,mother->GetFlow(2));
	daughter1->SetMEFlow(2,mother->GetMEFlow(2));
	daughter1->SetFlow(1,-1);
	daughter2->SetFlow(2,daughter1->GetFlow(1));
	daughter2->SetFlow(1,mother->GetFlow(1));
	daughter2->SetMEFlow(1,mother->GetMEFlow(1));
 	newc[0]=daughter1->GetFlow(1);
	daughter2->UpdateColours();
	SetColours(daughter1,daughter2,oldc,newc);
	mother->SetLeftOf(daughter2);
	daughter2->SetRight(daughter1);
	daughter1->SetLeft(daughter2);
	mother->SetRightOf(daughter1);
	return true;
      }
    }
    else if (mo.StrongCharge()==0) {
      if (abs(d1.StrongCharge())==3 && 
	  abs(d2.StrongCharge())==3) {  
	daughter1->SetFlow(1,-1);
	daughter1->SetFlow(2,0);
	daughter2->SetFlow(2,daughter1->GetFlow(1));
	daughter2->SetFlow(1,0);
	daughter1->SetLeft(daughter2);
	daughter2->SetRight(daughter1);
	return true;
      }
      else if (d1.StrongCharge()==0 && 
	       d2.StrongCharge()==0) {
	daughter1->SetFlow(1,0);
	daughter1->SetFlow(2,0);
	daughter2->SetFlow(1,0);
	daughter2->SetFlow(2,0);
	return true;
      }
    }
  }
  else if (daughter1->GetType()==pst::IS) {
    if (d1.StrongCharge()==-3) {
      if (mo.StrongCharge()==-3) {
	if (d2.StrongCharge()==8) {
	  mother->SetFlow(1,-1);
	  daughter2->SetFlow(1,daughter1->GetFlow(1));
	  daughter2->SetMEFlow(1,daughter1->GetMEFlow(1));
	  daughter2->SetFlow(2,mother->GetFlow(1));
	  newc[0]=mother->GetFlow(1);
	  daughter2->UpdateColours();
	  SetColours(mother,daughter2,oldc,newc);
	  daughter1->SetLeftOf(daughter2);
	  daughter2->SetRight(mother);
	  mother->SetLeft(daughter2);
	  return true;
	}
	else if (d2.StrongCharge()==0) {
	  mother->SetFlow(1,daughter1->GetFlow(1));
	  mother->SetMEFlow(1,daughter1->GetMEFlow(1));
	  daughter2->SetFlow(1,0);
	  daughter2->SetFlow(2,0);
	  daughter1->SetLeftOf(mother);
	  return true;
	}
      }
      else if (d2.StrongCharge()==3) {
	if (mo.StrongCharge()==8) {
	  mother->SetFlow(1,daughter1->GetFlow(1));
	  mother->SetMEFlow(1,daughter1->GetMEFlow(1));
	  mother->SetFlow(2,-1);
	  daughter2->SetFlow(1,mother->GetFlow(2));
	  newc[1]=mother->GetFlow(2);
	  daughter2->UpdateColours();
	  SetColours(mother,daughter2,oldc,newc);
	  daughter1->SetLeftOf(mother);
	  mother->SetRight(daughter2);
	  daughter2->SetLeft(mother);
	  return true;
	}
	if (mo.StrongCharge()==0) {
	  daughter2->SetFlow(1,daughter1->GetFlow(1));
	  daughter2->SetMEFlow(1,daughter1->GetMEFlow(1));
	  mother->SetFlow(1,0);
	  mother->SetFlow(2,0);
	  daughter1->SetLeftOf(daughter2);
	  return true;
	}
      }
    }
    else if (d1.StrongCharge()==3) {
      if (mo.StrongCharge()==3) {
	if (d2.StrongCharge()==8) {
	  mother->SetFlow(2,-1);
	  daughter2->SetFlow(1,mother->GetFlow(2));
	  daughter2->SetFlow(2,daughter1->GetFlow(2));
	  daughter2->SetMEFlow(2,daughter1->GetMEFlow(2));
	  newc[1]=mother->GetFlow(2);
	  daughter2->UpdateColours();
	  SetColours(mother,daughter2,oldc,newc);
	  daughter1->SetRightOf(daughter2);
	  daughter2->SetLeft(mother);
	  mother->SetRight(daughter2);
	  return true;
	}
	else if (d2.StrongCharge()==0) {
	  mother->SetFlow(2,daughter1->GetFlow(2));
	  mother->SetMEFlow(2,daughter1->GetMEFlow(2));
	  daughter2->SetFlow(1,0);
	  daughter2->SetFlow(2,0);
	  daughter1->SetRightOf(mother);
	  return true;
	}
      }
      else if (d2.StrongCharge()==-3) {
	if (mo.StrongCharge()==8) {
	  mother->SetFlow(1,-1);
	  mother->SetFlow(2,daughter1->GetFlow(2));
	  mother->SetMEFlow(2,daughter1->GetMEFlow(2));
	  daughter2->SetFlow(2,mother->GetFlow(1));
	  newc[0]=mother->GetFlow(1);
	  daughter2->UpdateColours();
	  SetColours(mother,daughter2,oldc,newc);
	  daughter1->SetRightOf(mother);
	  mother->SetLeft(daughter2);
	  daughter2->SetRight(mother);
	  return true;
	}
	else if (mo.StrongCharge()==0) {
	  daughter2->SetFlow(2,daughter1->GetFlow(2));
	  daughter2->SetMEFlow(2,daughter1->GetMEFlow(2));
	  mother->SetFlow(1,0);
	  mother->SetFlow(2,0);
	  daughter1->SetRightOf(daughter2);
	  return true;
	}
      }
    }
    else if (d1.StrongCharge()==8) {
      if (abs(mo.StrongCharge())==3) {
	if (d2.StrongCharge()==-3) {
	  daughter2->SetFlow(1,0);
	  daughter2->SetFlow(2,daughter1->GetFlow(2));
	  daughter2->SetMEFlow(2,daughter1->GetMEFlow(2));
	  mother->SetFlow(2,0);
	  mother->SetFlow(1,daughter1->GetFlow(1));
	  mother->SetMEFlow(1,daughter1->GetMEFlow(1));
	  daughter1->SetLeftOf(mother);
	  daughter1->SetRightOf(daughter2);
	  return true;
	}
	else if (d2.StrongCharge()==3) {
	  daughter2->SetFlow(2,0);
	  daughter2->SetFlow(1,daughter1->GetFlow(1));
	  daughter2->SetMEFlow(1,daughter1->GetMEFlow(1));
	  mother->SetFlow(1,0);
	  mother->SetFlow(2,daughter1->GetFlow(2));
	  mother->SetMEFlow(2,daughter1->GetMEFlow(2));
	  daughter1->SetRightOf(mother);
	  daughter1->SetLeftOf(daughter2);
	  return true;
	}
      }
      else if (mo.StrongCharge()==8 && 
	       d2.StrongCharge()==8) {
	mother->SetFlow(2,daughter1->GetFlow(2));
	mother->SetMEFlow(2,daughter1->GetMEFlow(2));
	mother->SetFlow(1,-1);	
	daughter2->SetFlow(1,daughter1->GetFlow(1));
	daughter2->SetMEFlow(1,daughter1->GetMEFlow(1));
	daughter2->SetFlow(2,mother->GetFlow(1));
	newc[0]=mother->GetFlow(1);
	daughter2->UpdateColours();
	SetColours(mother,daughter2,oldc,newc);
	daughter1->SetRightOf(mother);
	mother->SetLeft(daughter2);
	daughter2->SetRight(mother);
	daughter1->SetLeftOf(daughter2);
	return true;
      }
    }
    else if (d1.StrongCharge()==0) {
      if (abs(mo.StrongCharge())==3) {
	if (d2.StrongCharge()==-3) {
	  daughter2->SetFlow(1,0);
	  daughter2->SetFlow(2,-1);
	  mother->SetFlow(1,daughter2->GetFlow(2));
	  mother->SetFlow(2,0);
	  daughter2->SetRight(mother);
	  mother->SetLeft(daughter2);
	  return true;
	}
	else if (d2.StrongCharge()==3) {
	  daughter2->SetFlow(2,0);
	  daughter2->SetFlow(1,-1);
	  mother->SetFlow(2,daughter2->GetFlow(1));
	  mother->SetFlow(1,0);
	  daughter2->SetLeft(mother);
	  mother->SetRight(daughter2);
	  return true;
	}
      }
    }
  }
  return false;
} 

void Singlet::BoostAllFS(Parton *l,Parton *r,Parton *s,Parton *f,
			 const Flavour &mo,const int mode)
{
  if (p_all==NULL) return;
  if (mode&2) {
    if (mode&1) {
      Vec4D pa(l->Momentum()), pb(s->Momentum());
      double ma2(p_ms->Mass2(l->GetFlavour()));
      double mb2(p_ms->Mass2(s->GetFlavour()));
      ZAlign lt(pa,pb,ma2,mb2);
      l->SetMomentum(lt.PaNew());
      for (All_Singlets::const_iterator asit(p_all->begin());
	   asit!=p_all->end();++asit) {
	for (PLiter plit((*asit)->begin());plit!=(*asit)->end();++plit) {
	  if (*plit==l) continue;
	  (*plit)->SetMomentum(lt.Align((*plit)->Momentum()));
	}
      }
    }
    else {
      Parton *b(NULL);
      for (PLiter plit(s->GetSing()->begin());plit!=s->GetSing()->end();++plit)
	if ((*plit)->GetType()==pst::IS && *plit!=s) {
	  b=*plit;
	  break;
	}
      if (b==NULL) THROW(fatal_error,"Corrupted singlet");
      Vec4D pb(b->Momentum()), ps(s->Momentum());
      double ma2(p_ms->Mass2(s->GetFlavour()));
      double mb2(p_ms->Mass2(b->GetFlavour()));
      if (ma2==0.0) return;
      ZAlign lt(ps,pb,ma2,mb2);
      for (All_Singlets::const_iterator asit(p_all->begin());
	   asit!=p_all->end();++asit) {
	for (PLiter plit((*asit)->begin());plit!=(*asit)->end();++plit) {
	  (*plit)->SetMomentum(lt.Align((*plit)->Momentum()));
	}
      }
    }
  }
  else {
    if (mode&1) {
      if (f->Kin()==1) {
      Parton *b(NULL);
      for (PLiter plit(f->GetSing()->begin());plit!=f->GetSing()->end();++plit)
	if ((*plit)->GetType()==pst::IS && *plit!=f) {
	  b=*plit;
	  break;
	}
      if (b==NULL) THROW(fatal_error,"Corrupted singlet");
      Vec4D pb(-b->Momentum()), pl(-l->Momentum());
      double ma2(p_ms->Mass2(l->GetFlavour())), maj2(p_ms->Mass2(mo));
      double mb2(p_ms->Mass2(b->GetFlavour()));
      if (ma2==0.0 && maj2==0.0) return;
      ZAlign lt(-pl,-pb,ma2,mb2);
      l->SetMomentum(-lt.PaNew());
      for (All_Singlets::const_iterator asit(p_all->begin());
	   asit!=p_all->end();++asit) {
	for (PLiter plit((*asit)->begin());plit!=(*asit)->end();++plit) {
	  if (*plit==l) continue;
	  (*plit)->SetMomentum(lt.Align((*plit)->Momentum()));
	}
      }
      }
      else {
      Parton *b(NULL);
      for (PLiter plit(f->GetSing()->begin());plit!=f->GetSing()->end();++plit)
	if ((*plit)->GetType()==pst::IS && *plit!=f) {
	  b=*plit;
	  break;
	}
      if (b==NULL) THROW(fatal_error,"Corrupted singlet");
      Vec4D pb(b->Momentum()), pl(l->Momentum());
      double ma2(p_ms->Mass2(l->GetFlavour()));
      double mb2(p_ms->Mass2(b->GetFlavour()));
      ZAlign lt(pl,pb,ma2,mb2);
      l->SetMomentum(lt.PaNew());
      for (All_Singlets::const_iterator asit(p_all->begin());
	   asit!=p_all->end();++asit) {
	for (PLiter plit((*asit)->begin());plit!=(*asit)->end();++plit) {
 	  if (*plit==l) continue;
	  (*plit)->SetMomentum(lt.Align((*plit)->Momentum()));
	}
      }
      }
    }
  }
}

void Singlet::BoostBackAllFS(Parton *l,Parton *r,Parton *s,Parton *f,
			     const Flavour &mo,const int mode)
{
  if (p_all==NULL) return;
  if (mode&2) {
    if (mode&1) {
      Vec4D pa(-l->Momentum()), pb(-s->Momentum());
      Vec4D pi(r->Momentum()), pai(pa+pi);
      double sai(pai.Abs2()), mb2(p_ms->Mass2(s->GetFlavour()));
      ZAlign lt(-pai,-pb,sai,mb2);
      for (All_Singlets::const_iterator asit(p_all->begin());
	   asit!=p_all->end();++asit) {
	for (PLiter plit((*asit)->begin());plit!=(*asit)->end();++plit) {
	  (*plit)->SetMomentum(lt.Align((*plit)->Momentum()));
	}
      }
    }
    else {
      Parton *b(NULL);
      for (PLiter plit(s->GetSing()->begin());plit!=s->GetSing()->end();++plit)
	if ((*plit)->GetType()==pst::IS && *plit!=s) {
	  b=*plit;
	  break;
	}
      if (b==NULL) THROW(fatal_error,"Corrupted singlet");
      Vec4D pb(-b->Momentum()), pa(-s->Momentum());
      Vec4D pi(l->Momentum()), pj(r->Momentum()), Q(pi+pj+pa);
      double ma2(p_ms->Mass2(s->GetFlavour())), mij2(p_ms->Mass2(mo));
      double mb2(p_ms->Mass2(b->GetFlavour()));
      if (ma2==0.0) return;
      double sij((pi+pj).Abs2()), Q2(Q.Abs2());
      double lrat((sqr(Q2-mij2-ma2)-4.0*mij2*ma2)/
		  (sqr(Q2-sij-ma2)-4.0*sij*ma2));
      Vec4D pat(sqrt(lrat)*(pa-(Q*pa/Q2)*Q)+(Q2+ma2-mij2)/(2.*Q2)*Q);
      ZAlign lt(-pat,-pb,ma2,mb2);
      for (All_Singlets::const_iterator asit(p_all->begin());
	   asit!=p_all->end();++asit) {
	for (PLiter plit((*asit)->begin());plit!=(*asit)->end();++plit) {
	  (*plit)->SetMomentum(lt.Align((*plit)->Momentum()));
	}
      }
    }
  }
  else {
    if (mode&1) {
      if (f->Kin()==1) {
      Parton *b(NULL);
      for (PLiter plit(f->GetSing()->begin());plit!=f->GetSing()->end();++plit)
	if ((*plit)->GetType()==pst::IS && *plit!=f) {
	  b=*plit;
	  break;
	}
      if (b==NULL) THROW(fatal_error,"Corrupted singlet");
      Vec4D pb(-b->Momentum()), pa(-l->Momentum());
      Vec4D pj(r->Momentum()), pk(s->Momentum()), Q(pa+pj+pk);
      double ma2(p_ms->Mass2(l->GetFlavour())), maj2(p_ms->Mass2(mo));
      double mb2(p_ms->Mass2(b->GetFlavour()));
      double mk2(p_ms->Mass2(s->GetFlavour()));
      if (ma2==0.0 && maj2==0.0) return;
      double sjk((pj+pk).Abs2()), Q2(Q.Abs2());
      double lrat((sqr(Q2-mk2-ma2)-4.0*mk2*ma2)/
		  (sqr(Q2-sjk-ma2)-4.0*sjk*ma2));
      Vec4D pat(sqrt(lrat)*(pa-(Q*pa/Q2)*Q)+(Q2+ma2-mk2)/(2.*Q2)*Q);
      ZAlign lt(-pat,-pb,ma2,mb2);
      for (All_Singlets::const_iterator asit(p_all->begin());
	   asit!=p_all->end();++asit) {
	for (PLiter plit((*asit)->begin());plit!=(*asit)->end();++plit) {
	  (*plit)->SetMomentum(lt.Align((*plit)->Momentum()));
	}
      }
      }
      else {
      Parton *b(NULL);
      for (PLiter plit(f->GetSing()->begin());plit!=f->GetSing()->end();++plit)
	if ((*plit)->GetType()==pst::IS && *plit!=f) {
	  b=*plit;
	  break;
	}
      if (b==NULL) THROW(fatal_error,"Corrupted singlet");
      Vec4D pb(-b->Momentum()), pa(-l->Momentum());
      Vec4D pj(r->Momentum()), pk(s->Momentum()), Q(pa+pj+pk);
      double mb2(p_ms->Mass2(b->GetFlavour())), maj2(p_ms->Mass2(mo));
      double mk2(p_ms->Mass2(s->GetFlavour()));
      double saj((pa+pj).Abs2()), Q2(Q.Abs2());
      double lrat((sqr(Q2-maj2-mk2)-4.0*maj2*mk2)/
		  (sqr(Q2-saj-mk2)-4.0*saj*mk2));
      Vec4D pikt(sqrt(lrat)*(pk-(Q*pk/Q2)*Q)+(Q2+mk2-maj2)/(2.*Q2)*Q);
      Vec4D pat(Q-pikt); 
      ZAlign lt(-pat,-pb,maj2,mb2);
      for (All_Singlets::const_iterator asit(p_all->begin());
	   asit!=p_all->end();++asit) {
	for (PLiter plit((*asit)->begin());plit!=(*asit)->end();++plit) {
	  (*plit)->SetMomentum(lt.Align((*plit)->Momentum()));
	}
      }
      }
    }
  }
}

void Singlet::UpdateDaughters()
{
  for (PLiter plit(begin());plit!=end();++plit)
    (*plit)->UpdateDaughters();
}
