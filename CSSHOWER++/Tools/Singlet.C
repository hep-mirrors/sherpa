#include "CSSHOWER++/Tools/Singlet.H"
#include "CSSHOWER++/Tools/Parton.H"
#include "CSSHOWER++/Showers/Sudakov.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
#include "ATOOLS/Math/ZAlign.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
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
  
}

bool Singlet::JetVeto
(Sudakov *const sud,PHASIC::Jet_Finder *const jf,const double &q2,
 const Flavour &fj,const Vec4D &pj) const
{
  for (const_iterator eit(begin());eit!=end();++eit) {
    for (const_iterator sit(begin());sit!=end();++sit) {
      if (eit==sit) continue;
      Parton *e(*eit), *s(*sit);
      bool ei(e->GetType()==pst::IS), si(s->GetType()==pst::IS);
      cstp::code et(ei?(si?cstp::II:cstp::IF):(si?cstp::FI:cstp::FF));
      const SF_E_Map *sfs(sud->HasKernel(e->GetFlavour(),fj,et));
      if (sfs==NULL) continue;
      for (SF_E_Map::const_iterator
	     kit(sfs->begin());kit!=sfs->end();++kit) {
	if (kit->second->Coupling()->AllowSpec(s->GetFlavour())) {
	  Flavour fi(ei?e->GetFlavour().Bar():e->GetFlavour());
	  Flavour fk(si?s->GetFlavour().Bar():s->GetFlavour());
	  Vec4D pi(ei?-e->Momentum():e->Momentum());
	  Vec4D pk(si?-s->Momentum():s->Momentum());
	  if (jf->Qij2(pi,pj,pk,fi,fj)<q2) return false;
	  break;
	}
      }
    }
  }
  return true;
}

void Singlet::SetColours(Singlet *const sing,
			 const int oldc[2],const int newc[2])
{
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Indent();
  msg_Debugging()<<"Shift ("<<oldc[0]<<","<<oldc[1]
		 <<") -> ("<<newc[0]<<","<<newc[1]<<")\n";
  msg_Debugging()<<*sing;
  for (PLiter pit(sing->begin());pit!=sing->end();++pit)
    for (int i(0);i<2;++i)
      if ((*pit)->GetFlow(i+1)==oldc[i]) 
	(*pit)->SetFlow(i+1,newc[i]);
  msg_Debugging()<<*sing;
  if (sing->GetLeft()) {
    SetColours(sing->GetLeft()->GetSing(),oldc,newc);
    if (sing->GetRight()->GetSing()!=sing->GetLeft()->GetSing())
      SetColours(sing->GetRight()->GetSing(),oldc,newc);
  }
  msg_Debugging()<<"}\n";
}

void Singlet::SetColours(Parton *const pa,
			 const int oldc[2],const int newc[2])
{
  msg_Debugging()<<METHOD<<"(): ("<<oldc[0]<<","<<oldc[1]
		 <<") -> ("<<newc[0]<<","<<newc[1]<<")\n";
  if (pa==NULL || pa->GetSing()->GetLeft()==NULL) return;
  SetColours(pa->GetSing()->GetLeft()->GetSing(),oldc,newc);
  if (pa->GetSing()->GetLeft()->GetSing()!=
      pa->GetSing()->GetRight()->GetSing())
    SetColours(pa->GetSing()->GetRight()->GetSing(),oldc,newc);
}

int Singlet::SplitParton(PLiter & plit, Parton * part1, Parton * part2) 
{
  Parton * mother((*plit));

  Flavour flav    = mother->GetFlavour(), flav1 = part1->GetFlavour(), flav2 = part2->GetFlavour();

  PLiter pos1,pos2;
  plit = insert(plit,part1);
  pos1 = plit;
  plit++;
  plit = insert(plit,part2);
  pos2 = plit;

  part1->SetSing(this);
  part2->SetSing(this);

  if (!ArrangeColours(mother,part1,part2)) {
    msg_Error()<<"ERROR in Singlet::SplitParton."<<std::endl
	       <<"   Do not know how to handle this colour flow: "<<std::endl
	       <<"   "<<(*mother)<<"   --> "<<(*part1)<<"       "<<(*part2)<<std::endl
	       <<"   Will abort."<<std::endl;
    abort();
  }
  /*
  std::cout<<"ARRANGED FOR : "<<mother->GetType()<<" : "
	   <<mother->GetFlavour()<<"("<<mother->GetFlow(1)<<","<<mother->GetFlow(2)<<") --> "
  	   <<part1->GetFlavour()<<"("<<part1->GetFlow(1)<<","<<part1->GetFlow(2)<<") + "
  	   <<part2->GetFlavour()<<"("<<part2->GetFlow(1)<<","<<part2->GetFlow(2)<<")"<<std::endl;
  */
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

bool Singlet::ArrangeColours(Parton * mother, Parton * daughter1, Parton * daughter2)
{
  Parton * swap;
  int oldc[2]={mother->GetFlow(1),mother->GetFlow(2)}, newc[2]={oldc[0],oldc[1]};
  if (mother->GetType()==pst::IS) {
    swap = mother; mother = daughter1; daughter1 = swap;
  }
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
	  SetColours(daughter1,oldc,newc);
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
	  SetColours(daughter1,oldc,newc);
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
	  SetColours(daughter1,oldc,newc);
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
	  SetColours(daughter1,oldc,newc);
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
	SetColours(daughter1,oldc,newc);
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
	  SetColours(mother,oldc,newc);
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
	  SetColours(mother,oldc,newc);
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
	  SetColours(mother,oldc,newc);
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
	  SetColours(mother,oldc,newc);
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
	SetColours(mother,oldc,newc);
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
	  Vec4D p((*plit)->Momentum());
 	  lt.Align(p);
	  (*plit)->SetMomentum(p);
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
	  Vec4D p((*plit)->Momentum());
	  lt.Align(p);
	  (*plit)->SetMomentum(p);
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
	  Vec4D p((*plit)->Momentum());
	  lt.Align(p);
	  (*plit)->SetMomentum(p);
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
	  Vec4D p((*plit)->Momentum());
  	  lt.Align(p);
	  (*plit)->SetMomentum(p);
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
	  Vec4D p((*plit)->Momentum());
	  lt.Align(p);
	  (*plit)->SetMomentum(p);
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
	  Vec4D p((*plit)->Momentum());
	  lt.Align(p);
 	  (*plit)->SetMomentum(p);
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
	  Vec4D p((*plit)->Momentum());
	  lt.Align(p);
	  (*plit)->SetMomentum(p);
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
	  Vec4D p((*plit)->Momentum());
	  lt.Align(p);
	  (*plit)->SetMomentum(p);
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
