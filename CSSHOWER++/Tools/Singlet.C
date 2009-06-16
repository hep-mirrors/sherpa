#include "CSSHOWER++/Tools/Singlet.H"
#include "CSSHOWER++/Tools/Parton.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <list>

using namespace CSSHOWER;
using namespace ATOOLS;
using namespace std;

std::ostream& CSSHOWER::operator<<(std::ostream& str, Singlet & singlet) {
  str<<"Singlet parton list from CS_Shower : "<<&singlet<<" Type "<<singlet.GetType()<<endl;
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

bool Singlet::CheckColours() {
  if (m_stype==stype::neutral) return true;
  int flow1((*begin())->GetFlow(1)), flow2((*begin())->GetFlow(2));
  for (PLiter plit=begin();plit!=end();plit++) { 
    if (plit==begin()) plit++;
    if ((*plit)->GetFlow(2)!=flow1) {
      if (m_stype!=stype::ring) return false;
      else {
	if ((*plit)->GetFlow(1)!=flow2) return false;
	else flow2=(*plit)->GetFlow(2);
      }
    }
    else flow1=(*plit)->GetFlow(1);
  }
  return true;
}

void Singlet::BoostAllFS(Parton *l,Parton *r,Parton *s,Parton *f,
			 const Flavour &mo,const int mode)
{
  if (p_all==NULL) return;
  if (mode&2) {
    if (mode&1) {
      Vec4D pr(r->Momentum()), pl(l->Momentum()), ps(s->Momentum());
      double papb(pl*ps), pipa(-pr*pl), pipb(-pr*ps);
      double xiab((papb+pipa+pipb)/papb), mi2(p_ms->Mass2(r->GetFlavour()));
      double mai2(p_ms->Mass2(mo)), ma2(p_ms->Mass2(l->GetFlavour()));
      double mb2(p_ms->Mass2(s->GetFlavour())), Q2((pl+ps-pr).Abs2());
      double ttau(Q2-mai2-mb2), tau(Q2-ma2-mi2-mb2);
      double xiiab(xiab*(ttau+sqrt(ttau*ttau-4.*mai2*mb2))/
		   (tau+sqrt(tau*tau-4.*ma2*mb2*xiab*xiab)));
      Vec4D pait(xiiab*(Q2*Q2-mai2*mb2)/(Q2*Q2-ma2*mb2)*(pl-ma2/Q2*ps)+mai2/Q2*ps);
      Vec4D K(pl+ps-pr), Kt(pait+ps), KKt(K+Kt);
      for (All_Singlets::const_iterator asit(p_all->begin());
	   asit!=p_all->end();++asit)
	for (PLiter plit((*asit)->begin());plit!=(*asit)->end();++plit)
	  if ((*plit)->GetType()!=pst::IS) {
	    Vec4D p((*plit)->Momentum());
	    p=p-2.0*(p*KKt)/KKt.Abs2()*KKt+2.0*(p*Kt)/K.Abs2()*K;
	    (*plit)->SetMomentum(p);
	  }
      r->SetMomentum(pr);
    }
    else {
      Parton *b(NULL);
      for (PLiter plit(s->GetSing()->begin());plit!=s->GetSing()->end();++plit)
	if ((*plit)->GetType()==pst::IS && *plit!=s) {
	  b=*plit;
	  break;
	}
      if (b==NULL) THROW(fatal_error,"Corrupted singlet");
      Vec4D pb(-b->Momentum()), ps(-s->Momentum());
      double ma2(p_ms->Mass2(s->GetFlavour()));
      double mb2(p_ms->Mass2(b->GetFlavour())), pspb(ps*pb);
      double del(pspb+sqrt(pspb*pspb-ma2*mb2));
      double pbm(0.5*(pb[0]-dabs(pb[3]))), sb(Sign(pb[3])), pap(0.25*del/pbm);
      Vec4D pm(pbm,0.,0.,-sb*pbm), pp(pap,0.,0.,sb*pap);
      Vec4D pan(pp+ma2/del*pm), pbn(pm+mb2/del*pp);
      Poincare cmso(-ps+b->Momentum());
      cmso.Boost(ps);
      Poincare zrot(ps,-sb*Vec4D::ZVEC);
      Poincare cmsn(-pan-pbn);
      for (All_Singlets::const_iterator asit(p_all->begin());
	   asit!=p_all->end();++asit) {
	Vec4D bsum, asum;
	for (PLiter plit((*asit)->begin());plit!=(*asit)->end();++plit) {
	  Vec4D p((*plit)->Momentum());
	  bsum+=(*plit)->GetType()==pst::IS?-p:p;
	  cmso.Boost(p);
	  zrot.Rotate(p);
	  cmsn.BoostBack(p);
	  (*plit)->SetMomentum(p);
	  asum+=(*plit)->GetType()==pst::IS?-p:p;
	}
	if (!IsEqual(bsum,asum,1.0e-6))
	  msg_Error()<<METHOD<<"(): Momentum not conserved FI:\n"
		     <<"  before "<<bsum<<"\n"
		     <<"  after  "<<asum<<"\n"<<**asit;
      }
    }
  }
  else {
    if (mode&1) {
      Parton *b(NULL);
      for (PLiter plit(f->GetSing()->begin());plit!=f->GetSing()->end();++plit)
	if ((*plit)->GetType()==pst::IS && *plit!=f) {
	  b=*plit;
	  break;
	}
      if (b==NULL) THROW(fatal_error,"Corrupted singlet");
      Vec4D pb(-b->Momentum()), pl(-l->Momentum());
      double ma2(p_ms->Mass2(l->GetFlavour()));
      double mb2(p_ms->Mass2(b->GetFlavour())), plpb(pl*pb);
      double del(plpb+sqrt(plpb*plpb-ma2*mb2));
      double pbm(0.5*(pb[0]-dabs(pb[3]))), sb(Sign(pb[3])), pap(0.25*del/pbm);
      Vec4D pm(pbm,0.,0.,-sb*pbm), pp(pap,0.,0.,sb*pap);
      Vec4D pan(pp+ma2/del*pm), pbn(pm+mb2/del*pp);
      Poincare cmso(-pl+b->Momentum());
      cmso.Boost(pl);
      Poincare zrot(pl,-sb*Vec4D::ZVEC);
      Poincare cmsn(-pan-pbn);
      for (All_Singlets::const_iterator asit(p_all->begin());
	   asit!=p_all->end();++asit) {
	Vec4D bsum, asum;
	for (PLiter plit((*asit)->begin());plit!=(*asit)->end();++plit) {
	  Vec4D p((*plit)->Momentum());
	  bsum+=(*plit)->GetType()==pst::IS?-p:p;
	  cmso.Boost(p);
	  zrot.Rotate(p);
	  cmsn.BoostBack(p);
	  (*plit)->SetMomentum(p);
	  asum+=(*plit)->GetType()==pst::IS?-p:p;
	}
	if (!IsEqual(bsum,asum,1.0e-6))
	  msg_Error()<<METHOD<<"(): Momentum not conserved IF:\n"
		     <<"  before "<<bsum<<"\n"
		     <<"  after  "<<asum<<"\n"<<**asit;
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
      Vec4D pr(r->Momentum()), pl(l->Momentum()), ps(s->Momentum());
      double papb(pl*ps), pipa(-pr*pl), pipb(-pr*ps);
      double xiab((papb+pipa+pipb)/papb), mi2(p_ms->Mass2(r->GetFlavour()));
      double mai2(p_ms->Mass2(mo)), ma2(p_ms->Mass2(l->GetFlavour()));
      double mb2(p_ms->Mass2(s->GetFlavour())), Q2((pl+ps-pr).Abs2());
      double ttau(Q2-mai2-mb2), tau(Q2-ma2-mi2-mb2);
      double xiiab(xiab*(ttau+sqrt(ttau*ttau-4.*mai2*mb2))/
		   (tau+sqrt(tau*tau-4.*ma2*mb2*xiab*xiab)));
      Vec4D pait(xiiab*(Q2*Q2-mai2*mb2)/(Q2*Q2-ma2*mb2)*(pl-ma2/Q2*ps)+mai2/Q2*ps);
      Vec4D K(pl+ps-pr), Kt(pait+ps), KKt(K+Kt);
      for (All_Singlets::const_iterator asit(p_all->begin());
	   asit!=p_all->end();++asit)
	for (PLiter plit((*asit)->begin());plit!=(*asit)->end();++plit)
	  if ((*plit)->GetType()!=pst::IS) {
	    Vec4D p((*plit)->Momentum());
	    p=p-2.0*(p*KKt)/KKt.Abs2()*KKt+2.0*(p*K)/K.Abs2()*Kt;
	    (*plit)->SetMomentum(p);
	  }
      r->SetMomentum(pr);
    }
    else {
      Parton *b(NULL);
      for (PLiter plit(s->GetSing()->begin());plit!=s->GetSing()->end();++plit)
	if ((*plit)->GetType()==pst::IS && *plit!=s) {
	  b=*plit;
	  break;
	}
      if (b==NULL) THROW(fatal_error,"Corrupted singlet");
      Vec4D pb(-b->Momentum()), ps(-s->Momentum());
      double ma2(p_ms->Mass2(s->GetFlavour()));
      double mb2(p_ms->Mass2(b->GetFlavour())), pspb(ps*pb);
      double del(pspb+sqrt(pspb*pspb-ma2*mb2));
      double pbm(0.5*(pb[0]-dabs(pb[3]))), sb(Sign(pb[3])), pap(0.25*del/pbm);
      Vec4D pm(pbm,0.,0.,-sb*pbm), pp(pap,0.,0.,sb*pap);
      Vec4D pan(pp+ma2/del*pm), pbn(pm+mb2/del*pp);
      Poincare cmso(-ps+b->Momentum());
      cmso.Boost(ps);
      Poincare zrot(ps,-sb*Vec4D::ZVEC);
      Poincare cmsn(-pan-pbn);
      for (All_Singlets::const_iterator asit(p_all->begin());
	   asit!=p_all->end();++asit) {
	Vec4D bsum, asum;
	for (PLiter plit((*asit)->begin());plit!=(*asit)->end();++plit) {
	  Vec4D p((*plit)->Momentum());
	  bsum+=(*plit)->GetType()==pst::IS?-p:p;
	  cmso.Boost(p);
	  zrot.Rotate(p);
	  cmsn.BoostBack(p);
 	  (*plit)->SetMomentum(p);
	  asum+=(*plit)->GetType()==pst::IS?-p:p;
	}
	if (!IsEqual(bsum,asum,1.0e-6))
	  msg_Error()<<METHOD<<"(): Momentum not conserved FI:\n"
		     <<"  before "<<bsum<<"\n"
		     <<"  after  "<<asum<<"\n"<<**asit;
      }
    }
  }
  else {
    if (mode&1) {
      Parton *b(NULL);
      for (PLiter plit(f->GetSing()->begin());plit!=f->GetSing()->end();++plit)
	if ((*plit)->GetType()==pst::IS && *plit!=f) {
	  b=*plit;
	  break;
	}
      if (b==NULL) THROW(fatal_error,"Corrupted singlet");
      Vec4D pb(-b->Momentum()), pl(-l->Momentum());
      double ma2(p_ms->Mass2(l->GetFlavour()));
      double mb2(p_ms->Mass2(b->GetFlavour())), plpb(pl*pb);
      double del(plpb+sqrt(plpb*plpb-ma2*mb2));
      double pbm(0.5*(pb[0]-dabs(pb[3]))), sb(Sign(pb[3])), pap(0.25*del/pbm);
      Vec4D pm(pbm,0.,0.,-sb*pbm), pp(pap,0.,0.,sb*pap);
      Vec4D pan(pp+ma2/del*pm), pbn(pm+mb2/del*pp);
      Poincare cmso(-pl+b->Momentum());
      cmso.Boost(pl);
      Poincare zrot(pl,-sb*Vec4D::ZVEC);
      Poincare cmsn(-pan-pbn);
      for (All_Singlets::const_iterator asit(p_all->begin());
	   asit!=p_all->end();++asit) {
	Vec4D bsum, asum;
	for (PLiter plit((*asit)->begin());plit!=(*asit)->end();++plit) {
	  Vec4D p((*plit)->Momentum());
	  bsum+=(*plit)->GetType()==pst::IS?-p:p;
	  cmso.Boost(p);
	  zrot.Rotate(p);
	  cmsn.BoostBack(p);
	  (*plit)->SetMomentum(p);
	  asum+=(*plit)->GetType()==pst::IS?-p:p;
	}
	if (!IsEqual(bsum,asum,1.0e-6))
	  msg_Error()<<METHOD<<"(): Momentum not conserved IF:\n"
		     <<"  before "<<bsum<<"\n"
		     <<"  after  "<<asum<<"\n"<<**asit;
      }
    }
  }
}

void Singlet::UpdateDaughters()
{
  for (PLiter plit(begin());plit!=end();++plit)
    (*plit)->UpdateDaughters();
}
