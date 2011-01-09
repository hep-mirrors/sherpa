#include "POWHEG/Tools/Singlet.H"
#include "POWHEG/Tools/Parton.H"
#include "POWHEG/Showers/Sudakov.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
#include "ATOOLS/Math/ZAlign.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <list>

using namespace POWHEG;
using namespace ATOOLS;
using namespace std;

std::ostream& POWHEG::operator<<(std::ostream& str, Singlet & singlet) {
  str<<"Singlet parton list from CS_POWHEG : "<<&singlet<<endl;
  Parton * part;
  for (PLiter plit=singlet.begin();plit!=singlet.end();plit++) {
    part = (*plit);
    str<<(*part);
  }
  str<<"-------------------------------------------------------------------------"<<endl;
  return str;
}

std::ostream& POWHEG::operator<<(std::ostream & str,All_Singlets & all) {
  str<<"Singlet list from CS_POWHEG : "<<endl;
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

bool Singlet::JetVeto(Sudakov *const sud) const
{
  if (p_jf==NULL) return false;
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
	  return true;
	}
	else if (d2.StrongCharge()==0) {
	  daughter1->SetFlow(2,mother->GetFlow(2));
	  daughter1->SetMEFlow(2,mother->GetMEFlow(2));
	  daughter2->SetFlow(1,0);
	  daughter2->SetFlow(2,0);
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
	  return true;
	}
	else if (d1.StrongCharge()==0) {
	  daughter2->SetFlow(2,mother->GetFlow(2));
	  daughter2->SetMEFlow(2,mother->GetMEFlow(2));
	  daughter1->SetFlow(1,0);
	  daughter1->SetFlow(2,0);
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
	  return true;
	}
	else if (d2.StrongCharge()==0) {
	  daughter1->SetFlow(1,mother->GetFlow(1));
	  daughter1->SetMEFlow(1,mother->GetMEFlow(1));
	  daughter2->SetFlow(1,0);
	  daughter2->SetFlow(2,0);
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
	  return true;
	}
	else if (d1.StrongCharge()==0) {
	  daughter2->SetFlow(1,mother->GetFlow(1));
	  daughter2->SetMEFlow(1,mother->GetMEFlow(1));
	  daughter1->SetFlow(1,0);
	  daughter1->SetFlow(2,0);
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
	  return true;
	}
	else if (d2.StrongCharge()==0) {
	  mother->SetFlow(1,daughter1->GetFlow(1));
	  mother->SetMEFlow(1,daughter1->GetMEFlow(1));
	  daughter2->SetFlow(1,0);
	  daughter2->SetFlow(2,0);
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
	  return true;
	}
	if (mo.StrongCharge()==0) {
	  daughter2->SetFlow(1,daughter1->GetFlow(1));
	  daughter2->SetMEFlow(1,daughter1->GetMEFlow(1));
	  mother->SetFlow(1,0);
	  mother->SetFlow(2,0);
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
	  return true;
	}
	else if (d2.StrongCharge()==0) {
	  mother->SetFlow(2,daughter1->GetFlow(2));
	  mother->SetMEFlow(2,daughter1->GetMEFlow(2));
	  daughter2->SetFlow(1,0);
	  daughter2->SetFlow(2,0);
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
	  return true;
	}
	else if (mo.StrongCharge()==0) {
	  daughter2->SetFlow(2,daughter1->GetFlow(2));
	  daughter2->SetMEFlow(2,daughter1->GetMEFlow(2));
	  mother->SetFlow(1,0);
	  mother->SetFlow(2,0);
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
	  return true;
	}
	else if (d2.StrongCharge()==3) {
	  daughter2->SetFlow(2,0);
	  daughter2->SetFlow(1,daughter1->GetFlow(1));
	  daughter2->SetMEFlow(1,daughter1->GetMEFlow(1));
	  mother->SetFlow(1,0);
	  mother->SetFlow(2,daughter1->GetFlow(2));
	  mother->SetMEFlow(2,daughter1->GetMEFlow(2));
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
	  return true;
	}
	else if (d2.StrongCharge()==3) {
	  daughter2->SetFlow(2,0);
	  daughter2->SetFlow(1,-1);
	  mother->SetFlow(2,daughter2->GetFlow(1));
	  mother->SetFlow(1,0);
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
      Vec4D pr(r->Momentum()), pl(l->Momentum()), ps(s->Momentum());
      double papb(pl*ps), pipa(-pr*pl), pipb(-pr*ps);
      double xiab((papb+pipa+pipb)/papb), mi2(p_ms->Mass2(r->GetFlavour()));
      double mai2(p_ms->Mass2(mo)), ma2(p_ms->Mass2(l->GetFlavour()));
      double mb2(p_ms->Mass2(s->GetFlavour())), Q2((pl+ps-pr).Abs2());
      double ttau(Q2-mai2-mb2), tau(Q2-ma2-mi2-mb2);
      double xiiab(xiab*(ttau+sqrt(ttau*ttau-4.*mai2*mb2))/
		   (tau+sqrt(tau*tau-4.*ma2*mb2*xiab*xiab)));
      double gam(papb+sqrt(papb*papb-ma2*mb2)), sb(Sign(ps[3]));
      Vec4D pait(xiiab*(1.0-mai2*mb2/sqr(gam*xiiab))/(1.0-ma2*mb2/sqr(gam))
		 *(pl-ma2/gam*ps)+mai2/(xiiab*gam)*ps);
      Poincare cmso(pait+ps), cmsn(pl+ps-pr);
      cmsn.Boost(pait);
      Poincare zrot(pait,-sb*Vec4D::ZVEC);
      for (All_Singlets::const_iterator asit(p_all->begin());
	   asit!=p_all->end();++asit) {
	for (PLiter plit((*asit)->begin());plit!=(*asit)->end();++plit) {
	  if ((*plit)==f || (*plit)==r || (*plit)==s) continue;
	  Vec4D p((*plit)->Momentum());
	  cmso.Boost(p);
	  zrot.RotateBack(p);
	  cmsn.BoostBack(p);
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
      Vec4D pr(r->Momentum()), pl(l->Momentum()), ps(s->Momentum());
      double papb(pl*ps), pipa(-pr*pl), pipb(-pr*ps);
      double xiab((papb+pipa+pipb)/papb), mi2(p_ms->Mass2(r->GetFlavour()));
      double mai2(p_ms->Mass2(mo)), ma2(p_ms->Mass2(l->GetFlavour()));
      double mb2(p_ms->Mass2(s->GetFlavour())), Q2((pl+ps-pr).Abs2());
      double ttau(Q2-mai2-mb2), tau(Q2-ma2-mi2-mb2);
      double xiiab(xiab*(ttau+sqrt(ttau*ttau-4.*mai2*mb2))/
		   (tau+sqrt(tau*tau-4.*ma2*mb2*xiab*xiab)));
      double gam(papb+sqrt(papb*papb-ma2*mb2)), sb(Sign(ps[3]));
      Vec4D pait(xiiab*(1.0-mai2*mb2/sqr(gam*xiiab))/(1.0-ma2*mb2/sqr(gam))
		 *(pl-ma2/gam*ps)+mai2/(xiiab*gam)*ps);
      Poincare cmso(pl+ps-pr), cmsn(pait+ps);
      cmso.Boost(pait);
      Poincare zrot(pait,-sb*Vec4D::ZVEC);
      for (All_Singlets::const_iterator asit(p_all->begin());
	   asit!=p_all->end();++asit) {
	for (PLiter plit((*asit)->begin());plit!=(*asit)->end();++plit) {
	  Vec4D p((*plit)->Momentum());
	  cmso.Boost(p);
	  zrot.Rotate(p);
	  cmsn.BoostBack(p);
 	  (*plit)->SetMomentum(p);
	}
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
