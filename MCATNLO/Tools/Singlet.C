#include "MCATNLO/Tools/Singlet.H"
#include "MCATNLO/Tools/Parton.H"
#include "PDF/Main/Shower_Base.H"
#include "ATOOLS/Math/ZAlign.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include <list>

using namespace MCATNLO;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

std::ostream& MCATNLO::operator<<(std::ostream& str, Singlet & singlet) {
  str<<"Singlet parton list from CS_MCATNLO : "<<&singlet<<endl;
  Parton * part;
  for (PLiter plit=singlet.begin();plit!=singlet.end();plit++) {
    part = (*plit);
    str<<(*part);
  }
  str<<"-------------------------------------------------------------------------"<<endl;
  return str;
}

std::ostream& MCATNLO::operator<<(std::ostream & str,All_Singlets & all) {
  str<<"Singlet list from CS_MCATNLO : "<<endl;
  Singlet * sing;
  for (ASiter asit=all.begin();asit!=all.end();asit++) {
    sing = (*asit);
    str<<sing<<" "<<sing->size()<<" "<<(*sing);
  }
  str<<"-------------------------------------------------------------------------"<<endl;
  return str;
}


Singlet::~Singlet() {
  for (Parton_List::const_iterator
	 pit(m_dels.begin());pit!=m_dels.end();++pit) delete *pit;
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
  Cluster_Amplitude *ampl(Cluster_Amplitude::New()); 
  ampl->SetJF(p_jf);
  size_t nin(0);
  for (const_iterator iit(begin());iit!=end();++iit)
    if ((*iit)->GetType()==pst::IS) ++nin;
  ampl->SetNIn(nin);
  for (const_iterator iit(begin());iit!=end();++iit)
    if ((*iit)->GetType()==pst::IS)
      ampl->CreateLeg(-(*iit)->Momentum(),
		      (*iit)->GetFlavour().Bar(),ColorID());
  for (const_iterator iit(begin());iit!=end();++iit)
    if ((*iit)->GetType()==pst::FS)
      ampl->CreateLeg((*iit)->Momentum(),
		      (*iit)->GetFlavour(),ColorID());
  bool veto(p_jf?p_shower->JetVeto(ampl):false);
  ampl->Delete();
  return veto;
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
  m_dels.push_back(mother);
  plit = erase(plit);
  if ((flav.IsGluon()  || flav.IsGluino()) && 
      (flav1.IsQuark() || flav1.IsSquark()) && 
      (flav2.IsQuark() || flav2.IsSquark())) { return 1; }
  return 0;
}

bool Singlet::ArrangeColours(Parton * mother, Parton * daughter1, Parton * daughter2)
{
  daughter1->SetSing(this);
  daughter2->SetSing(this);
  if (mother->GetType()==pst::IS) std::swap<Parton*>(mother,daughter1);
  Flavour mo(mother->GetFlavour()), d1(daughter1->GetFlavour()), d2(daughter2->GetFlavour());
  msg_Tracking()<<METHOD<<" for "<<mo<<" --> "<<d1<<" & "<<d2<<"\n";
  if (mother->GetType()==pst::FS) {
    if (mo.StrongCharge()==-3) {
      if (d1.StrongCharge()==-3) {
	if (d2.StrongCharge()==8) {
	  daughter2->SetFlow(2,mother->GetFlow(2));
	  daughter2->SetFlow(1,-1);
	  daughter1->SetFlow(2,daughter2->GetFlow(1));
	  return true;
	}
	else if (d2.StrongCharge()==0) {
	  daughter1->SetFlow(2,mother->GetFlow(2));
	  daughter2->SetFlow(1,0);
	  daughter2->SetFlow(2,0);
	  return true;
	}
      }
      else if (d2.StrongCharge()==-3) {
	if (d1.StrongCharge()==8) {
	  daughter1->SetFlow(2,mother->GetFlow(2));
	  daughter1->SetFlow(1,-1);
	  daughter2->SetFlow(2,daughter1->GetFlow(1));
	  return true;
	}
	else if (d1.StrongCharge()==0) {
	  daughter2->SetFlow(2,mother->GetFlow(2));
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
	  daughter2->SetFlow(2,-1);
	  daughter1->SetFlow(1,daughter2->GetFlow(2));
	  return true;
	}
	else if (d2.StrongCharge()==0) {
	  daughter1->SetFlow(1,mother->GetFlow(1));
	  daughter2->SetFlow(1,0);
	  daughter2->SetFlow(2,0);
	  return true;
	}
      }
      else if (d2.StrongCharge()==3) {
	if (d1.StrongCharge()==8) {
	  daughter1->SetFlow(1,mother->GetFlow(1));
	  daughter1->SetFlow(2,-1);
	  daughter2->SetFlow(1,daughter1->GetFlow(2));
	  return true;
	}
	else if (d1.StrongCharge()==0) {
	  daughter2->SetFlow(1,mother->GetFlow(1));
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
	daughter1->SetFlow(2,0);
	daughter2->SetFlow(1,0);
	daughter2->SetFlow(2,mother->GetFlow(2));
	return true;
      }
      else if (d1.StrongCharge()==-3 && 
	       d2.StrongCharge()==3) {  
	daughter2->SetFlow(1,mother->GetFlow(1));
	daughter2->SetFlow(2,0);
	daughter1->SetFlow(1,0);
	daughter1->SetFlow(2,mother->GetFlow(2));
	return true;
      }
      else if (d1.StrongCharge()==8 && 
	       d2.StrongCharge()==8) {
	if (mother->Col()<0) {
	  // if (mother->GetRight()==mother->GetSpect()) {
	  //   daughter2->SetFlow(1,mother->GetFlow(1));
	  //   daughter2->SetFlow(2,-1);
	  //   daughter1->SetFlow(1,daughter2->GetFlow(2));
	  //   daughter1->SetFlow(2,mother->GetFlow(2));
	  //   return true;
	  // }
	  daughter2->SetFlow(2,mother->GetFlow(2));
	  daughter2->SetFlow(1,-1);
	  daughter1->SetFlow(2,daughter2->GetFlow(1));
	  daughter1->SetFlow(1,mother->GetFlow(1));
	  return true;
	}
	// if (mother->GetRight()==mother->GetSpect()) {
	//   daughter1->SetFlow(1,mother->GetFlow(1));
	//   daughter1->SetFlow(2,-1);
	//   daughter2->SetFlow(1,daughter1->GetFlow(2));
	//   daughter2->SetFlow(2,mother->GetFlow(2));
	//   return true;
	// }
	daughter1->SetFlow(2,mother->GetFlow(2));
	daughter1->SetFlow(1,-1);
	daughter2->SetFlow(2,daughter1->GetFlow(1));
	daughter2->SetFlow(1,mother->GetFlow(1));
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
	  daughter2->SetFlow(2,mother->GetFlow(1));
	  return true;
	}
	else if (d2.StrongCharge()==0) {
	  mother->SetFlow(1,daughter1->GetFlow(1));
	  daughter2->SetFlow(1,0);
	  daughter2->SetFlow(2,0);
	  return true;
	}
      }
      else if (d2.StrongCharge()==3) {
	if (mo.StrongCharge()==8) {
	  mother->SetFlow(1,daughter1->GetFlow(1));
	  mother->SetFlow(2,-1);
	  daughter2->SetFlow(1,mother->GetFlow(2));
	  return true;
	}
	if (mo.StrongCharge()==0) {
	  daughter2->SetFlow(1,daughter1->GetFlow(1));
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
	  return true;
	}
	else if (d2.StrongCharge()==0) {
	  mother->SetFlow(2,daughter1->GetFlow(2));
	  daughter2->SetFlow(1,0);
	  daughter2->SetFlow(2,0);
	  return true;
	}
      }
      else if (d2.StrongCharge()==-3) {
	if (mo.StrongCharge()==8) {
	  mother->SetFlow(1,-1);
	  mother->SetFlow(2,daughter1->GetFlow(2));
	  daughter2->SetFlow(2,mother->GetFlow(1));
	  return true;
	}
	else if (mo.StrongCharge()==0) {
	  daughter2->SetFlow(2,daughter1->GetFlow(2));
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
	  mother->SetFlow(2,0);
	  mother->SetFlow(1,daughter1->GetFlow(1));
	  return true;
	}
	else if (d2.StrongCharge()==3) {
	  daughter2->SetFlow(2,0);
	  daughter2->SetFlow(1,daughter1->GetFlow(1));
	  mother->SetFlow(1,0);
	  mother->SetFlow(2,daughter1->GetFlow(2));
	  return true;
	}
      }
      else if (mo.StrongCharge()==8 && 
	       d2.StrongCharge()==8) {
	if (daughter1->Col()<0) {
	  // if (daughter1->GetRight()==daughter1->GetSpect()) {
	  //   daughter2->SetFlow(1,daughter1->GetFlow(1));
	  //   daughter2->SetFlow(2,-1);	
	  //   mother->SetFlow(2,daughter1->GetFlow(2));
	  //   mother->SetFlow(1,daughter2->GetFlow(2));
	  //   return true;
	  // }
	  daughter2->SetFlow(2,daughter1->GetFlow(2));
	  daughter2->SetFlow(1,-1);	
	  mother->SetFlow(1,daughter1->GetFlow(1));
	  mother->SetFlow(2,daughter2->GetFlow(1));
	  return true;
	}
	// if (daughter1->GetRight()==daughter1->GetSpect()) {
	//   mother->SetFlow(1,daughter1->GetFlow(1));
	//   mother->SetFlow(2,-1);	
	//   daughter2->SetFlow(2,daughter1->GetFlow(2));
	//   daughter2->SetFlow(1,mother->GetFlow(2));
	//   return true;
	// }
	mother->SetFlow(2,daughter1->GetFlow(2));
	mother->SetFlow(1,-1);	
	daughter2->SetFlow(1,daughter1->GetFlow(1));
	daughter2->SetFlow(2,mother->GetFlow(1));
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
  if (l->LT().empty()) return;
  for (All_Singlets::const_iterator asit(p_all->begin());
       asit!=p_all->end();++asit) {
    for (PLiter plit((*asit)->begin());plit!=(*asit)->end();++plit) {
      if (*plit==f || *plit==l || *plit==r || *plit==s) continue;
      (*plit)->SetMomentum(l->LT()*(*plit)->Momentum());
    }
  }
}

void Singlet::BoostBackAllFS(Parton *l,Parton *r,Parton *s,Parton *f,
			     const Flavour &mo,const int mode)
{
  if (p_all==NULL) return;
  Vec4D pa(l->Momentum()), pk(s->Momentum()), pi(r->Momentum());
  double ma2(p_ms->Mass2(l->GetFlavour())), mk2(p_ms->Mass2(s->GetFlavour()));
  double mi2(p_ms->Mass2(r->GetFlavour())), mai2(p_ms->Mass2(f->GetFlavour()));
  Kin_Args lp;
  if (mode&2) {
    if (mode&1) {
      lp=ClusterIIDipole(ma2,mi2,mai2,mk2,pa,pi,pk,2|(((mode&4)?f:l)->Kin()?4:0));
    }
  }
  else {
    if (mode&1) {
      Parton *b(NULL);
      for (PLiter pit(f->GetSing()->begin());pit!=f->GetSing()->end();++pit)
	if ((*pit)->GetType()==pst::IS && *pit!=f) {
	  b=*pit;
	  break;
	}
      if (b==NULL) THROW(fatal_error,"Corrupted singlet");
      double mb2(p_ms->Mass2(b->GetFlavour()));
      lp=ClusterIFDipole(ma2,mi2,mai2,mk2,mb2,pa,pi,pk,
			 b->Momentum(),2|(((mode&4)?f:l)->Kin()?4:0));
    }
  }
  if (lp.m_lam.empty()) return;
  for (All_Singlets::const_iterator asit(p_all->begin());
       asit!=p_all->end();++asit) {
    for (PLiter plit((*asit)->begin());plit!=(*asit)->end();++plit) {
      (*plit)->SetMomentum(lp.m_lam*(*plit)->Momentum());
    }
  }
}
