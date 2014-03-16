#include "MCATNLO/Tools/Singlet.H"
#include "MCATNLO/Tools/Parton.H"
#include "PDF/Main/Shower_Base.H"
#include "PDF/Main/Jet_Criterion.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
#include "ATOOLS/Math/ZAlign.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Random.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include <list>

using namespace MCATNLO;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

std::ostream& MCATNLO::operator<<(std::ostream& str, Singlet & singlet) {
  str<<"Singlet parton list from CS_MCatNLO : "<<&singlet<<endl;
  Parton * part;
  for (PLiter plit=singlet.begin();plit!=singlet.end();plit++) {
    part = (*plit);
    str<<(*part);
  }
  str<<"-------------------------------------------------------------------------"<<endl;
  return str;
}

std::ostream& MCATNLO::operator<<(std::ostream & str,All_Singlets & all) {
  str<<"Singlet list from CS_MCatNLO : "<<endl;
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
  bool veto(p_jf?p_jf->JC()->Jets(ampl):false);
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
  Flavour mo(mother->GetFlavour()), d1(daughter1->GetFlavour()), d2(daughter2->GetFlavour());
  msg_Tracking()<<METHOD<<" for "<<mo<<" --> "<<d1<<" & "<<d2<<"\n";
  if (mother->GetType()==pst::IS) { mo=mo.Bar(); d1=d1.Bar(); }
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
      int rc(mother->GetFlow(2) && mother->GetFlow(2) 
	     ==mother->GetSpect()->GetFlow(1));
      int lc(mother->GetFlow(1) && mother->GetFlow(1) 
	     ==mother->GetSpect()->GetFlow(2));
      if (rc==lc) {
	rc=mother->GetSpect()->GetFlow(1)?1:0;
	lc=mother->GetSpect()->GetFlow(2)?1:0;
	if (rc==lc) {
	  rc=ran->Get()>0.5?1:0;
	  lc=1-rc;
	}
      }
      if (rc==lc) THROW(fatal_error,"Internal error");
      if (mother->Col()<0) {
	if (rc) {
	  daughter2->SetFlow(1,mother->GetFlow(1));
	  daughter2->SetFlow(2,-1);
	  daughter1->SetFlow(1,daughter2->GetFlow(2));
	  daughter1->SetFlow(2,mother->GetFlow(2));
	  return true;
	}
	else {
	  daughter2->SetFlow(2,mother->GetFlow(2));
	  daughter2->SetFlow(1,-1);
	  daughter1->SetFlow(2,daughter2->GetFlow(1));
	  daughter1->SetFlow(1,mother->GetFlow(1));
	  return true;
	}
      }
      else {
	if (rc) {
	  daughter1->SetFlow(1,mother->GetFlow(1));
	  daughter1->SetFlow(2,-1);
	  daughter2->SetFlow(1,daughter1->GetFlow(2));
	  daughter2->SetFlow(2,mother->GetFlow(2));
	  return true;
	}
	else {
	  daughter1->SetFlow(2,mother->GetFlow(2));
	  daughter1->SetFlow(1,-1);
	  daughter2->SetFlow(2,daughter1->GetFlow(1));
	  daughter2->SetFlow(1,mother->GetFlow(1));
	  return true;
	}
      }
    }
  }
  else if (mo.StrongCharge()==0) {
    if (d1.StrongCharge()==3 &&
	d2.StrongCharge()==-3) {
      daughter1->SetFlow(1,-1);
      daughter1->SetFlow(2,0);
      daughter2->SetFlow(2,daughter1->GetFlow(1));
      daughter2->SetFlow(1,0);
      return true;
    }
    else if (d1.StrongCharge()==-3 &&
	     d2.StrongCharge()==3) {
      daughter2->SetFlow(1,-1);
      daughter2->SetFlow(2,0);
      daughter1->SetFlow(2,daughter2->GetFlow(1));
      daughter1->SetFlow(1,0);
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
  return false;
} 

void Singlet::BoostAllFS(Parton *l,Parton *r,Parton *s)
{
  if (l->LT().empty()) return;
  for (All_Singlets::const_iterator asit(p_all->begin());
       asit!=p_all->end();++asit) {
    for (PLiter plit((*asit)->begin());plit!=(*asit)->end();++plit) {
      if (*plit==l || *plit==r || *plit==s) continue;
      (*plit)->SetMomentum(l->LT()*(*plit)->Momentum());
    }
  }
}

void Singlet::BoostBackAllFS(Parton *l,Parton *r,Parton *s)
{
  if (p_all==NULL) return;
  Poincare_Sequence lt(l->LT());
  if (lt.size()) lt.Invert();
  if (lt.empty()) return;
  for (All_Singlets::const_iterator asit(p_all->begin());
       asit!=p_all->end();++asit) {
    for (PLiter plit((*asit)->begin());plit!=(*asit)->end();++plit) {
      (*plit)->SetMomentum(lt*(*plit)->Momentum());
    }
  }
}
