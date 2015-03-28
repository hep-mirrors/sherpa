#include "CSSHOWER++/Tools/Singlet.H"
#include "CSSHOWER++/Tools/Parton.H"
#include "CSSHOWER++/Showers/Sudakov.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
#include "PDF/Main/Jet_Criterion.H"
#include "ATOOLS/Math/ZAlign.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_Limits.H"
#include <list>

using namespace CSSHOWER;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

std::ostream& CSSHOWER::operator<<(std::ostream& str, Singlet & singlet) {
  Vec4D sum;
  str<<"Singlet parton list from CS_Shower : "<<&singlet<<", jf = "<<singlet.JF()<<endl;
  Parton * part;
  for (PLiter plit=singlet.begin();plit!=singlet.end();plit++) {
    part = (*plit);
    sum+=part->GetType()==pst::IS?-part->Momentum():part->Momentum();
    str<<(*part);
  }
  if (singlet.GetSplit() || singlet.GetLeft() || singlet.GetRight() || singlet.GetSpec()) {
    if (singlet.GetSplit()) str<<"Split:  "<<singlet.GetSplit()<<"  ";
    if (singlet.GetLeft()) str<<"Left:  "<<singlet.GetLeft()<<"  ";
    if (singlet.GetRight()) str<<"Right:  "<<singlet.GetRight()<<"  ";
    if (singlet.GetSpec()) str<<"Spec:  "<<singlet.GetSpec()<<"  ";
    str<<"\n";
  }
  str<<"mom sum "<<sum
     <<", k_T,next = "<<sqrt(singlet.KtNext())
     <<", nlo = "<<singlet.NLO()<<", K = "<<singlet.LKF()<<"\n";
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
  Cluster_Amplitude *ampl(Cluster_Amplitude::New());
  for (const_iterator iit(begin());iit!=end();++iit) {
    if ((*iit)->GetType()==pst::FS) continue;
    ampl->CreateLeg(-(*iit)->Momentum(),(*iit)->GetFlavour().Bar(),
		    ColorID((*iit)->GetFlow(1),(*iit)->GetFlow(2)),
		    1<<ampl->Legs().size());
  }
  ampl->SetNIn(ampl->Legs().size());
  for (const_iterator iit(begin());iit!=end();++iit) {
    if ((*iit)->GetType()==pst::IS) continue;
    ampl->CreateLeg((*iit)->Momentum(),(*iit)->GetFlavour(),
		    ColorID((*iit)->GetFlow(1),(*iit)->GetFlow(2)),
		    1<<ampl->Legs().size());
  }
  ampl->SetJF(p_jf);
  ampl->SetMS(p_ms);
  ampl->Decays()=m_decs;
  bool res(p_jf->JC()->Jets(ampl));
  ampl->Delete();
  if (res) msg_Debugging()<<"--- Jet veto ---\n";
  return res;
}

int Singlet::SplitParton(Parton * mother, Parton * part1, Parton * part2) 
{
  iterator plit(begin());
  for (;plit!=end();++plit) if (*plit==mother) break;
  if (plit==end()) THROW(fatal_error,"Internal error");

  if (mother->GetLeft()==part1->GetLeft()) part1->SetSoft(0,mother->KtSoft(0));
  else if (mother->GetLeft()==part2->GetLeft()) part2->SetSoft(0,mother->KtSoft(0));
  if (mother->GetRight()==part1->GetRight()) part1->SetSoft(1,mother->KtSoft(1));
  else if (mother->GetRight()==part2->GetRight()) part2->SetSoft(1,mother->KtSoft(1));

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

  plit++;
  delete mother; 
  plit = erase(plit);
  if (flav.StrongCharge()==8 && 
      abs(flav1.StrongCharge())==3 && 
      abs(flav2.StrongCharge())==3) { return 1; }
  return 0;
}

void Singlet::ExtractPartons
(ATOOLS::Blob * blob,ATOOLS::Mass_Selector *const ms) 
{
  Particle * part;
  for (PLiter plit=begin();plit!=end();plit++) {
    if ((*plit)->Stat()&1) continue;
    part = new Particle(-1,(*plit)->GetFlavour(),(*plit)->Momentum(),'F');
    part->SetNumber(0);
    for (size_t i(0);i<m_decs.size();++i)
      if (m_decs[i]->m_id&(*plit)->Id()) part->SetMEId((*plit)->Id());
    if ((*plit)->GetType()==pst::IS) {
      part->SetBeam((*plit)->Beam());
      part->SetInfo('I');
      blob->AddToInParticles(part);
    } 
    else {
      blob->AddToOutParticles(part);
      if (rpa->gen.SoftSC()) {
	size_t j=2;
	for (size_t i=0; i<blob->NInP(); ++i) {
	  if (blob->InParticle(i)->ProductionBlob() &&
	      blob->InParticle(i)->ProductionBlob()->Type()!=btp::Beam) {
	    if ((*plit)->Id()==(1<<j)) {
	      part->SetOriginalPart(blob->InParticle(i));
	    }
	    ++j;
	  }
	}
      }
    }
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
      np->SetMass2(p->Mass2());
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
  daughter1->SetSing(this);
  for (iterator pit(begin());pit!=end();++pit)
    if (*pit==mother) {
      *pit=daughter1;
      break;
    }
  daughter1->SetFlow(1,mother->GetFlow(1));
  daughter1->SetFlow(2,mother->GetFlow(2));
  daughter1->SetPrev(mother);
  daughter1->UpdateColours();
  daughter1->SetLeftOf(mother);
  daughter1->SetRightOf(mother);
  for (iterator pit(begin());pit!=end();++pit)
    if (*pit==daughter1) *pit=mother;
  return true;
}


bool Singlet::ArrangeColours(Parton * mother, Parton * daughter1, Parton * daughter2)
{
  daughter1->SetSing(this);
  daughter2->SetSing(this);
  for (iterator pit(begin());pit!=end();++pit)
    if (*pit==mother) {
      *pit=daughter1;
      break;
    }
  daughter1->SetMEFlow(1,mother->GetFlow(1));
  daughter1->SetMEFlow(2,mother->GetFlow(2));
  daughter1->SetPrev(mother);
  daughter2->SetFlow(1,0);
  daughter2->SetFlow(2,0);
  Flavour mo(mother->GetFlavour()), d1(daughter1->GetFlavour()), d2(daughter2->GetFlavour());
  if (mother->GetType()==pst::IS) { mo=mo.Bar(); d1=d1.Bar(); }
  if (mo.StrongCharge()==-3) {
    if (d1.StrongCharge()==-3) {
      if (d2.StrongCharge()==8) {
	daughter2->SetFlow(2,mother->GetFlow(2));
	daughter2->SetFlow(1,-1);
	daughter1->SetFlow(2,daughter2->GetFlow(1));
      }
      else if (d2.StrongCharge()==0) {
	daughter1->SetFlow(2,mother->GetFlow(2));
      }
    }
    else if (d2.StrongCharge()==-3) {
      if (d1.StrongCharge()==8) {
	daughter1->SetFlow(2,mother->GetFlow(2));
	daughter1->SetFlow(1,-1);
	daughter2->SetFlow(2,daughter1->GetFlow(1));
      }
      else if (d1.StrongCharge()==0) {
	daughter2->SetFlow(2,mother->GetFlow(2));
      }
    }
  }
  else if (mo.StrongCharge()==3) {
    if (d1.StrongCharge()==3) {
      if (d2.StrongCharge()==8) {
	daughter2->SetFlow(1,mother->GetFlow(1));
	daughter2->SetFlow(2,-1);
	daughter1->SetFlow(1,daughter2->GetFlow(2));
      }
      else if (d2.StrongCharge()==0) {
	daughter1->SetFlow(1,mother->GetFlow(1));
      }
    }
    else if (d2.StrongCharge()==3) {
      if (d1.StrongCharge()==8) {
	daughter1->SetFlow(1,mother->GetFlow(1));
	daughter1->SetFlow(2,-1);
	daughter2->SetFlow(1,daughter1->GetFlow(2));
      }
      else if (d1.StrongCharge()==0) {
	daughter2->SetFlow(1,mother->GetFlow(1));
      }
    }
  }
  else if (mo.StrongCharge()==8) {
    if (d1.StrongCharge()==3 && 
	d2.StrongCharge()==-3) {  
      daughter1->SetFlow(1,mother->GetFlow(1));
      daughter1->SetFlow(2,0);
      daughter2->SetFlow(2,mother->GetFlow(2));
    }
    else if (d1.StrongCharge()==-3 && 
	     d2.StrongCharge()==3) {  
      daughter2->SetFlow(1,mother->GetFlow(1));
      daughter1->SetFlow(1,0);
      daughter1->SetFlow(2,mother->GetFlow(2));
    }
    else if (d1.StrongCharge()==8 && 
	     d2.StrongCharge()==8) {
      if (mother->Col()<0) {
	if (mother->GetRight()==mother->GetSpect()) {
	  daughter2->SetFlow(1,mother->GetFlow(1));
	  daughter2->SetFlow(2,-1);
	  daughter1->SetFlow(1,daughter2->GetFlow(2));
	  daughter1->SetFlow(2,mother->GetFlow(2));
	}
	else {
	daughter2->SetFlow(2,mother->GetFlow(2));
	daughter2->SetFlow(1,-1);
	daughter1->SetFlow(2,daughter2->GetFlow(1));
	daughter1->SetFlow(1,mother->GetFlow(1));
	}
      }
      else {
      if (mother->GetRight()==mother->GetSpect()) {
	daughter1->SetFlow(1,mother->GetFlow(1));
	daughter1->SetFlow(2,-1);
	daughter2->SetFlow(1,daughter1->GetFlow(2));
	daughter2->SetFlow(2,mother->GetFlow(2));
      }
      else {
      daughter1->SetFlow(2,mother->GetFlow(2));
      daughter1->SetFlow(1,-1);
      daughter2->SetFlow(2,daughter1->GetFlow(1));
      daughter2->SetFlow(1,mother->GetFlow(1));
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
    }
    else if (d1.StrongCharge()==-3 &&
	     d2.StrongCharge()==3) {
      daughter2->SetFlow(1,-1);
      daughter1->SetFlow(1,0);
      daughter1->SetFlow(2,daughter2->GetFlow(1));
    }
    else if (d1.StrongCharge()==0 && 
	     d2.StrongCharge()==0) {
      daughter1->SetFlow(1,0);
      daughter1->SetFlow(2,0);
    }
  }
  daughter1->UpdateColours();
  daughter2->UpdateColours();
  for (iterator pit(begin());pit!=end();++pit)
    if (*pit==daughter1) *pit=mother;
  return true;
} 

void Singlet::BoostAllFS(Parton *l,Parton *r,Parton *s)
{
  if (l->LT().empty()) return;
  for (All_Singlets::const_iterator asit(p_all->begin());
       asit!=p_all->end();++asit) {
    for (PLiter plit((*asit)->begin());plit!=(*asit)->end();++plit) {
      if ((*plit)->FixSpec()!=Vec4D()) {
	(*plit)->SetFixSpec(l->LT()*(*plit)->FixSpec());
	(*plit)->SetOldMomentum(l->LT()*(*plit)->OldMomentum());
      }
      Vec4D p(l->LT()*(*plit)->Momentum());
      if ((*plit)->GetType()==pst::IS &&
	  IsZero(p.PPerp2())) p[1]=p[2]=0.0;
      if ((*plit)->Mass2()==0.0) p[0]=p.PSpat();
      (*plit)->SetMomentum(p);
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
      Vec4D p(lt*(*plit)->Momentum());
      if ((*plit)->GetType()==pst::IS &&
	  IsZero(p.PPerp2())) p[1]=p[2]=0.0;
      if ((*plit)->Mass2()==0.0) p[0]=p.PSpat();
      (*plit)->SetMomentum(p);
      if ((*plit)->FixSpec()!=Vec4D()) {
	(*plit)->SetFixSpec(lt*(*plit)->FixSpec());
	(*plit)->SetOldMomentum(lt*(*plit)->OldMomentum());
      }
    }
  }
}

void Singlet::UpdateDaughters()
{
  for (PLiter plit(begin());plit!=end();++plit)
    (*plit)->UpdateDaughters();
}
