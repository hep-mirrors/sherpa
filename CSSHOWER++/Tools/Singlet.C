#include "CSSHOWER++/Tools/Singlet.H"
#include "CSSHOWER++/Tools/Parton.H"
#include "CSSHOWER++/Showers/Sudakov.H"
#include "CSSHOWER++/Showers/Shower.H"
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
  str<<"Singlet parton list from CS_Shower:"<<endl;
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
  str<<"k_T,next = "<<sqrt(singlet.KtNext())
     <<", mu_R = "<<sqrt(singlet.MuR2())
     <<", nlo = "<<singlet.NLO()<<", nmax = "<<singlet.NMax()
     <<", K = "<<singlet.LKF()<<"\n";
  str<<"mom sum "<<sum<<"\n";
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


Singlet::~Singlet()
{
  for (PLiter plit(begin());plit!=end();++plit) delete *plit;
}

Parton *Singlet::IdParton(const size_t &id) const
{
  for (const_iterator it(begin());it!=end();++it)
    if ((*it)->Id()==id) return *it;
  return NULL;
}

double Singlet::JetVeto(Sudakov *const sud) const
{
  DEBUG_FUNC("jf = "<<p_jf);
  if (p_jf==NULL) return 0.0;
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
  double jcv(p_jf->JC()->Value(ampl));
  ampl->Delete();
  msg_Debugging()<<"--- "<<jcv<<" ---\n";
  return jcv;
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
  if (mother==p_split) p_split=part1;
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
	    if ((*plit)->FromDec()==0 && (*plit)->Id()==(1<<j)) {
	      part->SetOriginalPart(blob->InParticle(i));
	    }
	    ++j;
	  }
	}
      }
    }
    part->SetFromDec((*plit)->FromDec());
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
  daughter1->SetPrev(mother);
  daughter1->UpdateColours(mother->GetFlow(1),mother->GetFlow(2));
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
  int newr(daughter1->GetFlow(1)), newa(daughter1->GetFlow(2));
  daughter1->SetFlow(1,mother->GetFlow(1));
  daughter1->SetFlow(2,mother->GetFlow(2));
  daughter1->UpdateColours(newr,newa);
  daughter2->UpdateColours(daughter2->GetFlow(1),daughter2->GetFlow(2));
  for (iterator pit(begin());pit!=end();++pit)
    if (*pit==daughter1) *pit=mother;
  return true;
}

void Singlet::BoostAllFS(Parton *l,Parton *r,Parton *s,bool onlyFS)
{
  if (l->LT().empty()) return;
  for (PLiter plit(begin());plit!=end();++plit) {
    if (onlyFS && ((*plit)->GetType()!=pst::FS || (*plit)==r)) continue;
    Vec4D p(l->LT()*(*plit)->Momentum());
    if ((*plit)->GetType()==pst::IS &&
	IsZero(p.PPerp2())) p[1]=p[2]=0.0;
    if ((*plit)->Mass2()==0.0) p[0]=p.PSpat();
    (*plit)->SetMomentum(p);
  }
}


void Singlet::BoostBackAllFS(Parton *l,Parton *r,Parton *s,bool onlyFS)
{
  if (p_all==NULL) return;
  Poincare_Sequence lt(l->LT());
  if (lt.size()) lt.Invert();
  if (lt.empty()) return;
  for (PLiter plit(begin());plit!=end();++plit) {
    if (onlyFS && ((*plit)->GetType()!=pst::FS || (*plit)==r)) continue;
    Vec4D p(lt*(*plit)->Momentum());
    if ((*plit)->GetType()==pst::IS &&
	IsZero(p.PPerp2())) p[1]=p[2]=0.0;
    if ((*plit)->Mass2()==0.0) p[0]=p.PSpat();
    (*plit)->SetMomentum(p);
  }
}

void Singlet::Reduce()
{
  if (p_split==NULL) return;
  for (const_iterator it(begin());it!=end();++it) {
    (*it)->SetNext(NULL);
    (*it)->SetStat((*it)->Stat()&~1);
    if ((*it)->Stat()) THROW(fatal_error,"Cannot reduce singlet");
  }
  Singlet *next(p_left->GetSing());
  while (next) {
    bool found(false);
    for (All_Singlets::iterator asit(p_all->begin());
	 asit!=p_all->end();++asit)
      if (*asit==next) {
	Parton *split(next->GetSplit());
	if (split && (split->Stat()&2)) {
	  bool fdec(false);
	  for (iterator cit(begin());cit!=end();++cit)
	    if ((*cit)->Id()==split->Id()) {
	      split->SetPrev(*cit);
	      double jcv(0.0);
	      p_shower->ReconstructDaughters(next,jcv,NULL,NULL);
	      Parton *left(new Parton(next->GetLeft()->GetFlavour(),
				      next->GetLeft()->Momentum(),
				      next->GetLeft()->GetType()));
	      *left=*next->GetLeft();
	      Parton *right(new Parton(next->GetRight()->GetFlavour(),
				       next->GetRight()->Momentum(),
				       next->GetRight()->GetType()));
	      *right=*next->GetRight();
	      SplitParton(*cit,left,right);
	      fdec=true;
	      break;
	    }
	  if (!fdec) THROW(fatal_error,"Invalid tree structure");
	}
	next=next->GetLeft()?next->GetLeft()->GetSing():NULL;
	delete *asit;
	p_all->erase(asit);
	found=true;
	break;
      }
    if (!found) THROW(fatal_error,"Invalid tree structure");
  }
  p_left=p_right=p_split=p_spec=NULL;
  m_kt2_next=0.0;
}
