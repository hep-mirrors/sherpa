#include "CSSHOWER++/Showers/Shower.H"
#include "CSSHOWER++/Tools/Parton.H"
#include "PDF/Remnant/Remnant_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Exception.H"

using namespace CSSHOWER;
using namespace ATOOLS;
using namespace std;

Shower::Shower(PDF::ISR_Handler * isr,const int qed,
	       Data_Reader *const dataread) : 
  p_actual(NULL), m_sudakov(isr,qed), p_isr(isr)
{
  m_sudakov.InitSplittingFunctions(MODEL::s_model);
  double k0sq   = dataread->GetValue<double>("CSS_PT2MIN",1.0);
  double is_fac = dataread->GetValue<double>("CSS_AS_IS_FAC",0.25);
  double fs_fac = dataread->GetValue<double>("CSS_AS_FS_FAC",1.0);
  m_sudakov.SetCoupling(MODEL::s_model,k0sq,is_fac,fs_fac);
  m_sudakov.SetShower(this);
}

Shower::~Shower() 
{
}

bool Shower::EvolveShower(Singlet * actual,const size_t &maxem,size_t &nem)
{
  return EvolveSinglet(actual,maxem,nem);
}

int Shower::RemnantTest(Parton *const p)
{
  if (p->Momentum()[0]>rpa.gen.PBeam(p->Beam())[0]) return -1;
  return p_isr->GetRemnant(p->Beam())->
    TestExtract(p->GetFlavour(),p->Momentum())?1:-1;
}

bool Shower::ReconstructDaughters(Singlet *const split,const bool one)
{
  if (split==NULL || split->GetLeft()==NULL) return true;
  if (split->GetRight()==NULL) THROW(fatal_error,"Invalid tree structure");
  msg_Debugging()<<METHOD<<"("<<split<<"): {\n";
  msg_Indent();
  Parton *l(split->GetLeft()), *r(split->GetRight());
  Parton *c(split->GetSplit()->FollowUp()), *s(split->GetSpec());
  Flavour fli(l->GetFlavour());
  s->SetMomentum(s->GetPrev()->Momentum());
  l->SetMomentum(c->Momentum());
  l->SetFlavour(c->GetFlavour());
  l->SetSpect(s);
  msg_Debugging()<<"before: c: "<<*l<<"        s: "<<*s<<"\n";
  msg_Debugging()<<"kt = "<<sqrt(l->KtTest())<<", z = "
		 <<l->ZTest()<<", y = "<<l->YTest()
		 <<", phi = "<<l->Phi()<<"\n\n";
  int stat=0;
  if (c->GetType()==pst::FS) {
    if (s->GetPrev()->GetType()==pst::FS) {
      m_kinFF.SetJF(NULL);
      stat=m_kinFF.MakeKinematics(l,fli,r->GetFlavour(),r);
    }
    else {
      m_kinFI.SetJF(NULL);
      stat=m_kinFI.MakeKinematics(l,fli,r->GetFlavour(),r);
      if (stat>0) {
	split->BoostAllFS(l,r,s,split->GetSplit(),
			  split->GetSplit()->GetFlavour(),2);
	stat=RemnantTest(s);
	if (stat<=0) split->BoostBackAllFS
	  (l,r,s,split->GetSplit(),split->GetSplit()->GetFlavour(),2);
      }
      s->SetXbj(s->GetPrev()->Xbj()/(1.0-l->YTest()));
    }
  }
  else {
    if (s->GetPrev()->GetType()==pst::FS) {
      m_kinIF.SetJF(NULL);
      stat=m_kinIF.MakeKinematics(l,fli,r->GetFlavour(),r);
      if (stat>0) {
	split->BoostAllFS(l,r,s,split->GetSplit(),
			  split->GetSplit()->GetFlavour(),1);
	stat=RemnantTest(l);
	if (stat<=0) split->BoostBackAllFS
	  (l,r,s,split->GetSplit(),split->GetSplit()->GetFlavour(),1);
      }
    }
    else {
      m_kinII.SetJF(NULL);
      stat=m_kinII.MakeKinematics(l,fli,r->GetFlavour(),r);
      if (stat>0) {
	split->BoostAllFS(l,r,s,split->GetSplit(),
			  split->GetSplit()->GetFlavour(),3);
	stat=RemnantTest(l);
	if (stat<=0) split->BoostBackAllFS
	  (l,r,s,split->GetSplit(),split->GetSplit()->GetFlavour(),3);
      }
    }
    l->SetXbj(c->Xbj()/l->ZTest());
  }
  l->SetFlavour(fli);
  msg_Debugging()<<"after: l: "<<*l<<"       r: "<<*r<<"       s: "<<*s<<"\n";
  if (stat<0) return false;
  l->GetSing()->UpdateDaughters();
  if(l->GetSing()!=r->GetSing()) r->GetSing()->UpdateDaughters();
  if (one) return true;
  bool nres(ReconstructDaughters(l->GetSing()));
  if (nres && l->GetSing()!=r->GetSing())
    nres=ReconstructDaughters(r->GetSing());
  split->BoostBackAllFS
    (l,r,s,split->GetSplit(),split->GetSplit()->GetFlavour(),
     s->GetType()==pst::IS?(c->GetType()==pst::IS?3:2):
     (c->GetType()==pst::IS?1:0));
  msg_Debugging()<<"}\n";
  return nres;
}

bool Shower::UpdateDaughters(Parton *const split,Parton *const newpB,
			     Parton *const newpC,const bool del)
{
  newpB->SetStart(split->KtTest());
  newpC->SetStart(split->KtTest());
  newpB->SetKtMax(split->KtMax());
  newpC->SetKtMax(split->KtMax());
  newpB->SetVeto(split->KtVeto());
  newpC->SetVeto(split->KtVeto());
  newpB->SetKtPrev(split->KtPrev());
  newpC->SetKtPrev(split->KtPrev());
  newpB->SetNext(split->GetNext());
  newpB->SetKtNext(split->KtNext());
  newpB->SetStat(split->Stat());
  if (split->GetNext()) {
    split->GetNext()->SetPrev(newpB);
    split->SetNext(newpB);
    newpB->SetPrev(split);
  }
  newpB->UpdateDaughters();
  split->GetSpect()->UpdateDaughters();
  if (!ReconstructDaughters(split->GetSing())) {
    if (split->GetNext()) {
      newpB->GetNext()->SetPrev(split);
      split->SetNext(newpB->GetNext());
    }
    if (del) {
      delete newpB;
      delete newpC;
    }
    return false;
  }
  newpB->SetPrev(split->GetPrev());
  if (split==split->GetSing()->GetSplit()) {
    split->GetSing()->SetSplit(newpB);
    split->GetSing()->GetLeft()->SetPrev(newpB);
    split->GetSing()->GetRight()->SetPrev(newpB);
  }
  return true;
}

bool Shower::EvolveSinglet(Singlet * act,const size_t &maxem,size_t &nem)
{
  p_actual        = act;
  bool mustshower = true;
  PLiter   splitter;
  Parton * newpB, * newpC, *newpA(NULL), *spect(NULL);
  Vec4D mom;
  double kt2win, kt2old(std::numeric_limits<double>::max());
  int mustsplit(0);
  
  if (nem>=maxem) return true;
  while (mustshower) {
    kt2win = 0.;
    splitter = SelectSplitting(kt2win,kt2old);
    kt2old=kt2win;
    //no shower anymore 
    if (splitter==p_actual->end()) return true;
    else {
      msg_Debugging()<<"Emission "<<m_flavB<<" -> "<<m_flavB<<" "<<m_flavC
		     <<" at kt = "<<sqrt((*splitter)->KtTest())
		     <<"( "<<sqrt((*splitter)->KtNext())<<" .. "
		     <<sqrt((*splitter)->KtPrev())<<" ), z = "<<(*splitter)->ZTest()<<", y = "
		     <<(*splitter)->YTest()<<" for\n"<<**splitter
		     <<*(*splitter)->GetSpect()<<"\n";
      if ((*splitter)->GetSing()->GetLeft()) {
	if ((*splitter)->GetType()==pst::IS) {
	  if (m_flavA!=(*splitter)->GetFlavour()) continue;
	}
	else {
	  if (m_flavB!=(*splitter)->GetFlavour()) continue;
	}
      }
      if (kt2win<(*splitter)->KtNext()) {
	msg_Debugging()<<"... Defer split ...\n\n";
	return true;
      }
      if (kt2win>(*splitter)->KtPrev()) {
	(*splitter)->SetStart(kt2win);
	msg_Debugging()<<"... Veto split ...\n\n";
	continue;
      }
      Vec4D splitorig((*splitter)->Momentum()), spectorig((*splitter)->GetSpect()->Momentum());
      //the FF case 
      if ((*splitter)->GetType()==pst::FS && (*splitter)->GetSpect()->GetType()==pst::FS) {
	newpC=NULL;
	if ((*splitter)->KtTest()<=(*splitter)->KtMax()) m_kinFF.SetJF(NULL);
	else m_kinFF.SetJF((*splitter)->GetSing()->JF());
	msg_Debugging()<<sqrt((*splitter)->KtTest())<<" vs. "<<sqrt((*splitter)->KtMax())
		       <<" -> "<<(*splitter)->GetSing()->JF()<<" vs. "<<m_kinFF.JF()<<"\n";
	int stat(m_kinFF.MakeKinematics((*splitter),m_flavB,m_flavC,newpC));
	if (stat==-1) {
	  (*splitter)->SetMomentum(splitorig);
	  (*splitter)->GetSpect()->SetMomentum(spectorig);
	  continue;
	}
	if (stat==0) return false;
	mom       = (*splitter)->Momentum();
	newpB     = new Parton(m_flavB,mom,(*splitter)->GetType());
	newpB->SetId((*splitter)->Id());
	spect     = (*splitter)->GetSpect();
	p_actual->push_back(newpC);
	bool ustat(UpdateDaughters(*splitter,newpB,newpC));
	p_actual->erase(--p_actual->end());
	if (!ustat) {
	  (*splitter)->SetMomentum(splitorig);
	  spect->SetMomentum(spectorig);
	  msg_Debugging()<<"Save history for\n"<<**splitter
		     <<*(*splitter)->GetSpect()<<"\n";
	  (*splitter)->UpdateDaughters();
	  spect->UpdateDaughters();
	  (*splitter)->SetStart(kt2win);
	  if (!ReconstructDaughters((*splitter)->GetSing())) {
	    msg_Error()<<METHOD<<"(): Invalid history."<<std::endl;
	    return false;
	  }
	  continue;
	}
	mustsplit = p_actual->SplitParton(splitter,newpB,newpC);
	m_last[0]=newpB;
	m_last[1]=newpC;
	m_last[2]=spect;
      }
      //the FI case 
      else if ((*splitter)->GetType()==pst::FS && (*splitter)->GetSpect()->GetType()==pst::IS) {
	newpC=NULL;
	if ((*splitter)->KtTest()<=(*splitter)->KtMax()) m_kinFI.SetJF(NULL);
	else m_kinFI.SetJF((*splitter)->GetSing()->JF());
	msg_Debugging()<<sqrt((*splitter)->KtTest())<<" vs. "<<sqrt((*splitter)->KtMax())
		       <<" -> "<<(*splitter)->GetSing()->JF()<<" vs. "<<m_kinFI.JF()<<"\n";
	int stat(m_kinFI.MakeKinematics((*splitter),m_flavB,m_flavC,newpC));
	if (stat>0) {
	  stat=RemnantTest((*splitter)->GetSpect());
	  if (stat==-1) delete newpC;
	}
	if (stat==-1) {
	  (*splitter)->SetMomentum(splitorig);
	  (*splitter)->GetSpect()->SetMomentum(spectorig);
	  continue;
	}
	if (stat==0) return false;
	mom       = (*splitter)->Momentum();
	newpB     = new Parton(m_flavB,mom,(*splitter)->GetType());
	newpB->SetId((*splitter)->Id());
	newpB->SetSing((*splitter)->GetSing());
	spect     = (*splitter)->GetSpect();
	// Boost the full thing into the c.m. frame
 	p_actual->BoostAllFS(newpB,newpC,spect,*splitter,
			     (*splitter)->GetFlavour(),2);
	p_actual->push_back(newpC);
	bool ustat(UpdateDaughters(*splitter,newpB,newpC,false));
	p_actual->erase(--p_actual->end());
	if (!ustat) {
	  p_actual->BoostBackAllFS(newpB,newpC,spect,*splitter,
				   (*splitter)->GetFlavour(),2);
	  delete newpB;
	  delete newpC;
	  (*splitter)->SetMomentum(splitorig);
	  spect->SetMomentum(spectorig);
	  msg_Debugging()<<"Save history for\n"<<**splitter
		     <<*(*splitter)->GetSpect()<<"\n";
	  (*splitter)->UpdateDaughters();
	  spect->UpdateDaughters();
	  (*splitter)->SetStart(kt2win);
	  if (!ReconstructDaughters((*splitter)->GetSing())) {
	    msg_Error()<<METHOD<<"(): Invalid history."<<std::endl;
	    return false;
	  }
	  stat=-1;
	}
	if (stat>0) {
	  spect->SetXbj(spect->Xbj()/(1.0-(*splitter)->YTest()));
	  mustsplit = p_actual->SplitParton(splitter,newpB,newpC);
	  m_last[0]=newpB;
	  m_last[1]=newpC;
	  m_last[2]=spect;
	}
      }
      //the IF case
      else if ((*splitter)->GetType()==pst::IS && (*splitter)->GetSpect()->GetType()==pst::FS) {
	newpC=NULL;
	if ((*splitter)->KtTest()<=(*splitter)->KtMax()) m_kinIF.SetJF(NULL);
	else m_kinIF.SetJF((*splitter)->GetSing()->JF());
	msg_Debugging()<<sqrt((*splitter)->KtTest())<<" vs. "<<sqrt((*splitter)->KtMax())
		       <<" -> "<<(*splitter)->GetSing()->JF()<<" vs. "<<m_kinIF.JF()<<"\n";
	int stat(m_kinIF.MakeKinematics((*splitter),m_flavA,m_flavC,newpC));
	if (stat>0) {
	  stat=RemnantTest(*splitter);
	  if (stat==-1) delete newpC;
	}
	if (stat==-1) {
	  (*splitter)->SetMomentum(splitorig);
	  (*splitter)->GetSpect()->SetMomentum(spectorig);
	  continue;
	}
	if (stat==0) return false;
	mom       = (*splitter)->Momentum();
	newpA     = new Parton(m_flavA,mom,(*splitter)->GetType());
	newpA->SetId((*splitter)->Id());
	newpA->SetXbj((*splitter)->Xbj()/(*splitter)->ZTest());
	newpA->SetBeam((*splitter)->Beam());
	newpA->SetSing((*splitter)->GetSing());
	spect     = (*splitter)->GetSpect();
	// Boost the full thing into the c.m. frame
 	p_actual->BoostAllFS(newpA,newpC,spect,*splitter,
			     (*splitter)->GetFlavour(),1);
	p_actual->push_back(newpC);
	bool ustat(UpdateDaughters(*splitter,newpA,newpC,false));
	p_actual->erase(--p_actual->end());
	if (!ustat) {
	  p_actual->BoostBackAllFS(newpA,newpC,spect,*splitter,
				   (*splitter)->GetFlavour(),1);
	  delete newpA;
	  delete newpC;
	  (*splitter)->SetMomentum(splitorig);
	  spect->SetMomentum(spectorig);
	  msg_Debugging()<<"Save history for\n"<<**splitter
		     <<*(*splitter)->GetSpect()<<"\n";
	  (*splitter)->UpdateDaughters();
	  spect->UpdateDaughters();
	  (*splitter)->SetStart(kt2win);
	  if (!ReconstructDaughters((*splitter)->GetSing())) {
	    msg_Error()<<METHOD<<"(): Invalid history."<<std::endl;
	    return false;
	  }
	  stat=-1;
	}
	if (stat>0) {
	  mustsplit = p_actual->SplitParton(splitter,newpA,newpC);
	  m_last[0]=newpA;
	  m_last[1]=newpC;
	  m_last[2]=spect;
	}
      }
      //the II case
      else if ((*splitter)->GetType()==pst::IS && (*splitter)->GetSpect()->GetType()==pst::IS) {
	splitorig=(*splitter)->Momentum();
	spectorig=(*splitter)->GetSpect()->Momentum();
	newpC=NULL;
	if ((*splitter)->KtTest()<=(*splitter)->KtMax()) m_kinII.SetJF(NULL);
	else m_kinII.SetJF((*splitter)->GetSing()->JF());
	msg_Debugging()<<sqrt((*splitter)->KtTest())<<" vs. "<<sqrt((*splitter)->KtMax())
		       <<" -> "<<(*splitter)->GetSing()->JF()<<" vs. "<<m_kinII.JF()<<"\n";
	int stat(m_kinII.MakeKinematics((*splitter),m_flavA,m_flavC,newpC));
	if (stat>0) {
	  stat=RemnantTest(*splitter);
	  if (stat==-1) delete newpC;
	}
	if (stat==-1) {
	  (*splitter)->SetMomentum(splitorig);
	  (*splitter)->GetSpect()->SetMomentum(spectorig);
	  continue;
	}
	if (stat==0) return false;
	mom       = (*splitter)->Momentum();
	newpA     = new Parton(m_flavA,mom,(*splitter)->GetType());
	newpA->SetId((*splitter)->Id());
	newpA->SetXbj((*splitter)->Xbj()/(*splitter)->ZTest());
	newpA->SetBeam((*splitter)->Beam());
	spect     = (*splitter)->GetSpect();
	// Boost the full thing into the c.m. frame
 	p_actual->BoostAllFS(newpA,newpC,spect,*splitter,
			     (*splitter)->GetFlavour(),3);
	p_actual->push_back(newpC);
	bool ustat(UpdateDaughters(*splitter,newpA,newpC,false));
	p_actual->erase(--p_actual->end());
	if (!ustat) {
	  p_actual->BoostBackAllFS(newpA,newpC,spect,*splitter,
				   (*splitter)->GetFlavour(),3);
	  delete newpA;
	  delete newpC;
	  (*splitter)->SetMomentum(splitorig);
	  spect->SetMomentum(spectorig);
	  msg_Debugging()<<"Save history for\n"<<**splitter
		     <<*(*splitter)->GetSpect()<<"\n";
	  (*splitter)->UpdateDaughters();
	  spect->UpdateDaughters();
	  (*splitter)->SetStart(kt2win);
	  if (!ReconstructDaughters((*splitter)->GetSing())) {
	    msg_Error()<<METHOD<<"(): Invalid history."<<std::endl;
	    return false;
	  }
	  stat=-1;
	}
	if (stat>0) {
	  mustsplit = p_actual->SplitParton(splitter,newpA,newpC);
	  m_last[0]=newpA;
	  m_last[1]=newpC;
	  m_last[2]=spect;
	}
      }
      else abort();
      msg_Debugging()<<"nem = "<<nem+1<<" vs. maxem = "<<maxem<<"\n";
      if (++nem>=maxem) return true;
    }
    //cout<<"-----------------------------------------------------------"<<endl<<(*p_actual);
  }
  return true;
}

PLiter Shower::SelectSplitting(double & kt2win,const double &kt2old) {
  PLiter winner = p_actual->end();
  for (PLiter splitter = p_actual->begin(); splitter!=p_actual->end();splitter++) {
    if (TrialEmission(kt2win,kt2old,*splitter)) winner = splitter;
  }
  return winner;
}

bool Shower::TrialEmission(double & kt2win,const double &kt2old,
			   Parton * split) 
{
  double kt2(0.),z(0.),y(0.),phi(0.);
  if (m_sudakov.Dice(split)) {
    m_sudakov.GetSplittingParameters(kt2,z,y,phi);
    if (kt2>kt2old) {
      split->SetStart(kt2);
      return false;
    }
    if (kt2>kt2win) {
      kt2win  = kt2;
      m_flavA = m_sudakov.GetFlavourA();
      m_flavB = m_sudakov.GetFlavourB();
      m_flavC = m_sudakov.GetFlavourC();
      split->SetTest(kt2,z,y,phi);
      return true;
    }
  }
  return false;
}

void Shower::SetMS(ATOOLS::Mass_Selector *const ms)
{
  m_sudakov.SetMS(ms);
  m_kinFF.SetMS(ms);
  m_kinFI.SetMS(ms);
  m_kinIF.SetMS(ms);
  m_kinII.SetMS(ms);
}

void Shower::SetRBMax(const double &rbmax)
{
  m_sudakov.SetRBMax(rbmax);
}
