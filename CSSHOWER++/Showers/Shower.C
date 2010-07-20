#include "CSSHOWER++/Showers/Shower.H"
#include "CSSHOWER++/Tools/Parton.H"
#include "PDF/Remnant/Remnant_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Phys/Cluster_Leg.H"

using namespace CSSHOWER;
using namespace ATOOLS;
using namespace std;

Shower::Shower(PDF::ISR_Handler * isr,const int qed,
	       Data_Reader *const dataread) : 
  p_actual(NULL), m_sudakov(isr,qed), p_isr(isr)
{
  int kfmode = dataread->GetValue<int>("CSS_KFACTOR_SCHEME",1);
  double k0sq   = dataread->GetValue<double>("CSS_PT2MIN",1);
  double is_fac = dataread->GetValue<double>("CSS_AS_IS_FAC",1.0);
  double fs_fac = dataread->GetValue<double>("CSS_AS_FS_FAC",1.0);
  m_kscheme = dataread->GetValue<int>("CSS_KIN_SCHEME",0);
  std::vector<std::vector<std::string> > helpsvv;
  dataread->MatrixFromFile(helpsvv,"CSS_ENHANCE");
  m_efac.clear();
  for (size_t i(0);i<helpsvv.size();++i)
    if (helpsvv[i].size()==2) {
      m_efac[helpsvv[i][0]]=ToType<double>(helpsvv[i][1]);
    }
  m_sudakov.SetShower(this);
  m_sudakov.InitSplittingFunctions(MODEL::s_model,kfmode);
  m_sudakov.SetCoupling(MODEL::s_model,k0sq,is_fac,fs_fac);
  m_kinFF.SetSudakov(&m_sudakov);
  m_kinFI.SetSudakov(&m_sudakov);
  m_kinIF.SetSudakov(&m_sudakov);
  m_kinII.SetSudakov(&m_sudakov);
  m_last[0]=NULL;
  m_last[1]=NULL;
  m_last[2]=NULL;
  p_old[0]=Cluster_Leg::New(NULL,Vec4D(),kf_none,ColorID());
  p_old[1]=Cluster_Leg::New(NULL,Vec4D(),kf_none,ColorID());
}

Shower::~Shower()
{
  p_old[0]->Delete();
  p_old[1]->Delete();
}

double Shower::EFac(const std::string &sfk) const 
{ 
  for (std::map<std::string,double,ATOOLS::String_Sort>::const_reverse_iterator
	 eit=m_efac.rbegin();eit!=m_efac.rend();++eit)
    if (sfk.find(eit->first)!=std::string::npos) return eit->second;
  return 1.0;
}

bool Shower::EvolveShower(Singlet * actual,const size_t &maxem,size_t &nem)
{
  m_weight=1.0;
  return EvolveSinglet(actual,maxem,nem);
}

double Shower::GetXBj(Parton *const p) const
{
  if (p->Beam()==0) return p->Momentum().PPlus()/rpa.gen.PBeam(0).PPlus();
  return p->Momentum().PMinus()/rpa.gen.PBeam(1).PMinus();
}

int Shower::SetXBj(Parton *const p) const
{
  double x(GetXBj(p));
  if (x>1.0) return -1;
  p->SetXbj(x);
  return 1;
}

int Shower::RemnantTest(Parton *const p)
{
  if (p->Momentum()[0]<0.0 || p->Momentum().Nan()) return -1;
  if (p->Momentum()[0]>rpa.gen.PBeam(p->Beam())[0] &&
      !IsEqual(p->Momentum()[0],rpa.gen.PBeam(p->Beam())[0],1.0e-6)) return -1;
  if (!m_sudakov.CheckPDF(GetXBj(p),p->GetFlavour(),p->Beam())) return -1;
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
  if (c->GetFlavour().Strong() &&
      l->GetFlavour().Strong() && r->GetFlavour().Strong()) {
    Parton *lr[2]={l,r};
    for (int n(1);n<=2;++n) {
      int sf(s->GetFlow(3-n)), smf(s->GetMEFlow(3-n));
      for (int m(0);m<2;++m) {
	int cf(lr[m]->GetFlow(n)), cmf(lr[m]->GetMEFlow(n));
	if (cmf && cmf==smf && (cf!=cmf || sf!=smf)) {
	  for (Singlet::const_iterator sit(l->GetSing()->begin());
	       sit!=l->GetSing()->end();++sit)
	    if ((*sit)->Id()!=s->Id() &&
		(*sit)->GetFlow(3-n)==cf) {
	      s->SetMomentum(s->GetPrev()->Momentum());
	      s->UpdateDaughters();
	      s=*sit;
	      msg_Debugging()<<"new spec for "<<cf<<" "<<*s<<"\n";
	      break;
	    }
	  break;
	}
      }
    }
  }
  Flavour fli(l->GetFlavour());
  int kin(l->Kin());
  s->SetMomentum(s->GetPrev()->Momentum());
  l->SetMomentum(c->Momentum());
  l->SetFlavour(c->GetFlavour());
  l->SetKin(c->Kin());
  l->SetSpect(s);
  msg_Debugging()<<"before: c: "<<*l<<"        s: "<<*s<<"\n";
  msg_Debugging()<<"kt = "<<sqrt(l->KtTest())<<", z = "
		 <<l->ZTest()<<", y = "<<l->YTest()
		 <<", phi = "<<l->Phi()<<", scheme = "<<l->Kin()<<"\n\n";
  int stat=0;
  if (c->GetType()==pst::FS) {
    if (s->GetPrev()->GetType()==pst::FS) {
      m_kinFF.SetJF(NULL);
      stat=m_kinFF.MakeKinematics(l,fli,r->GetFlavour(),r);
      l->SetFlavour(fli);
      l->SetKin(kin);
    }
    else {
      m_kinFI.SetJF(NULL);
      stat=m_kinFI.MakeKinematics(l,fli,r->GetFlavour(),r);
      l->SetFlavour(fli);
      l->SetKin(kin);
      if (stat>0) {
	split->BoostAllFS(l,r,s,split->GetSplit(),
			  split->GetSplit()->GetFlavour(),2);
	stat=RemnantTest(s);
	if (stat<=0) split->BoostBackAllFS
	  (l,r,s,split->GetSplit(),split->GetSplit()->GetFlavour(),2);
      }
    }
  }
  else {
    if (s->GetPrev()->GetType()==pst::FS) {
      m_kinIF.SetJF(NULL);
      stat=m_kinIF.MakeKinematics(l,fli,r->GetFlavour(),r);
      l->SetFlavour(fli);
      l->SetKin(kin);
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
      l->SetFlavour(fli);
      l->SetKin(kin);
      if (stat>0) {
	split->BoostAllFS(l,r,s,split->GetSplit(),
			  split->GetSplit()->GetFlavour(),3);
	stat=RemnantTest(l);
	if (stat>0) stat=RemnantTest(s);
	if (stat<=0) split->BoostBackAllFS
	  (l,r,s,split->GetSplit(),split->GetSplit()->GetFlavour(),3);
      }
    }
  }
  msg_Debugging()<<"after: l: "<<*l<<"       r: "<<*r<<"       s: "<<*s<<"\n";
  if (stat<0) {
    if (s!=split->GetSpec()) s->GetPrev()->UpdateDaughters();
    return false;
  }
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
  msg_Debugging()<<"} -> "<<nres<<"\n";
  if (s!=split->GetSpec()) s->GetPrev()->UpdateDaughters();
  return nres;
}

bool Shower::UpdateDaughters(Parton *const split,Parton *const newpB,
			     Parton *const newpC)
{
  newpB->SetStart(split->KtTest());
  newpC->SetStart(split->KtTest());
  newpB->SetKtMax(split->KtMax());
  newpC->SetKtMax(split->KtMax());
  newpB->SetVeto(split->KtVeto());
  newpC->SetVeto(split->KtVeto());
  newpB->SetKtPrev(split->KtPrev());
  newpC->SetKtPrev(split->KtPrev());
  newpB->SetKtNext(split->KtNext());
  newpB->SetStat(split->Stat());
  if (split->GetNext()) {
    split->GetNext()->SetPrev(newpB);
    newpB->SetNext(split->GetNext());
  }
  newpB->UpdateDaughters();
  newpC->UpdateNewDaughters();
  split->GetSpect()->UpdateDaughters();
  p_actual->ArrangeColours(split,newpB,newpC);
  bool rd(ReconstructDaughters(split->GetSing()));
  if (rd) {
    if (newpB->GetType()==pst::IS &&
	RemnantTest(newpB)==-1) rd=false;
    if (split->GetSpect()->GetType()==pst::IS &&
	RemnantTest(split->GetSpect())==-1) rd=false;
  }
  newpC->UpdateDaughters();
  p_actual->RemoveParton(newpC);
  if (!rd) {
    p_actual->RearrangeColours(split,newpB,newpC);
    if (split->GetNext()) {
      newpB->GetNext()->SetPrev(split);
      split->SetNext(newpB->GetNext());
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

void Shower::ResetScales(Parton *const split)
{
  for (PLiter pit(p_actual->begin());pit!=p_actual->end();++pit)
    (*pit)->SetStart(split->KtTest());
  m_last[0]=m_last[1]=m_last[2]=NULL;
}

void Shower::SetSplitInfo
(const Vec4D &psplit,const Vec4D &pspect,Parton *const split,
 Parton *const newb,Parton *const newc,const int mode)
{
  p_old[0]->SetMom((mode&1)?-psplit:psplit);
  p_old[1]->SetMom((mode&2)?-pspect:pspect);
  p_old[0]->SetFlav(split->GetFlavour());
  p_old[0]->SetCol(ColorID(split->GetFlow((mode&1)?2:1),
			   split->GetFlow((mode&1)?1:2)));
  m_last[0]=newb;
  m_last[1]=newc;
  m_last[2]=split->GetSpect();
}

bool Shower::EvolveSinglet(Singlet * act,const size_t &maxem,size_t &nem)
{
  p_actual        = act;
  bool mustshower = true;
  Parton * split, * newpB, * newpC, *newpA(NULL), *spect(NULL);
  Vec4D mom;
  double kt2win, kt2old(std::numeric_limits<double>::max());
  int mustsplit(0);
  
  if (nem>=maxem) return true;
  while (mustshower) {
    for (Singlet::const_iterator it=p_actual->begin();it!=p_actual->end();++it)
      if ((*it)->GetType()==pst::IS) SetXBj(*it);
    kt2win = 0.;
    split = SelectSplitting(kt2win);
    //no shower anymore 
    if (split==NULL) {
      for (Singlet::const_iterator it=p_actual->begin(); it!=p_actual->end();
           ++it) {
        if ((*it)->Weight()!=1.0)
          msg_Debugging()<<"Add wt for "<<(**it)<<": "<<(*it)->Weight()<<"\n";
        m_weight*=(*it)->Weight();
      }
      return true;
    }
    else {
      msg_Debugging()<<"Emission "<<m_flavA<<" -> "<<m_flavB<<" "<<m_flavC
		     <<" at kt = "<<sqrt(split->KtTest())
		     <<"( "<<sqrt(split->KtNext())<<" .. "
		     <<sqrt(split->KtPrev())<<" ), z = "<<split->ZTest()<<", y = "
		     <<split->YTest()<<" for\n"<<*split
		     <<*split->GetSpect()<<"\n";
      m_last[0]=m_last[1]=m_last[2]=NULL;
      if (kt2win<split->KtNext()) {
	msg_Debugging()<<"... Defer split ...\n\n";
	return true;
      }
      if (kt2win>Min(kt2old,split->KtPrev())) {
	THROW(fatal_error,"Internal error");
      }
      if (split->GetSing()->GetLeft()) {
	if (split->GetType()==pst::IS) {
	  if (m_flavA!=split->GetFlavour()) {
	    msg_Debugging()<<"... Veto flavour change ...\n\n";
	    ResetScales(split);
	    continue;
	  }
	}
	else {
	  if (m_flavB!=split->GetFlavour()) {
	    msg_Debugging()<<"... Veto flavour change ...\n\n";
	    ResetScales(split);
	    continue;
	  }
	}
      }
      Vec4D splitorig(split->Momentum()), spectorig(split->GetSpect()->Momentum());
      //the FF case 
      if (split->GetType()==pst::FS && split->GetSpect()->GetType()==pst::FS) {
	newpC=NULL;
	if (split->KtTest()<=split->KtMax()) m_kinFF.SetJF(NULL);
	else m_kinFF.SetJF(split->GetSing()->JF());
	msg_Debugging()<<sqrt(split->KtTest())<<" vs. "<<sqrt(split->KtMax())
		       <<" -> "<<split->GetSing()->JF()<<" vs. "<<m_kinFF.JF()<<"\n";
	int stat(m_kinFF.MakeKinematics(split,m_flavB,m_flavC,newpC));
	if (stat==-1) {
	  split->SetMomentum(splitorig);
	  split->GetSpect()->SetMomentum(spectorig);
	  ResetScales(split);
	  continue;
	}
	if (stat==0) return false;
	mom       = split->Momentum();
	newpB     = new Parton(m_flavB,mom,split->GetType());
	newpB->SetId(split->Id());
	newpB->SetKin(m_kscheme);
	spect     = split->GetSpect();
	newpC->SetKin(m_kscheme);
	SetSplitInfo(splitorig,spectorig,split,newpB,newpC,0);
	p_actual->AddParton(newpC);
	bool ustat(UpdateDaughters(split,newpB,newpC));
	if (!ustat) {
	  delete newpB;
	  newpC->DeleteAll();
	  split->SetMomentum(splitorig);
	  spect->SetMomentum(spectorig);
	  msg_Debugging()<<"Save history for\n"<<*split
		     <<*split->GetSpect()<<"\n";
	  split->UpdateDaughters();
	  spect->UpdateDaughters();
	  ResetScales(split);
	  if (!ReconstructDaughters(split->GetSing())) {
	    msg_Error()<<METHOD<<"(): Reconstruction error. Reject event."<<std::endl;
	    return false;
	  }
	  continue;
	}
	m_weight*=split->Weight();
	msg_Debugging()<<"sw = "<<split->Weight()
		       <<", w = "<<m_weight<<"\n";
	mustsplit = p_actual->SplitParton(split,newpB,newpC);
      }
      //the FI case 
      else if (split->GetType()==pst::FS && split->GetSpect()->GetType()==pst::IS) {
	newpC=NULL;
	if (split->KtTest()<=split->KtMax()) m_kinFI.SetJF(NULL);
	else m_kinFI.SetJF(split->GetSing()->JF());
	msg_Debugging()<<sqrt(split->KtTest())<<" vs. "<<sqrt(split->KtMax())
		       <<" -> "<<split->GetSing()->JF()<<" vs. "<<m_kinFI.JF()<<"\n";
	int stat(m_kinFI.MakeKinematics(split,m_flavB,m_flavC,newpC));
	if (stat==-1) {
	  split->SetMomentum(splitorig);
	  split->GetSpect()->SetMomentum(spectorig);
	  ResetScales(split);
	  continue;
	}
	if (stat==0) return false;
	mom       = split->Momentum();
	newpB     = new Parton(m_flavB,mom,split->GetType());
	newpB->SetId(split->Id());
	newpB->SetSing(split->GetSing());
	newpB->SetKin(m_kscheme);
	spect     = split->GetSpect();
	// Boost the full thing into the c.m. frame
	newpC->SetKin(m_kscheme);
	SetSplitInfo(splitorig,spectorig,split,newpB,newpC,2);
	p_actual->AddParton(newpC);
 	p_actual->BoostAllFS(newpB,newpC,spect,split,
			     split->GetFlavour(),2);
	bool ustat(UpdateDaughters(split,newpB,newpC));
	if (!ustat) {
	  p_actual->BoostBackAllFS(newpB,newpC,spect,split,
				   split->GetFlavour(),2);
	  delete newpB;
	  newpC->DeleteAll();
	  split->SetMomentum(splitorig);
	  spect->SetMomentum(spectorig);
	  msg_Debugging()<<"Save history for\n"<<*split
		     <<*split->GetSpect()<<"\n";
	  split->UpdateDaughters();
	  spect->UpdateDaughters();
	  ResetScales(split);
	  if (!ReconstructDaughters(split->GetSing())) {
	    msg_Error()<<METHOD<<"(): Reconstruction error. Reject event."<<std::endl;
	    return false;
	  }
	  stat=-1;
	}
	if (stat>0) {
	  m_weight*=split->Weight();
	  msg_Debugging()<<"sw = "<<split->Weight()
			 <<", w = "<<m_weight<<"\n";
	  mustsplit = p_actual->SplitParton(split,newpB,newpC);
	}
      }
      //the IF case
      else if (split->GetType()==pst::IS && split->GetSpect()->GetType()==pst::FS) {
	newpC=NULL;
	if (split->KtTest()<=split->KtMax()) m_kinIF.SetJF(NULL);
	else m_kinIF.SetJF(split->GetSing()->JF());
	msg_Debugging()<<sqrt(split->KtTest())<<" vs. "<<sqrt(split->KtMax())
		       <<" -> "<<split->GetSing()->JF()<<" vs. "<<m_kinIF.JF()<<"\n";
	int stat(m_kinIF.MakeKinematics(split,m_flavA,m_flavC,newpC));
	if (stat==-1) {
	  split->SetMomentum(splitorig);
	  split->GetSpect()->SetMomentum(spectorig);
	  ResetScales(split);
	  continue;
	}
	if (stat==0) return false;
	mom       = split->Momentum();
	newpA     = new Parton(m_flavA,mom,split->GetType());
	newpA->SetId(split->Id());
	newpA->SetBeam(split->Beam());
	newpA->SetKin(m_kscheme);
	newpA->SetSing(split->GetSing());
	spect     = split->GetSpect();
	// Boost the full thing into the c.m. frame
	newpC->SetKin(m_kscheme);
	SetSplitInfo(splitorig,spectorig,split,newpA,newpC,1);
	p_actual->AddParton(newpC);
 	p_actual->BoostAllFS(newpA,newpC,spect,split,
			     split->GetFlavour(),1);
	bool ustat(UpdateDaughters(split,newpA,newpC));
	if (!ustat) {
	  p_actual->BoostBackAllFS(newpA,newpC,spect,split,
				   split->GetFlavour(),1);
	  delete newpA;
	  newpC->DeleteAll();
	  split->SetMomentum(splitorig);
	  spect->SetMomentum(spectorig);
	  msg_Debugging()<<"Save history for\n"<<*split
		     <<*split->GetSpect()<<"\n";
	  split->UpdateDaughters();
	  spect->UpdateDaughters();
	  ResetScales(split);
	  if (!ReconstructDaughters(split->GetSing())) {
	    msg_Error()<<METHOD<<"(): Reconstruction error. Reject event."<<std::endl;
	    return false;
	  }
	  stat=-1;
	}
	if (stat>0) {
	  m_weight*=split->Weight();
	  msg_Debugging()<<"sw = "<<split->Weight()
			 <<", w = "<<m_weight<<"\n";
	  mustsplit = p_actual->SplitParton(split,newpA,newpC);
	}
      }
      //the II case
      else if (split->GetType()==pst::IS && split->GetSpect()->GetType()==pst::IS) {
	splitorig=split->Momentum();
	spectorig=split->GetSpect()->Momentum();
	newpC=NULL;
	if (split->KtTest()<=split->KtMax()) m_kinII.SetJF(NULL);
	else m_kinII.SetJF(split->GetSing()->JF());
	msg_Debugging()<<sqrt(split->KtTest())<<" vs. "<<sqrt(split->KtMax())
		       <<" -> "<<split->GetSing()->JF()<<" vs. "<<m_kinII.JF()<<"\n";
	int stat(m_kinII.MakeKinematics(split,m_flavA,m_flavC,newpC));
	if (stat==-1) {
	  split->SetMomentum(splitorig);
	  split->GetSpect()->SetMomentum(spectorig);
	  ResetScales(split);
	  continue;
	}
	if (stat==0) return false;
	mom       = split->Momentum();
	newpA     = new Parton(m_flavA,mom,split->GetType());
	newpA->SetId(split->Id());
	newpA->SetBeam(split->Beam());
	newpA->SetKin(m_kscheme);
	spect     = split->GetSpect();
	// Boost the full thing into the c.m. frame
	newpC->SetKin(m_kscheme);
	SetSplitInfo(splitorig,spectorig,split,newpA,newpC,3);
	p_actual->AddParton(newpC);
 	p_actual->BoostAllFS(newpA,newpC,spect,split,
			     split->GetFlavour(),3);
	bool ustat(UpdateDaughters(split,newpA,newpC));
	if (!ustat) {
	  p_actual->BoostBackAllFS(newpA,newpC,spect,split,
				   split->GetFlavour(),3);
	  delete newpA;
	  newpC->DeleteAll();
	  split->SetMomentum(splitorig);
	  spect->SetMomentum(spectorig);
	  msg_Debugging()<<"Save history for\n"<<*split
		     <<*split->GetSpect()<<"\n";
	  split->UpdateDaughters();
	  spect->UpdateDaughters();
	  ResetScales(split);
	  if (!ReconstructDaughters(split->GetSing())) {
	    msg_Error()<<METHOD<<"(): Reconstruction error. Reject event."<<std::endl;
	    return false;
	  }
	  stat=-1;
	}
	if (stat>0) {
	  m_weight*=split->Weight();
	  msg_Debugging()<<"sw = "<<split->Weight()
			 <<", w = "<<m_weight<<"\n";
	  mustsplit = p_actual->SplitParton(split,newpA,newpC);
	}
      }
      else abort();
      msg_Debugging()<<"nem = "<<nem+1<<" vs. maxem = "<<maxem<<"\n";
      if (m_last[0]) {
        for (Singlet::const_iterator it=p_actual->begin();
             it!=p_actual->end();++it) {
          if ((*it)->Weight()!=1.0) {
            msg_Debugging()<<"Add wt for "<<(**it)<<": "
                           <<(*it)->Weight(m_last[0]->KtStart())<<"\n";
            m_weight*=(*it)->Weight(m_last[0]->KtStart());
            (*it)->Weights().clear();
          }
          (*it)->SetKtPrev(m_last[0]->KtStart());
        }
      }
      else {
        for (Singlet::const_iterator it=p_actual->begin();
             it!=p_actual->end();++it) {
          (*it)->SetKtPrev(kt2win);
        }
      }
      if (++nem>=maxem) return true;
      kt2old=kt2win;
    }
    //cout<<"-----------------------------------------------------------"<<endl<<(*p_actual);
  }
  return true;
}

Parton *Shower::SelectSplitting(double & kt2win) {
  Parton *winner(NULL);
  for (PLiter splitter = p_actual->begin(); splitter!=p_actual->end();splitter++) {
    if (TrialEmission(kt2win,*splitter)) winner = *splitter;
  }
  return winner;
}

bool Shower::TrialEmission(double & kt2win,Parton * split) 
{
  double kt2(0.),z(0.),y(0.),phi(0.);
  while (true) {
  if (m_sudakov.Dice(split)) {
    m_sudakov.GetSplittingParameters(kt2,z,y,phi);
    split->SetWeight(m_sudakov.Weight());
    if (kt2>split->KtNext() && kt2>split->KtPrev()) {
      split->SetStart(kt2);
      continue;
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
  else {
    split->SetWeight(m_sudakov.Weight());
  }
  return false;
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
