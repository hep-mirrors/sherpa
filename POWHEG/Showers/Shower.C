#include "POWHEG/Showers/Shower.H"
#include "POWHEG/Tools/Parton.H"
#include "POWHEG/Main/CS_Gamma.H"
#include "PDF/Remnant/Remnant_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Phys/Cluster_Leg.H"
#include "ATOOLS/Org/My_Limits.H"

using namespace POWHEG;
using namespace ATOOLS;
using namespace std;

Shower::Shower(PDF::ISR_Handler * isr,const int qed,
	       Data_Reader *const dataread) : 
  p_actual(NULL), m_sudakov(isr,qed), p_isr(isr)
{
  int kfmode = dataread->GetValue<int>("PH_CSS_KFACTOR_SCHEME",1);
  double k0sq   = dataread->GetValue<double>("PH_CSS_PT2MIN",1);
  double is_fac = dataread->GetValue<double>("PH_CSS_AS_IS_FAC",1.0);
  double fs_fac = dataread->GetValue<double>("PH_CSS_AS_FS_FAC",1.0);
  m_kscheme = dataread->GetValue<int>("PH_CSS_KIN_SCHEME",1);
  std::vector<std::vector<std::string> > helpsvv;
  m_sudakov.SetShower(this);
  m_sudakov.InitSplittingFunctions(MODEL::s_model,kfmode);
  m_sudakov.SetCoupling(MODEL::s_model,k0sq,is_fac,fs_fac);
  m_kinFF.SetSudakov(&m_sudakov);
  m_kinFI.SetSudakov(&m_sudakov);
  m_kinIF.SetSudakov(&m_sudakov);
  m_kinII.SetSudakov(&m_sudakov);
  m_last[0]=m_last[1]=m_last[2]=m_last[3]=NULL;
  p_old[0]=Cluster_Leg::New(NULL,Vec4D(),kf_none,ColorID());
  p_old[1]=Cluster_Leg::New(NULL,Vec4D(),kf_none,ColorID());
}

Shower::~Shower()
{
  p_old[0]->Delete();
  p_old[1]->Delete();
}

bool Shower::EvolveShower(Singlet * actual,const size_t &maxem,size_t &nem)
{
  return EvolveSinglet(actual,maxem,nem);
}

int Shower::SetXBj(Parton *const p) const
{
  double x=0.0;
  if (p->Beam()==0) x=p->Momentum().PPlus()/rpa.gen.PBeam(0).PPlus();
  else x=p->Momentum().PMinus()/rpa.gen.PBeam(1).PMinus();
  if (x>1.0) return -1;
  p->SetXbj(x);
  return 1;
}

int Shower::RemnantTest(Parton *const p)
{
  if (p->Momentum()[0]<0.0 || p->Momentum().Nan()) return -1;
  if (p->Momentum()[0]>rpa.gen.PBeam(p->Beam())[0] &&
      !IsEqual(p->Momentum()[0],rpa.gen.PBeam(p->Beam())[0],1.0e-6)) return -1;
  return p_isr->GetRemnant(p->Beam())->
    TestExtract(p->GetFlavour(),p->Momentum())?1:-1;
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
  bool rd(true);
  if (rd) {
    if (newpB->GetType()==pst::IS &&
	RemnantTest(newpB)==-1) rd=false;
    if (split->GetSpect()->GetType()==pst::IS &&
	RemnantTest(split->GetSpect())==-1) rd=false;
  }
  int sci[2]={split->GetFlow(1),split->GetFlow(2)};
  Flavour sfi(split->GetFlavour());
  split->SetFlavour(newpB->GetFlavour());
  split->SetFlow(1,newpB->GetFlow(1));
  split->SetFlow(2,newpB->GetFlow(2));
  if (rd) rd=!p_gamma->Reject();
  split->SetFlavour(sfi);
  split->SetFlow(1,sci[0]);
  split->SetFlow(2,sci[1]);
  return rd;
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
  m_last[3]=split;
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
    kt2win = 0.;
    split = SelectSplitting(kt2win);
    //no shower anymore 
    if (split==NULL) {
      return true;
    }
    else {
      msg_Debugging()<<"Emission "<<m_flavA<<" -> "<<m_flavB<<" "<<m_flavC
		     <<" at kt = "<<sqrt(split->KtTest())
		     <<", z = "<<split->ZTest()<<", y = "
		     <<split->YTest()<<" for\n"<<*split
		     <<*split->GetSpect()<<"\n";
      m_last[0]=m_last[1]=m_last[2]=m_last[3]=NULL;
      if (kt2win>kt2old) {
	THROW(fatal_error,"Internal error");
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
	if (stat==0) return false;
	if (stat==-1) {
	  split->SetMomentum(splitorig);
	  split->GetSpect()->SetMomentum(spectorig);
	  ResetScales(split);
	  continue;
	}
	mom       = split->Momentum();
	newpB     = new Parton(m_flavB,mom,split->GetType());
	newpB->SetId(split->Id());
	newpB->SetKin(m_kscheme);
	spect     = split->GetSpect();
	newpC->SetKin(m_kscheme);
	SetSplitInfo(splitorig,spectorig,split,newpB,newpC,0);
	p_actual->push_back(newpC);
	bool ustat(UpdateDaughters(split,newpB,newpC));
	p_actual->erase(--p_actual->end());
	if (!ustat) {
	  delete newpB;
	  delete newpC;
	  split->SetMomentum(splitorig);
	  spect->SetMomentum(spectorig);
	  msg_Debugging()<<"Save history for\n"<<*split
		     <<*split->GetSpect()<<"\n";
	  ResetScales(split);
	  continue;
	}
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
	if (stat==0) return false;
	if (stat==-1) {
	  split->SetMomentum(splitorig);
	  split->GetSpect()->SetMomentum(spectorig);
	  ResetScales(split);
	  continue;
	}
	stat=1;
	mom       = split->Momentum();
	newpB     = new Parton(m_flavB,mom,split->GetType());
	newpB->SetId(split->Id());
	newpB->SetSing(split->GetSing());
	newpB->SetKin(m_kscheme);
	spect     = split->GetSpect();
	// Boost the full thing into the c.m. frame
	newpC->SetKin(m_kscheme);
	SetSplitInfo(splitorig,spectorig,split,newpB,newpC,2);
	p_actual->push_back(newpC);
 	p_actual->BoostAllFS(newpB,newpC,spect,split,
			     split->GetFlavour(),2);
	bool ustat(UpdateDaughters(split,newpB,newpC));
	p_actual->erase(--p_actual->end());
	if (!ustat) {
	  p_actual->BoostBackAllFS(newpB,newpC,spect,split,
				   split->GetFlavour(),2);
	  delete newpB;
	  delete newpC;
	  split->SetMomentum(splitorig);
	  spect->SetMomentum(spectorig);
	  msg_Debugging()<<"Save history for\n"<<*split
		     <<*split->GetSpect()<<"\n";
	  ResetScales(split);
	  continue;
	}
	SetXBj(spect);
	mustsplit = p_actual->SplitParton(split,newpB,newpC);
      }
      //the IF case
      else if (split->GetType()==pst::IS && split->GetSpect()->GetType()==pst::FS) {
	newpC=NULL;
	if (split->KtTest()<=split->KtMax()) m_kinIF.SetJF(NULL);
	else m_kinIF.SetJF(split->GetSing()->JF());
	msg_Debugging()<<sqrt(split->KtTest())<<" vs. "<<sqrt(split->KtMax())
		       <<" -> "<<split->GetSing()->JF()<<" vs. "<<m_kinIF.JF()<<"\n";
	int stat(m_kinIF.MakeKinematics(split,m_flavA,m_flavC,newpC));
	if (stat==0) return false;
	if (stat==-1) {
	  split->SetMomentum(splitorig);
	  split->GetSpect()->SetMomentum(spectorig);
	  ResetScales(split);
	  continue;
	}
	stat=1;
	mom       = split->Momentum();
	newpA     = new Parton(m_flavA,mom,split->GetType());
	newpA->SetId(split->Id());
	newpA->SetBeam(split->Beam());
	newpA->SetKin(m_kscheme);
	SetXBj(newpA);
	newpA->SetSing(split->GetSing());
	spect     = split->GetSpect();
	// Boost the full thing into the c.m. frame
	newpC->SetKin(m_kscheme);
	SetSplitInfo(splitorig,spectorig,split,newpA,newpC,1);
	p_actual->push_back(newpC);
 	p_actual->BoostAllFS(newpA,newpC,spect,split,
			     split->GetFlavour(),1);
	bool ustat(UpdateDaughters(split,newpA,newpC));
	p_actual->erase(--p_actual->end());
	if (!ustat) {
	  p_actual->BoostBackAllFS(newpA,newpC,spect,split,
				   split->GetFlavour(),1);
	  delete newpA;
	  delete newpC;
	  split->SetMomentum(splitorig);
	  spect->SetMomentum(spectorig);
	  msg_Debugging()<<"Save history for\n"<<*split
		     <<*split->GetSpect()<<"\n";
	  ResetScales(split);
	  continue;
	}
	mustsplit = p_actual->SplitParton(split,newpA,newpC);
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
	if (stat==0) return false;
	if (stat==-1) {
	  split->SetMomentum(splitorig);
	  split->GetSpect()->SetMomentum(spectorig);
	  ResetScales(split);
	  continue;
	}
	stat=1;
	mom       = split->Momentum();
	newpA     = new Parton(m_flavA,mom,split->GetType());
	newpA->SetId(split->Id());
	newpA->SetBeam(split->Beam());
	newpA->SetKin(m_kscheme);
	SetXBj(newpA);
	spect     = split->GetSpect();
	// Boost the full thing into the c.m. frame
	newpC->SetKin(m_kscheme);
	SetSplitInfo(splitorig,spectorig,split,newpA,newpC,3);
	p_actual->push_back(newpC);
 	p_actual->BoostAllFS(newpA,newpC,spect,split,
			     split->GetFlavour(),3);
	bool ustat(UpdateDaughters(split,newpA,newpC));
	p_actual->erase(--p_actual->end());
	if (!ustat) {
	  p_actual->BoostBackAllFS(newpA,newpC,spect,split,
				   split->GetFlavour(),3);
	  delete newpA;
	  delete newpC;
	  split->SetMomentum(splitorig);
	  spect->SetMomentum(spectorig);
	  msg_Debugging()<<"Save history for\n"<<*split
		     <<*split->GetSpect()<<"\n";
	  ResetScales(split);
	  continue;
	}
	mustsplit = p_actual->SplitParton(split,newpA,newpC);
      }
      else abort();
      msg_Debugging()<<"nem = "<<nem+1<<" vs. maxem = "<<maxem<<"\n";
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
  if (m_sudakov.Dice(split)) {
    m_sudakov.GetSplittingParameters(kt2,z,y,phi);
    split->SetSF(m_sudakov.Selected());
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

void Shower::SetRB(ATOOLS::RB_Map *const rbmap)
{
  m_sudakov.SetRB(rbmap);
}
