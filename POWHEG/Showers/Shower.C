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
  int kfmode = dataread->GetValue<int>("CSS_KFACTOR_SCHEME",1);
  double k0sqf = dataread->GetValue<double>("CSS_FS_PT2MIN",1.0);
  double k0sqi = dataread->GetValue<double>("CSS_IS_PT2MIN",4.0);
  double as_fac = dataread->GetValue<double>("CSS_AS_FAC",1.0);
  m_kscheme = dataread->GetValue<int>("PH_CSS_KIN_SCHEME",1);
  std::vector<std::vector<std::string> > helpsvv;
  m_sudakov.SetShower(this);
  m_sudakov.InitSplittingFunctions(MODEL::s_model,kfmode);
  m_sudakov.SetCoupling(MODEL::s_model,k0sqi,k0sqf,as_fac,as_fac);
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
  if (p->Beam()==0) x=p->Momentum().PPlus()/rpa->gen.PBeam(0).PPlus();
  else x=p->Momentum().PMinus()/rpa->gen.PBeam(1).PMinus();
  if (x>1.0) return -1;
  p->SetXbj(x);
  return 1;
}

int Shower::RemnantTest(Parton *const p)
{
  if (p->Momentum()[0]<0.0 || p->Momentum().Nan()) return -1;
  if (p->Momentum()[0]>rpa->gen.PBeam(p->Beam())[0] &&
      !IsEqual(p->Momentum()[0],rpa->gen.PBeam(p->Beam())[0],1.0e-6)) return -1;
  return p_isr->GetRemnant(p->Beam())->
    TestExtract(p->GetFlavour(),p->Momentum())?1:-1;
}

int Shower::UpdateDaughters(Parton *const split,Parton *const newpB,
			    Parton *const newpC,const int mode)
{
  newpB->SetStart(split->KtTest());
  newpC->SetStart(split->KtTest());
  newpB->SetKtMax(split->KtMax());
  newpC->SetKtMax(split->KtMax());
  newpB->SetVeto(split->KtVeto());
  newpC->SetVeto(split->KtVeto());
  int rd(1);
  if (rd) {
    if (newpB->GetType()==pst::IS &&
	RemnantTest(newpB)==-1) rd=-1;
    if (split->GetSpect()->GetType()==pst::IS &&
	RemnantTest(split->GetSpect())==-1) rd=-1;
  }
  int sci[2]={split->GetFlow(1),split->GetFlow(2)};
  m_flav=split->GetFlavour();
  split->SetFlavour(newpB->GetFlavour());
  split->SetFlow(1,newpB->GetFlow(1));
  split->SetFlow(2,newpB->GetFlow(2));
  if (rd==1) rd=p_gamma->Reject()?-1:1;
  if (rd==1 && split->KtTest()>split->KtMax())
    rd=!split->GetSing()->JetVeto(&m_sudakov);
  split->SetFlavour(m_flav);
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

int Shower::MakeKinematics
(Parton *split,const Flavour &fla,const Flavour &flb,
 const Flavour &flc,const int mode)
{
  DEBUG_FUNC(mode);
  Parton *spect(split->GetSpect()), *pj(NULL);
  Vec4D peo(split->Momentum()), pso(spect->Momentum());
  int stype(-1), stat(-1);
  if (split->GetType()==pst::FS) {
    if (spect->GetType()==pst::FS) {
      stype=0;
      stat=m_kinFF.MakeKinematics(split,flb,flc,pj);
    }
    else {
      stype=2;
      stat=m_kinFI.MakeKinematics(split,flb,flc,pj);
    }
  }
  else {
    if (spect->GetType()==pst::FS) {
      stype=1;
      stat=m_kinIF.MakeKinematics(split,fla,flc,pj);
    }
    else {
      stype=3;
      stat=m_kinII.MakeKinematics(split,fla,flc,pj);
    }
  }
  if (stat==-1) {
    split->SetMomentum(peo);
    spect->SetMomentum(pso);
    if (mode==0) ResetScales(split);
    delete pj;
    return stat;
  }
  Parton *pi(new Parton((stype&1)?fla:flb,
			split->Momentum(),split->GetType()));
  pi->SetSing(split->GetSing());
  pi->SetId(split->Id());
  pi->SetKin(m_kscheme);
  pj->SetKin(m_kscheme);
  pi->SetLT(split->LT());
  if (stype&1) pi->SetBeam(split->Beam());
  if (mode==0) SetSplitInfo(peo,pso,split,pi,pj,stype);
  split->GetSing()->push_back(pj);
  if (stype) split->GetSing()->BoostAllFS
    (pi,pj,spect,split,split->GetFlavour(),stype);
  Flavour fls(split->GetFlavour());
  if (mode!=0) split->SetFlavour(pi->GetFlavour());
  int ustat(UpdateDaughters(split,pi,pj,mode));
  split->GetSing()->pop_back();
  if (ustat<=0 || mode!=0) {
    split->SetFlavour(fls);
    if (stype) split->GetSing()->BoostBackAllFS
      (pi,pj,spect,split,split->GetFlavour(),stype);
    delete pi;
    delete pj;
    msg_Debugging()<<"Save history for\n"<<*split<<*spect<<"\n";
    split->SetMomentum(peo);
    spect->SetMomentum(pso);
    if (mode==0) ResetScales(split);
    return ustat;
  }
  split->GetSing()->SplitParton(split,pi,pj);
  return 1;
}

bool Shower::EvolveSinglet(Singlet * act,const size_t &maxem,size_t &nem)
{
  p_actual        = act;
  Parton * split;
  Vec4D mom;
  double kt2win, kt2old(std::numeric_limits<double>::max());
  
  if (nem>=maxem) return true;
  while (true) {
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
      int kstat(MakeKinematics(split,m_flavA,m_flavB,m_flavC,0));
      if (kstat<0) continue;
      if (kstat==0) return false;
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
  if (m_sudakov.Generate(split)) {
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
