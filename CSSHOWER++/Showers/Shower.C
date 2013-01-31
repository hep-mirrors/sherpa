#include "CSSHOWER++/Showers/Shower.H"
#include "CSSHOWER++/Tools/Parton.H"
#include "PDF/Remnant/Remnant_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Phys/Cluster_Leg.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/My_Limits.H"

using namespace CSSHOWER;
using namespace ATOOLS;
using namespace std;

Shower::Shower(PDF::ISR_Handler * isr,const int qed,const double & ktfac,
	       Data_Reader *const dataread) : 
  p_actual(NULL), m_sudakov(isr,qed), p_isr(isr)
{
  int kfmode  = ToType<int>(rpa->gen.Variable("CSS_KFACTOR_SCHEME"));
  m_k0sqf     = ToType<double>(rpa->gen.Variable("CSS_FS_PT2MIN"))*sqr(ktfac);
  m_k0sqi     = ToType<double>(rpa->gen.Variable("CSS_IS_PT2MIN"))*sqr(ktfac);
  m_fs_as_fac = ToType<double>(rpa->gen.Variable("CSS_FS_AS_FAC"));
  m_is_as_fac = ToType<double>(rpa->gen.Variable("CSS_IS_AS_FAC"));
  msg_Out()<<METHOD<<"(asfacs: IS = "<<m_is_as_fac<<", FS = "<<m_fs_as_fac<<")"
	   <<" from |"<<this<<"|\n"; 
  m_kscheme   = dataread->GetValue<int>("CSS_KIN_SCHEME",1);
  m_noem      = dataread->GetValue<int>("CSS_NOEM",0);
  std::vector<std::vector<std::string> > helpsvv;
  dataread->MatrixFromFile(helpsvv,"CSS_ENHANCE");
  m_efac.clear();
  for (size_t i(0);i<helpsvv.size();++i)
    if (helpsvv[i].size()==2) {
      m_efac[helpsvv[i][0]]=ToType<double>(helpsvv[i][1]);
    }
  m_sudakov.SetShower(this);
  m_sudakov.InitSplittingFunctions(MODEL::s_model,kfmode);
  m_sudakov.SetCoupling(MODEL::s_model,m_k0sqi,m_k0sqf,m_is_as_fac,m_fs_as_fac);

  msg_Out()<<METHOD<<"(asfacs: IS = "<<m_is_as_fac<<", FS = "<<m_fs_as_fac<<")"
	   <<"\n"; 
  m_last[0]=m_last[1]=m_last[2]=m_last[3]=NULL;
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
  if (p->Beam()==0) return p->Momentum().PPlus()/rpa->gen.PBeam(0).PPlus();
  return p->Momentum().PMinus()/rpa->gen.PBeam(1).PMinus();
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
  if (p->Momentum()[0]>rpa->gen.PBeam(p->Beam())[0] &&
      !IsEqual(p->Momentum()[0],rpa->gen.PBeam(p->Beam())[0],1.0e-6)) return -1;
  if (!m_sudakov.CheckPDF(GetXBj(p),p->GetFlavour(),p->Beam())) return -1;
  return p_isr->GetRemnant(p->Beam())->
    TestExtract(p->GetFlavour(),p->Momentum())?1:-1;
}

int Shower::ReconstructDaughters(Singlet *const split,const int mode,
				 Parton *const pi,Parton *const pj)
{
  if (split==NULL) return 1;
  if (mode&2) return !split->JetVeto(&m_sudakov);
  if (split->GetLeft()==NULL) return 1;
  if (split->GetRight()==NULL) THROW(fatal_error,"Invalid tree structure");
  DEBUG_FUNC(split);
  Parton *l(split->GetLeft()), *r(split->GetRight());
  Parton *c(split->GetSplit()->FollowUp()), *s(split->GetSpec());
  int kin(l->Kin()), ckin(c->Kin());
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
  double mi2(l->Mass2()), mj2(r->Mass2());
  Vec4D spi(l->FixSpec()), spj(r->FixSpec()), spk(s->FixSpec());
  Vec4D opi(l->OldMomentum()), opj(r->OldMomentum()), opk(s->OldMomentum());
  Flavour fli(l->GetFlavour()), flj(r->GetFlavour());
  s->SetMomentum(s->GetPrev()->Momentum());
  l->SetMomentum(c->Momentum());
  l->SetFlavour(c->GetFlavour());
  l->SetMass2(c->Mass2());
  l->SetKin(ckin);
  l->SetSpect(s);
  msg_Debugging()<<"before: c: "<<*l<<"        s: "<<*s<<"\n";
  msg_Debugging()<<"kt = "<<sqrt(l->KtTest())<<", z = "
		 <<l->ZTest()<<", y = "<<l->YTest()
		 <<", phi = "<<l->Phi()<<", scheme = "<<l->Kin()<<"\n\n";
  int stat=0;
  if (c->GetType()==pst::FS) {
    if (s->GetPrev()->GetType()==pst::FS) {
      stat=m_kinFF.MakeKinematics(l,mi2,mj2,flj,r,1);
      l->SetFlavour(fli);
      l->SetMass2(mi2);
    }
    else {
      stat=m_kinFI.MakeKinematics(l,mi2,mj2,flj,r,1);
      l->SetFlavour(fli);
      l->SetMass2(mi2);
      if (stat>0) stat=RemnantTest(s);
      if (stat>0)
	split->BoostAllFS(l,r,s,c,c->GetFlavour(),4|2);
    }
  }
  else {
    if (s->GetPrev()->GetType()==pst::FS) {
      stat=m_kinIF.MakeKinematics(l,mi2,mj2,flj,r,1);
      l->SetFlavour(fli);
      l->SetMass2(mi2);
      if (stat>0) stat=RemnantTest(l);
      if (stat>0)
	split->BoostAllFS(l,r,s,c,c->GetFlavour(),4|1);
    }
    else {
      stat=m_kinII.MakeKinematics(l,mi2,mj2,flj,r,1);
      l->SetFlavour(fli);
      l->SetMass2(mi2);
      if (stat>0) stat=RemnantTest(l);
      if (stat>0)
	split->BoostAllFS(l,r,s,c,c->GetFlavour(),4|3);
    }
  }
  if (stat<=0) {
    l->SetFixSpec(spi);
    r->SetFixSpec(spj);
    s->SetFixSpec(spk);
    l->SetOldMomentum(opi);
    r->SetOldMomentum(opj);
    s->SetOldMomentum(opk);
  }
  l->SetKin(kin);
  msg_Debugging()<<"after: l: "<<*l<<"       r: "<<*r<<"       s: "<<*s<<"\n";
  if (stat<0) {
    if (s!=split->GetSpec()) s->GetPrev()->UpdateDaughters();
    return -1;
  }
  l->GetSing()->UpdateDaughters();
  if(l->GetSing()!=r->GetSing()) r->GetSing()->UpdateDaughters();
  if (mode&1) return 1;
  int nres(ReconstructDaughters(l->GetSing(),mode,pi,pj));
  if (nres>0 && l->GetSing()!=r->GetSing())
    nres=ReconstructDaughters(r->GetSing(),mode,pi,pj);
  if (stat>0) {
    l->SetKin(ckin);
    split->BoostBackAllFS
      (l,r,s,c,c->GetFlavour(),
       (s->GetType()==pst::IS?(c->GetType()==pst::IS?3:2):
	(c->GetType()==pst::IS?1:0))|4);
    l->SetKin(kin);
    l->SetFixSpec(spi);
    r->SetFixSpec(spj);
    s->SetFixSpec(spk);
    l->SetOldMomentum(opi);
    r->SetOldMomentum(opj);
    s->SetOldMomentum(opk);
  }
  msg_Debugging()<<"} -> "<<nres<<"\n";
  if (s!=split->GetSpec()) s->GetPrev()->UpdateDaughters();
  return nres;
}

int Shower::UpdateDaughters(Parton *const split,Parton *const newpB,
			    Parton *const newpC,int mode)
{
  DEBUG_FUNC("");
  newpB->SetStart(split->KtTest());
  newpC->SetStart(split->KtTest());
  newpB->SetKtMax(split->KtMax());
  newpC->SetKtMax(split->KtMax());
  newpB->SetVeto(split->KtVeto());
  newpC->SetVeto(split->KtVeto());
  newpB->SetStat(split->Stat());
  if (split->GetNext()) {
    split->GetNext()->SetPrev(newpB);
    newpB->SetNext(split->GetNext());
  }
  newpB->SetId(split->Id());
  newpC->SetId(split->Id());
  newpB->UpdateDaughters();
  newpC->UpdateNewDaughters(newpB);
  split->GetSpect()->UpdateDaughters();
  if (mode==0) split->GetSing()->ArrangeColours(split,newpB,newpC);
  else {
    newpB->SetFlow(1,split->GetFlow(1));
    newpB->SetFlow(2,split->GetFlow(2));
  }
  split->SetFlow(1,newpB->GetFlow(1));
  split->SetFlow(2,newpB->GetFlow(2));
  newpB->SetPrev(split->GetPrev());
  newpC->SetPrev(split->GetPrev());
  newpB->SetFixSpec(split->FixSpec());
  newpB->SetOldMomentum(split->OldMomentum());
  int rd(ReconstructDaughters(split->GetSing(),mode,newpB,newpC));
  split->GetSing()->RemoveParton(newpC);
  if (rd<=0 || (mode&2)) {
    if (mode==0) split->GetSing()->RearrangeColours(split,newpB,newpC);
    if (split->GetNext()) {
      newpB->GetNext()->SetPrev(split);
      split->SetNext(newpB->GetNext());
    }
    return rd;
  }
  if (split==split->GetSing()->GetSplit()) {
    split->GetSing()->SetSplit(newpB);
    split->GetSing()->GetLeft()->SetPrev(newpB);
    split->GetSing()->GetRight()->SetPrev(newpB);
  }
  return rd;
}

void Shower::ResetScales(const double &kt2)
{
  for (PLiter pit(p_actual->begin());pit!=p_actual->end();++pit)
    if ((*pit)->GetPrev()==NULL ||
	(*pit)->GetPrev()->KScheme()!=1) (*pit)->SetStart(kt2);
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
  Vec4D pem(split->OldMomentum()), psm(spect->OldMomentum());
  Vec4D pef(split->FixSpec()), psf(spect->FixSpec());
  int stype(-1), stat(-1);
  double mc2(m_kinFF.MS()->Mass2(flc)), mi2(0.0);
  if (split->GetType()==pst::FS) {
    mi2=m_kinFF.MS()->Mass2(flb);
    if (mode==0 && split->KScheme()) mi2=split->Mass2();
    if (spect->GetType()==pst::FS) {
      stype=0;
      stat=m_kinFF.MakeKinematics(split,mi2,mc2,flc,pj);
    }
    else {
      stype=2;
      stat=m_kinFI.MakeKinematics(split,mi2,mc2,flc,pj);
    }
  }
  else {
    mi2=m_kinFF.MS()->Mass2(fla);
    if (spect->GetType()==pst::FS) {
      stype=1;
      stat=m_kinIF.MakeKinematics(split,mi2,mc2,flc,pj);
    }
    else {
      stype=3;
      stat=m_kinII.MakeKinematics(split,mi2,mc2,flc,pj);
    }
  }
  if (stat==1) {
    if (split->GetType()==pst::IS &&
	RemnantTest(split)==-1) stat=-1;
    if (split->GetSpect()->GetType()==pst::IS &&
	RemnantTest(split->GetSpect())==-1) stat=-1;
  }
  if (stat==-1) {
    split->SetMomentum(peo);
    spect->SetMomentum(pso);
    split->SetOldMomentum(pem);
    spect->SetOldMomentum(psm);
    split->SetFixSpec(pef);
    spect->SetFixSpec(psf);
    delete pj;
    return stat;
  }
  Parton *pi(new Parton((stype&1)?fla:flb,
			split->Momentum(),split->GetType()));
  pi->SetMass2(mi2);
  pi->SetSing(split->GetSing());
  pi->SetId(split->Id());
  pi->SetKScheme(split->KScheme());
  pi->SetKin(split->Kin());
  pj->SetKin(m_kscheme);
  pi->SetLT(split->LT());
  if (stype&1) pi->SetBeam(split->Beam());
  if (mode==0) SetSplitInfo(peo,pso,split,pi,pj,stype);
  split->GetSing()->AddParton(pj);
  if (stype) split->GetSing()->BoostAllFS
    (pi,pj,spect,split,split->GetFlavour(),stype);
  Flavour fls(split->GetFlavour());
  if (mode!=0) split->SetFlavour(pi->GetFlavour());
  int ustat(UpdateDaughters(split,pi,pj,mode));
  if (ustat<=0 || mode!=0) {
    split->SetFlavour(fls);
    if (stype) split->GetSing()->BoostBackAllFS
      (pi,pj,spect,split,split->GetFlavour(),stype);
    delete pi;
    pj->DeleteAll();
    split->SetMomentum(peo);
    spect->SetMomentum(pso);
    split->SetOldMomentum(pem);
    spect->SetOldMomentum(psm);
    split->SetFixSpec(pef);
    spect->SetFixSpec(psf);
    msg_Debugging()<<"Save history for\n"<<*split<<*spect<<"\n";
    split->UpdateDaughters();
    spect->UpdateDaughters();
    if (mode==0) {
      if (!ReconstructDaughters(split->GetSing(),0)) {
	msg_Error()<<METHOD<<"(): Reconstruction error. Reject event."<<std::endl;
	return 0;
      }
    }
    return ustat;
  }
  m_weight*=split->Weight();
  msg_Debugging()<<"sw = "<<split->Weight()
		 <<", w = "<<m_weight<<"\n";
  split->GetSing()->SplitParton(split,pi,pj);
  return 1;
}

bool Shower::EvolveSinglet(Singlet * act,const size_t &maxem,size_t &nem)
{
  p_actual=act;
  Vec4D mom;
  double kt2win;
  if (p_actual->NLO()&4) {
    msg_Debugging()<<"Skip MC@NLO emission\nSet p_T = "
		   <<sqrt(p_actual->KtNext())<<"\n";
    ResetScales(p_actual->KtNext());
    return true;
  }
  if (p_actual->GetSplit() &&
      p_actual->GetSplit()->KScheme()==1) {
    return true;
  }
  if (nem>=maxem) return true;
  while (true) {
    for (Singlet::const_iterator it=p_actual->begin();it!=p_actual->end();++it)
      if ((*it)->GetType()==pst::IS) SetXBj(*it);
    kt2win = 0.;
    Parton *split=SelectSplitting(kt2win);
    //no shower anymore 
    if (split==NULL) {
      ResetScales(p_actual->KtNext());
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
		     <<"( "<<sqrt(split->GetSing()->KtNext())<<" .. "
		     <<sqrt(split->KtStart())<<" ), z = "<<split->ZTest()<<", y = "
		     <<split->YTest()<<" for\n"<<*split
		     <<*split->GetSpect()<<"\n";
      m_last[0]=m_last[1]=m_last[2]=m_last[3]=NULL;
      if (kt2win<split->GetSing()->KtNext()) {
	msg_Debugging()<<"... Defer split ...\n\n";
	ResetScales(split->GetSing()->KtNext());
	return true;
      }
      ResetScales(kt2win);
      if (split->GetSing()->GetLeft()) {
	static double maxt(sqrt(std::numeric_limits<double>::max()));
	if (split->GetType()==pst::IS) {
	  if (m_flavA!=split->GetFlavour()) {
	    msg_Debugging()<<"... Veto flavour change ...\n\n";
	    continue;
	  }
	}
	else {
	  if (m_flavB!=split->GetFlavour()) {
	    msg_Debugging()<<"... Veto flavour change ...\n\n";
	    continue;
	  }
	}
      }
      if (p_actual->JF() &&
	  split->KtTest()<split->KtMax()) {
	msg_Debugging()<<"Highest Multi -> Disable jet veto\n";
	Singlet *sing(p_actual);
	sing->SetJF(NULL);
	while (sing->GetLeft()) {
	  sing=sing->GetLeft()->GetSing();
	  sing->SetJF(NULL);
	}
      }
      if (p_actual->JF()) {
	int vstat(MakeKinematics(split,m_flavA,m_flavB,m_flavC,2));
	if (vstat==0) {
	  if (p_actual->NLO()&2) {
	    msg_Debugging()<<"Skip first truncated emission, K = "
			   <<p_actual->LKF()<<"\n";
	    m_weight*=1.0/p_actual->LKF();
	    p_actual->SetNLO(0);
	    continue;
	  }
	  msg_Debugging()<<"Shower jet veto\n";
	  return false;
	}
      }
      if (m_noem) continue;
      int kstat(MakeKinematics(split,m_flavA,m_flavB,m_flavC,0));
      if (kstat<0) continue;
      if (p_actual->JF()) {
	if (p_actual->GetSplit()) {
	  msg_Debugging()<<"Truncated shower veto\n";
	  if (!(p_actual->NLO()&8)) return false;
	  else msg_Debugging()<<"Skip truncated shower veto\n";
	}
	msg_Debugging()<<"Disable jet veto\n";
	Singlet *sing(p_actual);
	sing->SetJF(NULL);
	while (sing->GetLeft()) {
	  sing=sing->GetLeft()->GetSing();
	  sing->SetJF(NULL);
	}
      }
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
        }
      }
      if (++nem>=maxem) return true;
    }
  }
  return true;
}

Parton *Shower::SelectSplitting(double & kt2win) {
  Parton *winner(NULL);
  for (PLiter splitter = p_actual->begin(); 
       splitter!=p_actual->end();splitter++) {
    if (TrialEmission(kt2win,*splitter)) winner = *splitter;
  }
  return winner;
}

bool Shower::TrialEmission(double & kt2win,Parton * split) 
{
  if (split->KtStart()==0.0 ||
      split->KtStart()<split->GetSing()->KtNext()) return false;
  double kt2(0.),z(0.),y(0.),phi(0.);
  while (true) {
  if (m_sudakov.Generate(split)) {
    m_sudakov.GetSplittingParameters(kt2,z,y,phi);
    split->SetWeight(m_sudakov.Weight());
    if (kt2>kt2win) {
      kt2win  = kt2;
      m_flavA = m_sudakov.GetFlavourA();
      m_flavB = m_sudakov.GetFlavourB();
      m_flavC = m_sudakov.GetFlavourC();
      split->SetCol(m_sudakov.GetCol());
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
