#include "CSSHOWER++/Showers/Shower.H"
#include "CSSHOWER++/Tools/Parton.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
#include "REMNANTS/Main/Remnant_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Phys/Cluster_Leg.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace CSSHOWER;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

Shower::Shower(PDF::ISR_Handler * isr,
               const int qcd, const int qed,
               int type) :
  p_actual(NULL), m_sudakov(isr,qcd,qed), p_isr(isr),
  p_variationweights(NULL)
{
  Settings& s = Settings::GetMainSettings();
  const int evol{ s["CSS_EVOLUTION_SCHEME"].Get<int>() };
  int kfmode{ s["CSS_KFACTOR_SCHEME"].Get<int>() };
  const int scs{ s["CSS_SCALE_SCHEME"].Get<int>() };
  double k0sqf{ s["CSS_FS_PT2MIN"].Get<double>() };
  double k0sqi{ s["CSS_IS_PT2MIN"].Get<double>() };
  double fs_as_fac{ s["CSS_FS_AS_FAC"].Get<double>() };
  double is_as_fac{ s["CSS_IS_AS_FAC"].Get<double>() };
  const double mth{ s["CSS_MASS_THRESHOLD"].Get<double>() };
  const bool reweightalphas = s["CSS_REWEIGHT_ALPHAS"].Get<int>();
  const bool reweightpdfs = s["CSS_REWEIGHT_PDFS"].Get<int>();
  m_maxrewem  = s["REWEIGHT_MAXEM"].Get<unsigned int>();
  m_use_bbw   = s["CSS_USE_BBW"].Get<int>();
  m_kscheme   = s["CSS_KIN_SCHEME"].Get<int>();
  m_recdec    = s["CSS_RECO_DECAYS"].Get<int>();
  m_maxpart   = s["CSS_MAXPART"].Get<int>();
  if (type) {
    kfmode=s["MI_CSS_KFACTOR_SCHEME"].Get<int>();
    k0sqf=s["MI_CSS_FS_PT2MIN"].Get<double>();
    k0sqi=s["MI_CSS_IS_PT2MIN"].Get<double>();
    fs_as_fac=s["MI_CSS_FS_AS_FAC"].Get<double>();
    is_as_fac=s["MI_CSS_IS_AS_FAC"].Get<double>();
    m_kscheme = s["MI_CSS_KIN_SCHEME"].Get<int>();
  }
  std::vector<std::vector<std::string> > helpsvv{
    s["CSS_ENHANCE"].GetMatrix<std::string>() };
  m_efac.clear();
  for (size_t i(0);i<helpsvv.size();++i)
    if (helpsvv[i].size()==2) {
      m_efac[helpsvv[i][0]]=ToType<double>(helpsvv[i][1]);
    }
  m_sudakov.SetShower(this);
  m_sudakov.SetMassThreshold(mth);
  m_sudakov.SetScaleScheme(scs);
  std::pair<double, double> pdfmin;
  pdfmin.first = s["CSS_PDF_MIN"].Get<double>();
  pdfmin.second = s["CSS_PDF_MIN_X"].Get<double>();
  m_sudakov.SetPDFMin(pdfmin);
  m_sudakov.InitSplittingFunctions(MODEL::s_model,kfmode);
  m_sudakov.SetCoupling(MODEL::s_model,k0sqi,k0sqf,is_as_fac,fs_as_fac);
  m_sudakov.SetReweightAlphaS(reweightalphas);
  m_sudakov.SetReweightPDFs(reweightpdfs);
  m_sudakov.SetReweightScaleCutoff(
      s["CSS_REWEIGHT_SCALE_CUTOFF"].Get<double>());
  m_kinFF.SetEvolScheme(evol);
  m_kinFI.SetEvolScheme(evol);
  m_kinIF.SetEvolScheme(evol);
  m_kinII.SetEvolScheme(evol);
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
  if (nem < m_maxrewem) {
    m_sudakov.SetVariationWeights(p_variationweights);
  } else {
    m_sudakov.SetVariationWeights(NULL);
  }
  return EvolveSinglet(actual,maxem,nem);
}

double Shower::GetXBj(Parton *const p) const
{
  return p_isr->CalcX(p->Momentum());
}

int Shower::SetXBj(Parton *const p) const
{
  double x(GetXBj(p));
  if (x>1.0) return -1;
  p->SetXbj(x);
  return 1;
}

int Shower::RemnantTest(Parton *const p,const Poincare_Sequence *lt)
{
  Vec4D mom(p->Momentum());
  if (lt) mom=(*lt)*mom;
  if (mom[0]<0.0 || mom.Nan()) return -1;
  double x(p_isr->CalcX(mom));
  if (x>1.0 && !IsEqual(x,1.0,1.0e-6)) return -1;
  if (!m_sudakov.CheckPDF(mom[0]/rpa->gen.PBeam(p->Beam())[0],p->GetFlavour(),p->Beam())) return -1;
  return p_remnants->GetRemnant(p->Beam())->TestExtract(p->GetFlavour(),mom)?1:-1;
}

int Shower::ReconstructDaughters(Singlet *const split,double &jcv,
				 Parton *const pi,Parton *const pj)
{
  if (split->GetSplit()) {
    if (split->GetSplit()->Stat()&2) {
      msg_Debugging()<<"Decay. Skip truncated shower veto\n";
    }
    else {
      msg_Debugging()<<"Truncated shower veto\n";
      return 0;
    }
  }
  jcv=split->JetVeto(&m_sudakov);
  return 1;
}

int Shower::UpdateDaughters(Parton *const split,Parton *const newpB,
			    Parton *const newpC,double &jcv,int mode)
{
  DEBUG_FUNC("");
  newpB->SetStart(split->KtTest());
  newpC->SetStart(split->KtTest());
  newpB->SetVeto(split->KtVeto());
  newpC->SetVeto(split->KtVeto());
  newpB->SetStat(split->Stat());
  newpB->SetFromDec(split->FromDec());
  newpC->SetFromDec(split->FromDec());
  if (split->GetNext()) {
    split->GetNext()->SetPrev(newpB);
    newpB->SetNext(split->GetNext());
  }
  newpB->SetId(split->Id());
  newpC->SetId(split->Id());
  split->GetSing()->ArrangeColours(split,newpB,newpC);
  newpB->SetPrev(split->GetPrev());
  newpC->SetPrev(split->GetPrev());
  newpB->SetFixSpec(split->FixSpec());
  newpB->SetOldMomentum(split->OldMomentum());
  double m2=split->Mass2();
  split->SetMass2(newpB->Mass2());
  Flavour fls(split->GetFlavour());
  split->SetFlavour(newpB->GetFlavour());
  int rd=ReconstructDaughters(split->GetSing(),jcv,newpB,newpC);
  split->SetFlavour(fls);
  if (mode && rd==1) rd=-1;
  split->SetMass2(m2);
  split->GetSing()->RemoveParton(newpC);
  if (rd==1 && (p_actual->NLO()&16)) rd=0;
  if (rd<=0) {
    split->GetSing()->RearrangeColours(split,newpB,newpC);
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
    if ((*pit)->KtStart()>kt2) (*pit)->SetStart(kt2);
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
 const Flavour &flc,double &jcv,int mode)
{
  DEBUG_FUNC("");
  Parton *spect(split->GetSpect()), *pj(NULL);
  Vec4D peo(split->Momentum()), pso(spect->Momentum());
  Vec4D pem(split->OldMomentum()), psm(spect->OldMomentum());
  Vec4D pef(split->FixSpec()), psf(spect->FixSpec());


  // first step for the amplitude: create parton list
  Singlet *s_test = split->GetSing();
  Parton_List p_list;
  for (Singlet::const_iterator it=s_test->begin();it!=s_test->end();++it){
      Parton *parton = *it;
      if (parton == split) continue;
      if (parton == spect) continue;
      p_list.push_back(parton);
  }
  s_test =NULL;

  int stype(-1), stat(-1);
  double mc2(m_kinFF.MS()->Mass2(flc)), mi2(0.0);
  if (split->GetType()==pst::FS) {
    mi2=m_kinFF.MS()->Mass2(flb);
    if (split->KScheme()) mi2=split->Mass2();
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
	RemnantTest(split,&split->LT())==-1) stat=-1;
    if (split->GetSpect()->GetType()==pst::IS &&
	RemnantTest(split->GetSpect(),
		    split->GetType()==pst::IS?
		    &split->LT():NULL)==-1) stat=-1;
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
			split->LT()*split->Momentum(),
			split->GetType()));
  pi->SetMass2(mi2);
  pi->SetSing(split->GetSing());
  pi->SetId(split->Id());
  pi->SetKScheme(split->KScheme());
  pi->SetKin(split->Kin());
  pj->SetKin(m_kscheme);
  pi->SetLT(split->LT());
  if (stype&1) pi->SetBeam(split->Beam());
  SetSplitInfo(peo,pso,split,pi,pj,stype);
  split->GetSing()->AddParton(pj);
  if (stype) split->GetSing()->BoostAllFS(pi,pj,spect);
  int ustat(UpdateDaughters(split,pi,pj,jcv,mode));
  if (ustat<=0 || split->GetSing()->GetLeft()) {
    if (stype) split->GetSing()->BoostBackAllFS(pi,pj,spect);
    delete pi;
    pj->DeleteAll();
    split->SetMomentum(peo);
    spect->SetMomentum(pso);
    split->SetOldMomentum(pem);
    spect->SetOldMomentum(psm);
    split->SetFixSpec(pef);
    spect->SetFixSpec(psf);
    return ustat;
  }
  m_weight*=split->Weight();
  if (p_variationweights) *p_variationweights*=split->Weight();
  msg_Debugging()<<"sw = "<<split->Weight()
		 <<", w = "<<m_weight<<"\n";
  split->GetSing()->SplitParton(split,pi,pj);

  // add missing particles to parton list
  p_list.push_back(pi);
  p_list.push_back(pj);
  p_list.push_back(spect);

  /*  p_list now contains all relevant partons. add first the IS partons and
      afterwards the FS partons to the amplitude*/

  Cluster_Amplitude * tmp_ampl = Cluster_Amplitude::New();
  size_t count_is(0);
  for(Parton_List::const_iterator it=p_list.begin(); it!=p_list.end(); ++it){
      Parton *parton = *it;
      if (parton->GetType()==pst::IS)  {
          PartonToAmplitude(parton, tmp_ampl);
          count_is++;
      }
    }

  for(Parton_List::const_iterator it=p_list.begin(); it!=p_list.end(); ++it){
      Parton *parton = *it;
      if (parton->GetType()==pst::FS)  PartonToAmplitude(parton, tmp_ampl);
    }

  tmp_ampl->SetNIn(count_is);
  //CheckAmplitude(tmp_ampl);
  p_actual->UpdateAmplitude(tmp_ampl);

  msg_Debugging() << "Amplitude after Make kinematics: " << *tmp_ampl << "\n";
  tmp_ampl=NULL;
  // end amplitude

  return 1;
}

double Shower::GetWeight(Variation_Parameters *params,
			 Variation_Weights *weights,
			 std::vector<double> &v)
{
  return v[weights->CurrentParametersIndex()];
}

double Shower::VetoWeight(Variation_Parameters *params,
			  Variation_Weights *weights,
			  JetVeto_Args &args)
{
  msg_Debugging()<<METHOD<<"("<<weights<<"){\n";
  int stat(args.m_jcv>=0.0);
  Singlet *sing(args.p_sing);
  Jet_Finder *jf(sing->JF());
  if (stat && jf) {
    double fac(weights?params->m_Qcutfac:1.0);
    if (jf) {
      stat=args.m_jcv<sqr(jf->Qcut()*fac);
      msg_Debugging()<<"  jcv = "<<sqrt(args.m_jcv)<<" vs "
		     <<jf->Qcut()<<" * "<<fac
		     <<" = "<<jf->Qcut()*fac<<"\n";
    }
  }
  if (stat==1) {
    msg_Debugging()<<"} no jet veto\n";
    args.m_acc=1;
    return 1.0;
  }
  if (sing->NLO()&2) {
    msg_Debugging()<<"  skip emission\n";
    if (weights==NULL) args.m_skip.back()=1;
    else args.m_skip[weights->CurrentParametersIndex()]=1;
    args.m_acc=1;
    msg_Debugging()<<"} no jet veto\n";
    return 1.0;
  }
  msg_Debugging()<<"} jet veto\n";
  return 0.0;
}

bool Shower::EvolveSinglet(Singlet * act,const size_t &maxem,size_t &nem)
{
  p_actual=act;

  Cluster_Amplitude * cl_tmp = Cluster_Amplitude::New();
  cl_tmp = SingletToAmplitude(act, cl_tmp);
  msg_Debugging() << "Amplitude from Singlet: " << *cl_tmp<< "\n" ;
  p_actual->UpdateAmplitude(cl_tmp);
  cl_tmp=NULL;


  Vec4D mom;
  double kt2win;
  if (p_actual->NLO()&128) {
    p_actual->Reduce();
    p_actual->SetNLO(0);
    ResetScales(p_actual->KtNext());
  }
  if (p_actual->NLO()&(4|8)) {
    msg_Debugging()<<"Skip MC@NLO emission\nSet p_T = "
		   <<sqrt(p_actual->KtNext())<<"\n";
    ResetScales(p_actual->KtNext());
    return true;
  }
  if (p_actual->GetSplit() &&
      (p_actual->GetSplit()->Stat()&4) &&
      !(p_actual->GetSplit()->Stat()&2)) {
    msg_Debugging()<<"Skip EW clustering\n";
    return true;
  }
  if (p_actual->NME()+nem>m_maxpart) {
    if (p_actual->NLO()&32) {
      p_actual->Reduce();
      p_actual->SetNLO(0);
    }
    return true;
  }
  if (nem>=maxem) {
    if (p_actual->NLO()&32) {
      p_actual->Reduce();
      p_actual->SetNLO(0);
    }
    return true;
  }
  while (true) {
    for (Singlet::const_iterator it=p_actual->begin();it!=p_actual->end();++it)
      if ((*it)->GetType()==pst::IS) SetXBj(*it);
    kt2win = 0.;
    Parton *split=SelectSplitting(kt2win);
    //no shower anymore 
    if (split==NULL) {
      msg_Debugging()<<"No emission\n";
      ResetScales(p_actual->KtNext());
      for (Singlet::const_iterator it=p_actual->begin(); it!=p_actual->end();
           ++it) {
        if ((*it)->Weight()!=1.0)
          msg_Debugging()<<"Add wt for "<<(**it)<<": "<<(*it)->Weight()<<"\n";
        m_weight*=(*it)->Weight();
	if (p_variationweights) *p_variationweights*=(*it)->Weight();
      }
      if (p_actual->NLO()&32) {
	p_actual->Reduce();
	p_actual->SetNLO(0);
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
	if (p_actual->NLO()&32) {
	  p_actual->Reduce();
	  p_actual->SetNLO(0);
	  ResetScales(p_actual->KtNext());
	}
	return true;
      }
      ResetScales(kt2win);
      if (p_actual->NSkip()) {
	msg_Debugging()<<"Skip emissions "<<p_actual->NSkip()<<"\n";
	p_actual->SetNSkip(p_actual->NSkip()-1);
	continue;
      }
      if (p_actual->JF() && p_actual->NMax() &&
	  (p_actual->GetSplit()==NULL ||
	   (p_actual->GetSplit()->Stat()&2))) {
	msg_Debugging()<<"Highest Multi -> Disable jet veto\n";
	Singlet *sing(p_actual);
	sing->SetJF(NULL);
	while (sing->GetLeft()) {
	  sing=sing->GetLeft()->GetSing();
	  sing->SetJF(NULL);
	}
      }
      double jcv(0.0);
      int kstat(MakeKinematics(split,m_flavA,m_flavB,m_flavC,jcv,
			       p_actual->NME()+nem>=m_maxpart));
      msg_Debugging()<<"stat = "<<kstat<<"\n";
      if (kstat<0) continue;
      if (p_actual->NLO()&64) {
	msg_Debugging()<<"UNLOPS veto\n";
	p_actual->Reduce();
	p_actual->SetNLO(0);
	ResetScales(p_actual->KtNext());
	continue;
      }
      JetVeto_Args vwa(p_actual,kstat?jcv:-1.0,
		       p_variationweights?
		       p_variationweights->NumberOfParameters()+1:1);
      m_weight*=VetoWeight(NULL,NULL,vwa);
      if (p_variationweights)
	p_variationweights->UpdateOrInitialiseWeights
	  (&Shower::VetoWeight,*this,vwa);
      if (p_actual->NLO()&2) {
	int nskip(0);
	for (size_t i(0);i<vwa.m_skip.size();++i)
	  nskip+=vwa.m_skip[i];
	double wskip(nskip/double(vwa.m_skip.size()));
	if (ran->Get()<=wskip) {
	  double lkf(p_actual->LKF());
	  Singlet *sing(p_actual);
	  sing->SetLKF(1.0);
	  sing->SetNLO(sing->NLO()&~2);
	  while (sing->GetLeft()) {
	    sing=sing->GetLeft()->GetSing();
	    sing->SetLKF(1.0);
	    sing->SetNLO(sing->NLO()&~2);
	  }
	  m_weight*=vwa.m_skip.back()?1.0/lkf/wskip:0.0;
	  std::vector<double> swa(vwa.m_skip.size()-1,0.0);
	  for (size_t i(0);i<swa.size();++i)
	    if (vwa.m_skip[i]) swa[i]=1.0/lkf/wskip;
	  msg_Debugging()<<"skip -> "<<m_weight<<" "<<swa<<"\n";
          if (p_variationweights)
            p_variationweights->UpdateOrInitialiseWeights(
                &Shower::GetWeight, *this, swa);
	  continue;
	}
	else {
	  m_weight*=vwa.m_skip.back()?0.0:1.0/(1.0-wskip);
	  std::vector<double> swa(vwa.m_skip.size()-1,0.0);
	  for (size_t i(0);i<swa.size();++i)
	    if (!vwa.m_skip[i]) swa[i]=1.0/(1.0-wskip);
	  msg_Debugging()<<"no skip -> "<<m_weight<<" "<<swa<<"\n";
          if (p_variationweights)
            p_variationweights->UpdateOrInitialiseWeights(
                &Shower::GetWeight, *this, swa);
	}
      }
      if (vwa.m_acc==0) return false;
      Singlet *sing(p_actual);
      sing->SetJF(NULL);
      while (sing->GetLeft()) {
	sing=sing->GetLeft()->GetSing();
	sing->SetJF(NULL);
      }
      if (m_last[0]) {
        for (Singlet::const_iterator it=p_actual->begin();
             it!=p_actual->end();++it) {
          if ((*it)->Weight()!=1.0) {
            msg_Debugging()<<"Add wt for "<<(**it)<<": "
                           <<(*it)->Weight(m_last[0]->KtStart())<<"\n";
            m_weight*=(*it)->Weight(m_last[0]->KtStart());
	    if (p_variationweights)
	      *p_variationweights*=(*it)->Weight(m_last[0]->KtStart());
            (*it)->Weights().clear();
          }
        }
      }
      ++nem;
      if (nem >= m_maxrewem) m_sudakov.SetVariationWeights(NULL);
      if (p_actual->NME()+nem>m_maxpart) return true;
      if (nem >= maxem) return true;
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
	m_lastcpl = m_sudakov.Selected()->Coupling()->Last();
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

void Shower::SetMS(const ATOOLS::Mass_Selector *const ms)
{
  m_sudakov.SetMS(ms);
  m_kinFF.SetMS(ms);
  m_kinFI.SetMS(ms);
  m_kinIF.SetMS(ms);
  m_kinII.SetMS(ms);
}

ATOOLS::Cluster_Amplitude *  Shower::SingletToAmplitude(const Singlet * sing, ATOOLS::Cluster_Amplitude * ampl){
  size_t is(0);
  for (Singlet::const_iterator it=sing->begin();it!=sing->end();++it){
      Parton *parton = *it;
      if(parton->GetType()==pst::IS){
          PartonToAmplitude(parton, ampl);
          is++;
      }
  }
  for (Singlet::const_iterator it=sing->begin();it!=sing->end();++it){
      Parton *parton = *it;
      if(parton->GetType()==pst::FS)   PartonToAmplitude(parton, ampl);
  }
  ampl->SetNIn(is);
  //CheckAmplitude(ampl);
  return ampl;
}


void Shower::PartonToAmplitude(const Parton *parton, Cluster_Amplitude * ampl){
  ampl->CreateLeg(parton->Momentum(), parton->GetFlavour(),parton->Col(),parton->Id());
  ATOOLS::Cluster_Leg * leg = ampl->Legs().back();
  leg->SetFromDec(parton->FromDec());
  if (parton->GetType()==pst::IS) leg->SetMom(-leg->Mom());
  if(ampl->KT2()==0){
      ampl->SetKT2(parton->KtStart());
    }
}

void Shower::CheckAmplitude(const ATOOLS::Cluster_Amplitude *ampl){
  ATOOLS::Vec4D check;
  for(ClusterLeg_Vector::const_iterator it=ampl->Legs().begin(); it !=ampl->Legs().end(); ++it){
      ATOOLS::Cluster_Leg *leg = *it;
      check+=leg->Mom();
  }
  msg_Debugging() << " mom sum: " << check << "\n";
  //problematic: hadron decays, where the decaying hadron does not show up in the amplitude since its not a parton

}

