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

Shower::Shower(PDF::ISR_Handler* isr, const int qcd, const int qed, int type)
    : p_actual(NULL), m_sudakov(isr, qcd, qed), p_isr(isr)
{
  Settings& s = Settings::GetMainSettings();
  const int evol{ s["CSS_EVOLUTION_SCHEME"].Get<int>() };
  int kfmode{ s["CSS_KFACTOR_SCHEME"].Get<int>() };
  const int scs{ s["CSS_SCALE_SCHEME"].Get<int>() };
  double k0sqf{ s["CSS_FS_PT2MIN"].Get<double>() };
  double k0sqi{ s["CSS_IS_PT2MIN"].Get<double>() };
  double gsplit_fac{ s["CSS_PT2MIN_GSPLIT_FACTOR"].Get<double>() };
  double fs_as_fac{ s["CSS_FS_AS_FAC"].Get<double>() };
  double is_as_fac{ s["CSS_IS_AS_FAC"].Get<double>() };
  double is_pdf_fac{ s["CSS_PDF_FAC"].Get<double>() };
  const double mth{ s["CSS_MASS_THRESHOLD"].Get<double>() };
  bool   forced_decays{ s["CSS_FORCED_DECAYS"].Get<bool>() };
  double forced_decays_gluon_scaling{ s["CSS_FORCED_GLUON_SCALING"].Get<double>() };
  m_reweight          = s["CSS_REWEIGHT"].Get<bool>();
  m_maxreweightfactor = s["CSS_MAX_REWEIGHT_FACTOR"].Get<double>();
  m_kscheme           = s["CSS_KIN_SCHEME"].Get<int>();
  m_recdec            = s["CSS_RECO_DECAYS"].Get<int>();
  m_maxpart           = s["CSS_MAXPART"].Get<int>();
  if (type) {
    kfmode=s["MI_CSS_KFACTOR_SCHEME"].Get<int>();
    k0sqf=s["MI_CSS_FS_PT2MIN"].Get<double>();
    k0sqi=s["MI_CSS_IS_PT2MIN"].Get<double>();
    gsplit_fac=s["MI_CSS_PT2MIN_GSPLIT_FACTOR"].Get<double>();
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
  m_sudakov.SetFacScaleFactor(is_pdf_fac);
  std::pair<double, double> pdfmin;
  pdfmin.first = s["CSS_PDF_MIN"].Get<double>();
  pdfmin.second = s["CSS_PDF_MIN_X"].Get<double>();
  m_sudakov.SetPDFMin(pdfmin);
  m_sudakov.InitSplittingFunctions(MODEL::s_model,kfmode);
  m_sudakov.SetCoupling(MODEL::s_model,k0sqi,k0sqf,is_as_fac,fs_as_fac,gsplit_fac);
  m_sudakov.SetReweightScaleCutoff(
      s["CSS_REWEIGHT_SCALE_CUTOFF"].Get<double>());
  m_sudakov.SetForcedHQDecays(forced_decays,forced_decays_gluon_scaling);
  m_kinFF.SetEvolScheme(evol-100*(evol/100));
  m_kinFI.SetEvolScheme(evol-100*(evol/100));
  m_kinIF.SetEvolScheme(evol/100);
  m_kinII.SetEvolScheme(evol/100);
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
  m_weightsmap.Clear();
  m_weightsmap["PS"] = Weights {Variations_Type::qcd};
  m_weightsmap["PS_QCUT"] = Weights {Variations_Type::qcut};
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
  if (!m_sudakov.CheckPDF(mom[0]/rpa->gen.PBunch(p->Beam())[0],p->GetFlavour(),p->Beam())) return -1;
  return p_remnants->GetRemnant(p->Beam())->TestExtract(p->GetFlavour(),mom)?1:-1;
}

int Shower::ReconstructDaughters(Singlet *const split,double &jcv,
				 Parton *const pi,Parton *const pj)
{
  if (split->GetSplit()) {
    if (split->GetSplit()->Stat()&part_status::code::decayed) {
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
  if (split==split->GetSing()->GetSplit()) {
    split->GetSing()->SetSplit(newpB);
  }
  if (split->GetSing()->GetSplit()==NULL ||
      (split->GetSing()->GetSplit()->Stat()&part_status::code::decayed)) {
    split->GetSing()->ArrangeColours(split,newpB,newpC);
  }
  if (newpB==split->GetSing()->GetSplit())
    split->GetSing()->SetSplit(split);
  newpB->SetPrev(split->GetPrev());
  newpC->SetPrev(split->GetPrev());
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
    if (split==split->GetSing()->GetSplit()) {
      split->GetSing()->SetSplit(newpB);
    }
    if (split->GetSing()->GetSplit()==NULL ||
        (split->GetSing()->GetSplit()->Stat()&part_status::code::decayed)) {
      split->GetSing()->RearrangeColours(split,newpB,newpC);
    }
    if (newpB==split->GetSing()->GetSplit())
      split->GetSing()->SetSplit(split);
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
  Parton *pi(new Parton((stype&1)?fla:flb,
                        split->ForcedDecay()?split->Momentum():split->LT()*split->Momentum(),
                        split->GetType()));
  if (stype&1) pi->SetBeam(split->Beam());
  if (stat==1) {
    if (split->GetType()==pst::IS &&
	RemnantTest(pi,NULL)==-1) stat=-1;
    if (split->GetSpect()->GetType()==pst::IS &&
	RemnantTest(split->GetSpect(),
		    split->GetType()==pst::IS?
		    &split->LT():NULL)==-1) stat=-1;
  }
  if (stat==-1) {
    split->SetMomentum(peo);
    spect->SetMomentum(pso);
    delete pj;
    delete pi;
    return stat;
  }
  pi->SetMass2(mi2);
  pi->SetSing(split->GetSing());
  pi->SetId(split->Id());
  pi->SetKScheme(split->KScheme());
  pi->SetKin(split->Kin());
  pj->SetKin(m_kscheme);
  pi->SetLT(split->LT());
  SetSplitInfo(peo,pso,split,pi,pj,stype);
  split->GetSing()->AddParton(pj);
  if (stype) split->GetSing()->BoostAllFS(pi,pj,spect,split->ForcedDecay());
  int ustat(UpdateDaughters(split,pi,pj,jcv,mode));
  if (ustat<=0 || (split->GetSing()->GetLeft() &&
		   !(split->GetSing()->GetSplit()->Stat()&part_status::code::decayed))) {
    if (stype) split->GetSing()->BoostBackAllFS(pi,pj,spect,split->ForcedDecay());
    delete pi;
    pj->DeleteAll();
    split->SetMomentum(peo);
    spect->SetMomentum(pso);
    return ustat;
  }
  const double split_weight {split->Weight()};
  m_weightsmap["PS"] *= split_weight;
  msg_Debugging() << "sw = " << split_weight
                  << ", w = " << m_weightsmap["PS"].Nominal() << "\n";
  if (m_reweight) {
    ATOOLS::Reweight(m_weightsmap["PS"],
                     [this, split](double varweight,
                                   QCD_Variation_Params& varparams) -> double {
                       return varweight * Reweight(&varparams, *split);
                     });
  }
  Singlet * singlet = split->GetSing();
  split->GetSing()->SplitParton(split,pi,pj);
  for (PLiter plit(singlet->begin());
       plit!=singlet->end();++plit)
    (*plit)->UpdateDaughters();
  return 1;
}

bool Shower::EvolveSinglet(Singlet * act,const size_t &maxem,size_t &nem)
{
  p_actual=act;
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
      (p_actual->GetSplit()->Stat()&part_status::code::fragmented) &&
      !(p_actual->GetSplit()->Stat()&part_status::code::decayed)) {
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
  m_sudakov.SetKeepReweightingInfo(m_reweight);
  while (true) {
    for (Singlet::const_iterator it=p_actual->begin();it!=p_actual->end();++it)
      if ((*it)->GetType()==pst::IS) SetXBj(*it);
    kt2win = 0.;
    Parton *split=SelectSplitting(kt2win);
    //no shower anymore
    if (split==NULL) {
      msg_Debugging()<<"No emission\n";
      ResetScales(p_actual->KtNext());
      for (Singlet::const_iterator it=p_actual->begin(); it!=p_actual->end();++it) {
        const double singlet_weight {(*it)->Weight()};
        if (singlet_weight != 1.0)
          msg_Debugging() << "Add wt for " << (**it) << ": " << singlet_weight
                          << "\n";
        m_weightsmap["PS"] *= singlet_weight;
        if (m_reweight) {
          ATOOLS::Reweight(
              m_weightsmap["PS"],
              [this, it](double varweight, QCD_Variation_Params& varparams)
                  -> double { return varweight * Reweight(&varparams, **it); });
          (*it)->SudakovReweightingInfos().clear();
        }
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
	   (p_actual->GetSplit()->Stat()&part_status::code::decayed))) {
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
      jcv = kstat ? jcv : -1.0;
      const bool is_jcv_positive {jcv >= 0.0};
      bool all_vetoed {true};
      const int nqcuts = s_variations->Size(Variations_Type::qcut);
      std::vector<bool> skips (nqcuts + 1, false);
      int nskips {0};
      ATOOLS::ReweightAll(
          m_weightsmap["PS_QCUT"],
          [this, jcv, is_jcv_positive, &all_vetoed, &skips, &nskips](
              double varweight,
              size_t varindex,
              Qcut_Variation_Params* qcutparams) -> double {
            msg_Debugging()
                << "Applying veto weight to qcut var #" << varindex << " {\n";
            bool stat {is_jcv_positive};
            if (stat && p_actual->JF()) {
              const double fac {
                  qcutparams == nullptr ? 1.0 : qcutparams->m_scale_factor};
              stat = jcv < sqr(p_actual->JF()->Qcut() * fac);
              msg_Debugging() << "  jcv = " << sqrt(jcv) << " vs "
                              << p_actual->JF()->Qcut() << " * " << fac << " = "
                              << p_actual->JF()->Qcut() * fac << "\n";
            }
            if (stat) {
              msg_Debugging() << "} no jet veto\n";
              all_vetoed = false;
              return varweight;
            } else if (p_actual->NLO() & 2) {
              msg_Debugging() << "  skip emission\n";
              skips[varindex] = true;
              ++nskips;
              all_vetoed = false;
              msg_Debugging() << "} no jet veto\n";
              return varweight;
            } else {
              msg_Debugging() << "} jet veto\n";
              return 0.0;
            }
          });
      if (p_actual->NLO()&2) {
        const int nqcdvars = s_variations->Size(Variations_Type::qcd);
        if (skips[0])
          nskips += nqcdvars;
        const double wskip {nskips / double(nqcuts + nqcdvars + 1)};
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
          const double fac {1.0 / lkf / wskip};
          m_weightsmap["PS"] *= fac;
          m_weightsmap["PS_QCUT"] *= skips;
          continue;
        }
	else {
          const double fac {1.0 / (1.0 - wskip)};
          skips.flip();
          m_weightsmap["PS"] *= fac;
          m_weightsmap["PS_QCUT"] *= skips;
        }
      }
      if (all_vetoed) return false;
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
            const double singlet_weight {(*it)->Weight(m_last[0]->KtStart())};
            msg_Debugging()
                << "Add wt for " << (**it) << ": " << singlet_weight << "\n";
            m_weightsmap["PS"] *= singlet_weight;
            (*it)->Weights().clear();
          }
          if (m_reweight) {
            ATOOLS::Reweight(
                m_weightsmap["PS"],
                [this, it](double varweight,
                           QCD_Variation_Params& varparams) -> double {
                  return varweight * Reweight(&varparams, **it);
                });
            (*it)->SudakovReweightingInfos().clear();
          }
        }
      }
      ++nem;
      if (p_actual->NME()+nem>m_maxpart || nem >= maxem) {
	return true;
      }
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
    if (m_sudakov.Generate(split,kt2win)) {
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
	if (split->ForcedDecay()) {
	  split->IncForcedTrials();
	  if (split->ForcedTrials()>10) {
	    split->SetKtStart(0.0);
	    return false;
	  }
	}
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

double Shower::Reweight(QCD_Variation_Params* varparams,
                        Parton& splitter)
{
  const double kt2win {(m_last[0] == NULL) ? 0.0 : m_last[0]->KtStart()};
  Sudakov_Reweighting_Infos& infos = splitter.SudakovReweightingInfos();
  double overallrewfactor {1.0};

  for (auto info : infos) {

    // do not reweighting trial emissions below the scale of the accepted
    // emission
    // NOTE: contrary to what one would expect, the infos are not strictly
    // ordered descending in pT, which means that we can not use `break` here,
    // instead we must use `continue`
    if (info.scale < kt2win)
      continue;

    const double rejwgt {1.0 - info.accwgt};
    double rewfactor {1.0};
    double accrewfactor {1.0};

    // perform PDF reweighting
    // NOTE: also the Jacobians depend on the Running_AlphaS class, but only
    // through the number of flavours, which should not vary between AlphaS
    // variations anyway; therefore we do not insert AlphaS for the PDF
    // reweighting
    const cstp::code type {info.sf->GetType()};
    if (type == cstp::II || type == cstp::FI || type == cstp::IF) {
      // insert new PDF
      const Flavour swappedflspec {info.sf->Lorentz()->FlSpec()};
      info.sf->Lorentz()->SetFlSpec(info.flspec);
      PDF::PDF_Base** swappedpdf {info.sf->PDF()};
      PDF::PDF_Base* pdf[] = {varparams->p_pdf1, varparams->p_pdf2};
      info.sf->SetPDF(pdf);
      // calculate new J
      const double lastJ(info.sf->Lorentz()->LastJ());
      double newJ;
      double muF2fac {1.0};
      if (varparams->m_showermuF2enabled)
	muF2fac = varparams->m_muF2fac;
      switch (type) {
      case cstp::II:
        newJ = info.sf->Lorentz()->JII(
            info.z, info.y, info.x, muF2fac * info.scale);
        break;
      case cstp::IF:
        newJ = info.sf->Lorentz()->JIF(
            info.z, info.y, info.x, muF2fac * info.scale);
        break;
      case cstp::FI:
        newJ = info.sf->Lorentz()->JFI(
            info.y, info.x, muF2fac * info.scale);
        break;
      case cstp::FF:
      case cstp::none:
        THROW(fatal_error, "Unexpected splitting configuration");
      }
      // clean up
      info.sf->SetPDF(swappedpdf);
      info.sf->Lorentz()->SetLastJ(lastJ);
      info.sf->Lorentz()->SetFlSpec(swappedflspec);
      // validate and apply
      if (newJ == 0.0) {
        varparams->IncrementOrInitialiseWarningCounter(
            "different PDF cut-off");
        continue;
      } else {
        const double pdfrewfactor {newJ / info.lastj};
        if (pdfrewfactor < 0.25 || pdfrewfactor > 4.0) {
          varparams->IncrementOrInitialiseWarningCounter(
              "large PDF reweighting factor");
        }
        accrewfactor *= pdfrewfactor;
      }
    }

    // AlphaS reweighting
    if (info.sf->Coupling()->AllowsAlternativeCouplingUsage()) {
      // insert new AlphaS
      const double lastcpl {info.sf->Coupling()->Last()};
      double muR2fac {1.0};
      if (varparams->m_showermuR2enabled)
	muR2fac = varparams->m_muR2fac;
      info.sf->Coupling()->SetAlternativeUnderlyingCoupling(
          varparams->p_alphas, muR2fac);
      // calculate new coupling
      double newcpl {info.sf->Coupling()->Coupling(info.scale, 0)};
      // clean up
      info.sf->Coupling()->SetAlternativeUnderlyingCoupling(nullptr);
      info.sf->Coupling()->SetLast(lastcpl);
      // validate and apply
      const double alphasrewfactor {newcpl / info.lastcpl};
      if (alphasrewfactor < 0.5 || alphasrewfactor > 2.0) {
        varparams->IncrementOrInitialiseWarningCounter(
            "large AlphaS reweighting factor");
      }
      accrewfactor *= alphasrewfactor;
    }

    // calculate and apply overall factor
    if (info.accepted) {
      rewfactor = accrewfactor;
    } else {
      rewfactor = 1.0 + (1.0 - accrewfactor) * (1.0 - rejwgt) / rejwgt;
    }
    overallrewfactor *= rewfactor;
  }

  // guard against gigantic accumulated reweighting factors
  if (std::abs(overallrewfactor) > m_maxreweightfactor) {
    msg_Debugging() << "Veto large CSS Sudakov reweighting factor for parton: "
                    << splitter;
    varparams->IncrementOrInitialiseWarningCounter(
        "vetoed large reweighting factor for parton");
    return 1.0;
  }

  return overallrewfactor;
}

void Shower::SetMS(const ATOOLS::Mass_Selector *const ms)
{
  m_sudakov.SetMS(ms);
  m_kinFF.SetMS(ms);
  m_kinFI.SetMS(ms);
  m_kinIF.SetMS(ms);
  m_kinII.SetMS(ms);
}
