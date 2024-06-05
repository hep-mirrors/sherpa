#include "MCATNLO/Showers/Shower.H"
#include "MCATNLO/Tools/Parton.H"
#include "MCATNLO/Main/CS_Gamma.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
#include "REMNANTS/Main/Remnant_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Phys/Cluster_Leg.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace MCATNLO;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

Shower::Shower(PDF::ISR_Handler * isr,const int qed) :
  p_actual(NULL), m_sudakov(isr,qed), p_isr(isr)
{
  auto pss = Settings::GetMainSettings()["SHOWER"];
  auto nlopss = Settings::GetMainSettings()["MC@NLO"];
  const int evol{ pss["EVOLUTION_SCHEME"].Get<int>() };
  const int kfmode{ pss["KFACTOR_SCHEME"].Get<int>() };
  const int scs{ pss["SCALE_SCHEME"].Get<int>() };
  const double k0sqf{ pss["FS_PT2MIN"].Get<double>() };
  const double k0sqi{ pss["IS_PT2MIN"].Get<double>() };
  const double gsplit_fac{ pss["PT2MIN_GSPLIT_FACTOR"].Get<double>() };
  const double fs_as_fac{ pss["FS_AS_FAC"].Get<double>() };
  const double is_as_fac{ pss["IS_AS_FAC"].Get<double>() };
  const double is_pdf_fac{ pss["PDF_FAC"].Get<double>() };
  const double mth{ pss["MASS_THRESHOLD"].Get<double>() };
  m_reweight = pss["REWEIGHT"].Get<bool>();
  m_maxreweightfactor = pss["MAX_REWEIGHT_FACTOR"].Get<double>();
  m_kscheme = nlopss["KIN_SCHEME"].Get<int>();
  std::vector<size_t> disallowflavs{
    nlopss["DISALLOW_FLAVOUR"].GetVector<size_t>() };
  m_sudakov.SetShower(this);
  m_sudakov.SetMassThreshold(mth);
  m_sudakov.SetScaleScheme(scs);
  m_sudakov.SetFacScaleFactor(is_pdf_fac);
  std::pair<double, double> pdfmin;
  pdfmin.first = pss["PDF_MIN"].Get<double>();
  pdfmin.second = pss["PDF_MIN_X"].Get<double>();
  m_sudakov.SetPDFMin(pdfmin);
  m_sudakov.SetDisallowFlavour(disallowflavs);
  m_sudakov.InitSplittingFunctions(MODEL::s_model,kfmode);
  m_sudakov.SetCoupling(MODEL::s_model,k0sqi,k0sqf,is_as_fac,fs_as_fac,gsplit_fac);
  m_sudakov.SetReweightScaleCutoff(
      pss["REWEIGHT_SCALE_CUTOFF"].Get<double>());
  m_kinFF.SetSudakov(&m_sudakov);
  m_kinFI.SetSudakov(&m_sudakov);
  m_kinIF.SetSudakov(&m_sudakov);
  m_kinII.SetSudakov(&m_sudakov);
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

int Shower::SetXBj(Parton *const p) const
{
  double x(p_isr->CalcX(p->Momentum()));
  if (x>1.0) return -1;
  p->SetXbj(x);
  return 1;
}

int Shower::RemnantTest(Parton *const p)
{
  if (p->Momentum()[0]<0.0 || p->Momentum().Nan()) return -1;
  double x(p_isr->CalcX(p->Momentum()));
  if (x>1.0 && !IsEqual(x,1.0,1.0e-6)) return -1;
  return p_isr->GetRemnant(p->Beam())->
    TestExtract(p->GetFlavour(),p->Momentum())?1:-1;
}

int Shower::UpdateDaughters(Parton *const split,Parton *const newpB,
			    Parton *const newpC,double &jcv)
{
  newpB->SetStart(split->KtTest());
  newpC->SetStart(split->KtTest());
  newpB->SetKtMax(split->KtMax());
  newpC->SetKtMax(split->KtMax());
  newpB->SetVeto(split->KtVeto());
  newpC->SetVeto(split->KtVeto());
  newpB->SetSing(split->GetSing());
  newpC->SetSing(split->GetSing());
  newpB->SetMEFlow(1,split->Color().m_i[0]);
  newpB->SetMEFlow(2,split->Color().m_i[1]);
  newpC->SetMEFlow(1,split->Color().m_j[0]);
  newpC->SetMEFlow(2,split->Color().m_j[1]);
  int scol[2]={split->GetSpect()->GetMEFlow(1),
	       split->GetSpect()->GetMEFlow(2)};
  split->GetSpect()->SetMEFlow(1,split->Color().m_k[0]);
  split->GetSpect()->SetMEFlow(2,split->Color().m_k[1]);
  split->GetSing()->ArrangeColours(split,newpB,newpC);
  int rd(1);
  if (rd) {
    if (newpB->GetType()==pst::IS &&
	RemnantTest(newpB)==-1) rd=-1;
    if (split->GetSpect()->GetType()==pst::IS &&
	RemnantTest(split->GetSpect())==-1) rd=-1;
  }
  int sci[2]={split->GetFlow(1),split->GetFlow(2)};
  int scim[2]={split->GetMEFlow(1),split->GetMEFlow(2)};
  m_flav=split->GetFlavour();
  split->SetFlavour(newpB->GetFlavour());
  split->SetFlow(1,newpB->GetFlow(1));
  split->SetFlow(2,newpB->GetFlow(2));
  split->SetMEFlow(1,newpB->GetMEFlow(1));
  split->SetMEFlow(2,newpB->GetMEFlow(2));
  if (rd==1) rd=p_gamma->Reject()?-1:1;
  const double gamma_weight {p_gamma->Weight()};
  m_weightsmap["MC@NLO_PS"] *= gamma_weight;
  if (rd==1 && split->KtTest()>split->KtMax())
    jcv=split->GetSing()->JetVeto(&m_sudakov);
  split->SetFlavour(m_flav);
  split->SetFlow(1,sci[0]);
  split->SetFlow(2,sci[1]);
  split->SetMEFlow(1,scim[0]);
  split->SetMEFlow(2,scim[1]);
  if (rd==0) {
    split->GetSpect()->SetMEFlow(1,scol[0]);
    split->GetSpect()->SetMEFlow(2,scol[1]);
  }
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
  p_old[0]->SetFlav((mode&1)?split->GetFlavour().Bar():split->GetFlavour());
  p_old[0]->SetCol(ColorID(split->GetFlow((mode&1)?2:1),
			   split->GetFlow((mode&1)?1:2)));
  m_last[0]=newb;
  m_last[1]=newc;
  m_last[2]=split->GetSpect();
  m_last[3]=split;
}

int Shower::MakeKinematics
(Parton *split,const Flavour &fla,const Flavour &flb,
 const Flavour &flc,double &jcv)
{
  DEBUG_FUNC("");
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
  SetSplitInfo(peo,pso,split,pi,pj,stype);
  split->GetSing()->push_back(pj);
  if (stype) split->GetSing()->BoostAllFS(split,pj,spect);
  Flavour fls(split->GetFlavour());
  int ustat(UpdateDaughters(split,pi,pj,jcv));
  split->GetSing()->pop_back();
  if (ustat<=0) {
    split->SetFlavour(fls);
    if (stype) split->GetSing()->BoostBackAllFS(split,pj,spect);
    delete pi;
    delete pj;
    msg_Debugging()<<"Save history for\n"<<*split<<*spect<<"\n";
    split->SetMomentum(peo);
    spect->SetMomentum(pso);
    return ustat;
  }
  if (m_reweight) {
    ATOOLS::Reweight(m_weightsmap["MC@NLO_PS"],
                     [this, split](double varweight,
                                   QCD_Variation_Params& varparams) -> double {
                       return varweight * Reweight(&varparams, *split);
                     });
  }
  split->GetSing()->SplitParton(split,pi,pj);
  return 1;
}

bool Shower::EvolveShower(Singlet *act,const size_t &maxem,size_t &nem)
{
  m_weightsmap.Clear();
  m_weightsmap["MC@NLO_PS"] = Weights {Variations_Type::qcd};
  m_weightsmap["MC@NLO_QCUT"] = Weights {Variations_Type::qcut};
  p_actual=act;
  Parton * split;
  Vec4D mom;

  if (nem>=maxem) return true;

  m_sudakov.SetKeepReweightingInfo(m_reweight);

  while (true) {
    m_last[0]=m_last[1]=m_last[2]=NULL;
    double kt2win = 0.;
    split = SelectSplitting(kt2win);
    //no shower anymore 
    if (split==NULL) {
      if (m_reweight) {
        for (Singlet::const_iterator it = p_actual->begin();
             it != p_actual->end();
             ++it) {
          ATOOLS::Reweight(
              m_weightsmap["MC@NLO_PS"],
              [this, it](double varweight, QCD_Variation_Params& varparams)
                  -> double { return varweight * Reweight(&varparams, **it); });
          (*it)->SudakovReweightingInfos().clear();
        }
      }
      return true;
    }
    else {
      msg_Debugging()<<"Emission "<<m_flavA<<" -> "<<m_flavB<<" "<<m_flavC
		     <<" at kt = "<<sqrt(split->KtTest())
		     <<", z = "<<split->ZTest()<<", y = "
		     <<split->YTest()<<" for\n"<<*split
		     <<*split->GetSpect()<<"\n";
      m_last[0]=m_last[1]=m_last[2]=m_last[3]=NULL;
      ResetScales(split);
      double jcv(0.0);
      int kstat(MakeKinematics(split,m_flavA,m_flavB,m_flavC,jcv));
      msg_Debugging()<<"stat = "<<kstat<<"\n";
      if (kstat<0) continue;
      jcv = kstat ? jcv : -1.0;
      const bool is_jcv_positive {jcv >= 0.0};
      bool all_vetoed {true};
      ATOOLS::ReweightAll(
          m_weightsmap["MC@NLO_QCUT"],
          [this, jcv, is_jcv_positive, &all_vetoed](
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
            } else {
              msg_Debugging() << "} jet veto\n";
              return 0.0;
            }
          });
      if (all_vetoed)
        return false;
      msg_Debugging()<<"nem = "<<nem+1<<" vs. maxem = "<<maxem<<"\n";
      if (m_reweight) {
        for (Singlet::const_iterator it = p_actual->begin();
             it != p_actual->end();
             ++it) {
          ATOOLS::Reweight(
              m_weightsmap["MC@NLO_PS"],
              [this, it](double varweight, QCD_Variation_Params& varparams)
                  -> double { return varweight * Reweight(&varparams, **it); });
          (*it)->SudakovReweightingInfos().clear();
        }
      }
      if (++nem>=maxem) return true;
    }
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
      split->SetCol(m_sudakov.GetCol());
      split->SetTest(kt2,z,y,phi);
      return true;
    }
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
      switch (type) {
      case cstp::II:
        newJ = info.sf->Lorentz()->JII(
            info.z, info.y, info.x, varparams->m_muF2fac * info.scale, nullptr);
        break;
      case cstp::IF:
        newJ = info.sf->Lorentz()->JIF(
            info.z, info.y, info.x, varparams->m_muF2fac * info.scale, nullptr);
        break;
      case cstp::FI:
        newJ = info.sf->Lorentz()->JFI(
            info.y, info.x, varparams->m_muF2fac * info.scale, nullptr);
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
            "MCatNLO different PDF cut-off");
        continue;
      } else {
        const double pdfrewfactor {newJ / info.lastj};
        if (pdfrewfactor < 0.25 || pdfrewfactor > 4.0) {
          varparams->IncrementOrInitialiseWarningCounter(
              "MCatNLO large PDF reweighting factor");
        }
        accrewfactor *= pdfrewfactor;
      }
    }

    // AlphaS reweighting
    if (info.sf->Coupling()->AllowsAlternativeCouplingUsage()) {
      // insert new AlphaS
      const double lastcpl {info.sf->Coupling()->Last()};
      info.sf->Coupling()->SetAlternativeUnderlyingCoupling(
          varparams->p_alphas, varparams->m_muR2fac);
      // calculate new coupling
      double newcpl {info.sf->Coupling()->Coupling(info.scale, 0, nullptr)};
      // clean up
      info.sf->Coupling()->SetAlternativeUnderlyingCoupling(nullptr);
      info.sf->Coupling()->SetLast(lastcpl);
      // validate and apply
      const double alphasrewfactor {newcpl / info.lastcpl};
      if (alphasrewfactor < 0.5 || alphasrewfactor > 2.0) {
        varparams->IncrementOrInitialiseWarningCounter(
            "MCatNLO large AlphaS reweighting factor");
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
    msg_Debugging() << "Veto large MC@NLO Sudakov reweighting factor for parton: "
                    << splitter;
    varparams->IncrementOrInitialiseWarningCounter(
        "MCatNLOvetoed large reweighting factor for parton");
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
