#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Phys/Variations.H"

namespace PHASIC {

  class H_KFactor: public KFactor_Setter_Base {
  protected:

    double m_k0sq[2];
    int    m_fomode, m_mtmode, m_kfmode;

    PHASIC::Multi_Channel *p_fsmc;

  public:

    H_KFactor(const KFactor_Setter_Arguments &args);

  };// end of class H_KFactor

  class HNNLO_KFactor: public H_KFactor {
  public:

    inline HNNLO_KFactor(const KFactor_Setter_Arguments &args):
      H_KFactor(args) {}

    double KFactor(const int mode = 0);
    double KFactor(ATOOLS::QCD_Variation_Params* params, const int& mode);

  };// end of class HNNLO_KFactor

  class HHF1_KFactor: public H_KFactor {
  public:

    inline HHF1_KFactor(const KFactor_Setter_Arguments &args):
      H_KFactor(args) {}

    double KFactor(const int mode = 0);
    double KFactor(ATOOLS::QCD_Variation_Params* params, const int& mode);

  };// end of class HHF1_KFactor

  class HHF2_KFactor: public H_KFactor {
  public:

    inline HHF2_KFactor(const KFactor_Setter_Arguments &args):
      H_KFactor(args) {}

    double KFactor(const int mode = 0);
    double KFactor(ATOOLS::QCD_Variation_Params* params, const int& mode);

  };// end of class HHF2_KFactor

  class HNLO_KFactor: public H_KFactor {
  public:

    inline HNLO_KFactor(const KFactor_Setter_Arguments &args):
      H_KFactor(args) {}

    double KFactor(const int mode = 0);
    double KFactor(ATOOLS::QCD_Variation_Params* params, const int& mode);

  };// end of class HNLO_KFactor

  class HF1_KFactor: public H_KFactor {
  public:

    inline HF1_KFactor(const KFactor_Setter_Arguments &args):
      H_KFactor(args) {}

    double KFactor(const int mode = 0);
    double KFactor(ATOOLS::QCD_Variation_Params* params, const int& mode);

  };// end of class HF1_KFactor

}// end of namespace PHASIC

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#include "QT_Selector.H"
#include "coeffqt.H"
#include "hard.H"
#include "higgsfullsm.H"
#include "Tools.H"

using namespace SHNNLO;
using namespace PHASIC;
using namespace ATOOLS;

H_KFactor::H_KFactor
(const KFactor_Setter_Arguments &args):
  KFactor_Setter_Base(args)
{
  Settings& s = Settings::GetMainSettings();
  auto pss = s["SHOWER"];
  if (p_proc->Generator()->Name()!="Amegic")
    THROW(inconsistent_option,"Must use 'ME_Generator Amegic;'");
  if (s_pdf==NULL) {
    s_as=MODEL::as;
    s_pdf=p_proc->Integrator()->ISR()->PDF(0);
    s_pdfmin[0] = pss["PDF_MIN"].Get<double>();
    s_pdfmin[1] = pss["PDF_MIN_X"].Get<double>();
  }
  m_k0sq[0] = pss["FS_PT2MIN"].Get<double>();
  m_k0sq[1] = pss["IS_PT2MIN"].Get<double>();
  m_fomode = s["NNLOqT_FOMODE"].Get<int>();
  m_mtmode = s["HNNLO_MTOP_MODE"].Get<int>();
  m_kfmode = s["HNNLO_KF_MODE"].Get<int>();
  // initialize constants
  SHNNLO::Mt=Flavour(kf_t).Mass(true);
  SHNNLO::Mb=Flavour(kf_b).Mass(true);
  SHNNLO::Mc=Flavour(kf_c).Mass(true);
  nf=Flavour(kf_jet).Size()/2;
  updateparam();
  if (p_proc->Integrator()->PSHandler()!=NULL) {
    p_fsmc=p_proc->Integrator()->PSHandler()->FSRIntegrator();
    p_fsmc->AddERan("zeta_1");
    p_fsmc->AddERan("zeta_2");
    p_fsmc->AddERan("zeta_1'");
    p_fsmc->AddERan("zeta_2'");
    p_fsmc->AddERan("zeta_1''");
    p_fsmc->AddERan("zeta_2''");
  }
}

double HNNLO_KFactor::KFactor(const int mode) 
{
  DEBUG_FUNC(p_proc->Name()<<" "<<p_proc->Generator()->Name()<<" "<<mode);
  const int lmode(mode&~2);
  m_weight = KFactor(NULL, lmode);
  msg_Debugging()<<"Weight: "<<m_weight<<"\n";
  if (s_variations->Size()) {
    std::vector<double> &bkw(p_proc->GetMEwgtinfo()->m_bkw);
    if (mode&2) bkw.clear();
    size_t oldsize(bkw.size());
    s_variations->ForEach(
        [this, &lmode](size_t varindex,
                       QCD_Variation_Params& varparams) -> void {
          KFactor(&varparams, lmode);
        });
    msg_Debugging()<<"New K factors: "<<std::vector<double>
      (&bkw[oldsize],&bkw.back()+1)<<"\n";
    for (size_t i(oldsize);i<bkw.size();++i) bkw[i]*=m_weight?1.0/m_weight:0.0;
    msg_Debugging()<<"Weight variations: "<<bkw<<"\n";
  }
  return m_weight;
}

// NNLO K factor
double HNNLO_KFactor::KFactor(QCD_Variation_Params* params, const int& mode)
{
  if (params==NULL) {
    s_as=MODEL::as;
    s_pdf=p_proc->Integrator()->ISR()->PDF(0);
  }
  else {
    s_as=params->p_alphas;
    s_pdf=params->p_pdf1;
  }
  const Vec4D_Vector &p(p_proc->Integrator()->Momenta());
  Scale_Setter_Base *sc(p_proc->ScaleSetter());
  double mur(sqrt(sc->Scale(stp::ren)*(params?params->m_muR2fac:1.0))), Q(p[2].Mass());
  double muf(sqrt(sc->Scale(stp::fac)*(params?params->m_muF2fac:1.0)));
  double norm(m_mtmode?ggH1l(Q,Mt,Mb,Mc):1.);
  double as4pi((*s_as)(mur*mur)/(4.0*M_PI));
  if (m_kfmode&1) {
    norm*=sqr(1.0+as4pi*hf1tt(mur,0.,Mt)+sqr(as4pi)*hf2tt(mur,0.,Mt));
  }
  if (m_mtmode && 
     ((mode==1 && p_proc->NOut()>1) ||         // for B of H+1j
       p_proc->Info().Has(nlo_type::real)) ) { // for R of H+0j
    // add 1j mass dependence here
    Flavour fl0(p_proc->Flavours()[0]);
    Flavour fl1(p_proc->Flavours()[1]);
    int i0(fl0.IsGluon()?0:(long int)fl0);
    int i1(fl1.IsGluon()?0:(long int)fl1);
    Mh=Q;
    double oldnorm(norm);
    if (i0==0&&i1==0) norm*=ggHg1l(2.*(p[0]*p[1]),-2.*(p[0]*p[3]),-2.*(p[1]*p[3]),Mh,Mt,Mb,Mc);
    if (i0!=0&&i1!=0) norm*=qqHg1l(2.*(p[0]*p[1]),-2.*(p[0]*p[3]),-2.*(p[1]*p[3]),Mh,Mt,Mb,Mc);
    if (i0==0&&i1!=0) norm*=qqHg1l(-2.*(p[1]*p[3]),2.*(p[0]*p[1]),-2.*(p[0]*p[3]),Mh,Mt,Mb,Mc);
    if (i0!=0&&i1==0) norm*=qqHg1l(-2.*(p[0]*p[3]),-2.*(p[1]*p[3]),2.*(p[0]*p[1]),Mh,Mt,Mb,Mc);
    if (IsBad(norm)) norm=oldnorm;
    else norm/=ggH1l(Q,Mt,Mb,Mc);
  }
  if (p_proc->NOut()>
      (p_proc->Info().Has(nlo_type::real)?2:1)) {
    double weight(NNLODiffWeight(p_proc,norm,mur*mur,muf*muf,m_k0sq,mode,m_fomode,0,
				 params?params->Name():""));
    if (params) p_proc->GetMEwgtinfo()->m_bkw.push_back(weight);
    return params?1.0:weight;
  }
  // 1+2*corr1 reweights VIRS part of H+0jet @ NLO
  double corr1(0.0);
  if (!(m_kfmode&2)) corr1=as4pi*hf1tt(mur,0.,Mt);
  if (p_proc->Info().Has(nlo_type::real)) {
    double weight(NNLODeltaWeight(p_proc,(1+2.0*corr1)*norm,m_fomode));
    if (params) p_proc->GetMEwgtinfo()->m_bkw.push_back(weight);
    return params?1.0:weight;
  }
  Cluster_Amplitude *ampl(NULL);
  if (sc->Amplitudes().size()) ampl=sc->Amplitudes().front();
  if (ampl) ampl->SetNLO(4);
  if (mode!=1) {
    if (params) p_proc->GetMEwgtinfo()->m_bkw.push_back((1+2.0*corr1)*norm);
    return params?1.0:(1+2.0*corr1)*norm;
  }
  QT_Selector *jf=(QT_Selector*)
    p_proc->Selector()->GetSelector("NNLOqT_Selector");
  if (jf==NULL) THROW(fatal_error,"Must use selector \"NNLOqT\"");
  DEBUG_VAR(muf);
  DEBUG_VAR(Q);
  int swap(p[0][3]<p[1][3]);
  double x0(swap?p[0].PMinus()/rpa->gen.PBeam(1).PMinus():
	    p[0].PPlus()/rpa->gen.PBeam(0).PPlus());
  double x1(swap?p[1].PPlus()/rpa->gen.PBeam(0).PPlus():
	    p[1].PMinus()/rpa->gen.PBeam(1).PMinus());
  double z0(p_fsmc->ERan("zeta_1")), z1(p_fsmc->ERan("zeta_2"));
  s_z[0]=p_fsmc->ERan("zeta_1'");
  s_z[1]=p_fsmc->ERan("zeta_2'");
  s_z[2]=p_fsmc->ERan("zeta_1''");
  s_z[3]=p_fsmc->ERan("zeta_2''");
  double lnRF=(mur==muf?0.:log(mur/muf));
  double H1=2.*hf1tt(mur,0.,Mt);
  double H2=2.*hf2tt(mur,0.,Mt)+sqr(H1/2.);
  double K0(0.),K1(0.),K2(0.);
  K2+=Cgg2gg(x0,x1,z0,z1,jf->QTCut(),muf,Q);
  K2+=Cgg2qq(x0,x1,z0,z1,jf->QTCut(),muf,Q);
  K2+=Cgg2gq(x0,x1,z0,z1,jf->QTCut(),muf,Q)+Cgg2gq(x1,x0,z1,z0,jf->QTCut(),muf,Q);
  if (mur!=muf) { 
    K1+=Cgg1gg(x0,x1,z0,z1,jf->QTCut(),muf,Q);
    K1+=Cgg1gq(x0,x1,z0,z1,jf->QTCut(),muf,Q)+Cgg1gq(x1,x0,z1,z0,jf->QTCut(),muf,Q);
  }
  K0=Cgg0gg(x0,x1,z0,z1,jf->QTCut(),muf,Q);
  // based on hf1tt(muR)=hf1tt(muF), hf2tt(muR)=hf2tt(muF)+2*(beta0*hf1tt(muF)-beta1)*log(muR/muF)
  double weight=(K2+4.*beta0*lnRF*K1)/K0+4.*sqr(beta0*lnRF)+6.*beta1*lnRF;
  if (!(m_kfmode&2)) {
    weight+=H2-sqr(H1/2.);
  }
  if (m_kfmode==0) {
    weight-=3.0*sqr(H1/2.); // remove double counting because 2*hf1tt is multiplied by VIRS part of H+0jet @ NLO
  }
  if (m_kfmode&4) {
    weight-=hf2gg(muf,muf);
  }
  // weight reweights B part of H+0jet @ NLO
  weight=(1.+sqr(as4pi)*weight)*norm;
  // add two-loop ggH top mass dependence; note that hf0tt(mur,Q,Mt)==ggH1l(Q,Mt,0.,0.)
  if (m_mtmode) weight+=as4pi*(2.*hf1tt(mur,Q,Mt)-H1)*hf0tt(mur,Q,Mt); 
  DEBUG_VAR(weight);
  if (params) p_proc->GetMEwgtinfo()->m_bkw.push_back(weight);
  return params?1.0:weight;
}

double HHF1_KFactor::KFactor(const int mode) 
{
  DEBUG_FUNC(p_proc->Name()<<" "<<p_proc->Generator()->Name()<<" "<<mode);
  const int lmode(mode&~2);
  m_weight = KFactor(NULL, lmode);
  msg_Debugging()<<"Weight: "<<m_weight<<"\n";
  if (s_variations->Size()) {
    std::vector<double> &bkw(p_proc->GetMEwgtinfo()->m_bkw);
    if (mode&2) bkw.clear();
    size_t oldsize(bkw.size());
    s_variations->ForEach(
        [this, &lmode](size_t varindex,
                       QCD_Variation_Params& varparams) -> void {
          KFactor(&varparams, lmode);
        });
    msg_Debugging()<<"New K factors: "<<std::vector<double>
      (&bkw[oldsize],&bkw.back()+1)<<"\n";
    for (size_t i(oldsize);i<bkw.size();++i) bkw[i]*=m_weight?1.0/m_weight:0.0;
    msg_Debugging()<<"Weight variations: "<<bkw<<"\n";
  }
  return m_weight;
}

// for applying Higgs effective coupling in MC@NLO style
// provide 2*hf1tt to multiple a standalone Higgs NLO
double HHF1_KFactor::KFactor(QCD_Variation_Params* params, const int& mode)
{
  if (params==NULL) {
    s_as=MODEL::as;
    s_pdf=p_proc->Integrator()->ISR()->PDF(0);
  }
  else {
    s_as=params->p_alphas;
    s_pdf=params->p_pdf1;
  }
  const Vec4D_Vector &p(p_proc->Integrator()->Momenta());
  Scale_Setter_Base *sc(p_proc->ScaleSetter());
  double mur(sqrt(sc->Scale(stp::ren)*(params?params->m_muR2fac:1.0))), Q(p[2].Mass());
  double norm(m_mtmode?ggH1l(Q,Mt,Mb,Mc):1.);
  double K(0.0);
  if (!(m_kfmode&1)&&(m_kfmode&2)) K+=2.0*hf1tt(mur,0.,Mt);
  K*=(*s_as)(mur*mur)/(4.0*M_PI);
  if (m_mtmode && p_proc->Info().Has(nlo_type::real) ) { // for R of H+0j
    // add 1j mass dependence here
    Flavour fl0(p_proc->Flavours()[0]);
    Flavour fl1(p_proc->Flavours()[1]);
    int i0(fl0.IsGluon()?0:(long int)fl0);
    int i1(fl1.IsGluon()?0:(long int)fl1);
    Mh=Q;
    double oldnorm(norm);
    if (i0==0&&i1==0) norm*=ggHg1l(2.*(p[0]*p[1]),-2.*(p[0]*p[3]),-2.*(p[1]*p[3]),Mh,Mt,Mb,Mc);
    if (i0!=0&&i1!=0) norm*=qqHg1l(2.*(p[0]*p[1]),-2.*(p[0]*p[3]),-2.*(p[1]*p[3]),Mh,Mt,Mb,Mc);
    if (i0==0&&i1!=0) norm*=qqHg1l(-2.*(p[1]*p[3]),2.*(p[0]*p[1]),-2.*(p[0]*p[3]),Mh,Mt,Mb,Mc);
    if (i0!=0&&i1==0) norm*=qqHg1l(-2.*(p[0]*p[3]),-2.*(p[1]*p[3]),2.*(p[0]*p[1]),Mh,Mt,Mb,Mc);
    if (IsBad(norm)) norm=oldnorm;
    else norm/=ggH1l(Q,Mt,Mb,Mc);
  }
  if (params) p_proc->GetMEwgtinfo()->m_bkw.push_back(K*norm);
  return params?1.0:K*norm;
}

double HHF2_KFactor::KFactor(const int mode) 
{
  DEBUG_FUNC(p_proc->Name()<<" "<<p_proc->Generator()->Name()<<" "<<mode);
  const int lmode(mode&~2);
  m_weight = KFactor(NULL, lmode);
  msg_Debugging()<<"Weight: "<<m_weight<<"\n";
  if (s_variations->Size()) {
    std::vector<double> &bkw(p_proc->GetMEwgtinfo()->m_bkw);
    if (mode&2) bkw.clear();
    size_t oldsize(bkw.size());
    s_variations->ForEach(
        [this, &lmode](size_t varindex,
                       QCD_Variation_Params& varparams) -> void {
          KFactor(&varparams, lmode);
        });
    msg_Debugging()<<"New K factors: "<<std::vector<double>
      (&bkw[oldsize],&bkw.back()+1)<<"\n";
    for (size_t i(oldsize);i<bkw.size();++i) bkw[i]*=m_weight?1.0/m_weight:0.0;
    msg_Debugging()<<"Weight variations: "<<bkw<<"\n";
  }
  return m_weight;
}

// for applying Higgs effective coupling in MC@NLO style
// provide hf1tt^2+2*hf2tt to multiple a standalone Higgs LO
double HHF2_KFactor::KFactor(QCD_Variation_Params* params, const int& mode)
{
  if (params==NULL) {
    s_as=MODEL::as;
    s_pdf=p_proc->Integrator()->ISR()->PDF(0);
  }
  else {
    s_as=params->p_alphas;
    s_pdf=params->p_pdf1;
  }
  const Vec4D_Vector &p(p_proc->Integrator()->Momenta());
  Scale_Setter_Base *sc(p_proc->ScaleSetter());
  double mur(sqrt(sc->Scale(stp::ren)*(params?params->m_muR2fac:1.0))), Q(p[2].Mass());
  double norm(m_mtmode?ggH1l(Q,Mt,Mb,Mc):1.);
  double as4pi((*s_as)(mur*mur)/(4.0*M_PI));
  if (m_kfmode&1) {
    norm*=sqr(1.0+as4pi*hf1tt(mur,0.,Mt)+sqr(as4pi)*hf2tt(mur,0.,Mt));
  }
  double K(0.0);
  if (!(m_kfmode&1)&&(m_kfmode&2)) {
    K+=sqr(hf1tt(mur,0.,Mt));
    K+=2.0*hf2tt(mur,0.,Mt);
  }
  if (m_kfmode&4) K+=hf2gg(mur,mur); // no difference made using mur instead of muf as log(mu) drops
  K*=sqr(as4pi);
  double weight=K*norm;
  // add two-loop ggH top mass dependence; note that hf0tt(mur,Q,Mt)==ggH1l(Q,Mt,0.,0.)
  if (m_mtmode) weight+=as4pi*2.*(hf1tt(mur,Q,Mt)-hf1tt(mur,0.,Mt))*hf0tt(mur,Q,Mt); 
  if (params) p_proc->GetMEwgtinfo()->m_bkw.push_back(weight);
  return params?1.0:weight;
}

double HNLO_KFactor::KFactor(const int mode) 
{
  DEBUG_FUNC(p_proc->Name()<<" "<<p_proc->Generator()->Name()<<" "<<mode);
  m_weight = KFactor(NULL, mode);
  msg_Debugging()<<"Weight: "<<m_weight<<"\n";
  if (s_variations->Size()) {
    std::vector<double> &bkw(p_proc->GetMEwgtinfo()->m_bkw);
    bkw.clear();
    s_variations->ForEach(
        [this, &mode](size_t varindex,
                       QCD_Variation_Params& varparams) -> void {
          KFactor(&varparams, mode);
        });
    msg_Debugging()<<"New K factors: "<<bkw<<"\n";
    for (size_t i(0);i<bkw.size();++i) bkw[i]*=m_weight?1.0/m_weight:0.0;
    msg_Debugging()<<"Weight variations: "<<bkw<<"\n";
  }
  return m_weight;
}

// NLO K factor
double HNLO_KFactor::KFactor(QCD_Variation_Params* params, const int& mode)
{
  if (params==NULL) {
    s_as=MODEL::as;
    s_pdf=p_proc->Integrator()->ISR()->PDF(0);
  }
  else {
    s_as=params->p_alphas;
    s_pdf=params->p_pdf1;
  }
  const Vec4D_Vector &p(p_proc->Integrator()->Momenta());
  Scale_Setter_Base *sc(p_proc->ScaleSetter());
  double mur(sqrt(sc->Scale(stp::ren)*(params?params->m_muR2fac:1.0))), Q(p[2].Mass());
  double muf(sqrt(sc->Scale(stp::fac)*(params?params->m_muF2fac:1.0)));
  double norm(m_mtmode?ggH1l(Q,Mt,Mb,Mc):1.);
  double H1=2.*hf1tt(mur,0.,Mt);
  double as4pi((*s_as)(mur*mur)/(4.0*M_PI));
  if (m_mtmode && (mode==1 && p_proc->NOut()>1)) { // for B of H+1j
    // add 1j mass dependence here
    Flavour fl0(p_proc->Flavours()[0]);
    Flavour fl1(p_proc->Flavours()[1]);
    int i0(fl0.IsGluon()?0:(long int)fl0);
    int i1(fl1.IsGluon()?0:(long int)fl1);
    Mh=Q;
    double oldnorm(norm);
    if (i0==0&&i1==0) norm*=ggHg1l(2.*(p[0]*p[1]),-2.*(p[0]*p[3]),-2.*(p[1]*p[3]),Mh,Mt,Mb,Mc);
    if (i0!=0&&i1!=0) norm*=qqHg1l(2.*(p[0]*p[1]),-2.*(p[0]*p[3]),-2.*(p[1]*p[3]),Mh,Mt,Mb,Mc);
    if (i0==0&&i1!=0) norm*=qqHg1l(-2.*(p[1]*p[3]),2.*(p[0]*p[1]),-2.*(p[0]*p[3]),Mh,Mt,Mb,Mc);
    if (i0!=0&&i1==0) norm*=qqHg1l(-2.*(p[0]*p[3]),-2.*(p[1]*p[3]),2.*(p[0]*p[1]),Mh,Mt,Mb,Mc);
    if (IsBad(norm)) norm=oldnorm;
    else norm/=ggH1l(Q,Mt,Mb,Mc);
  }
  if (m_kfmode&1) norm*=1.0+as4pi*(H1+4.*sqr(M_PI));
  if (p_proc->NOut()>1) {
    double weight(NLODiffWeight(p_proc,norm,mur*mur,muf*muf,m_k0sq,m_fomode,0,
				params?params->Name():""));
    if (params) p_proc->GetMEwgtinfo()->m_bkw.push_back(weight);
    return params?1.0:weight;
  }
  DEBUG_FUNC(p_proc->Name());
  QT_Selector *jf=(QT_Selector*)
    p_proc->Selector()->GetSelector("NNLOqT_Selector");
  if (jf==NULL) THROW(fatal_error,"Must use selector \"NNLOqT\"");
  Cluster_Amplitude *ampl(NULL);
  if (sc->Amplitudes().size()) ampl=sc->Amplitudes().front();
  if (ampl) ampl->SetNLO(4);
  msg_Debugging()<<"\\mu_F = "<<muf<<", \\mu_R = "<<mur<<"\n";
  int swap(p[0][3]<p[1][3]);
  double x0(swap?p[0].PMinus()/rpa->gen.PBeam(1).PMinus():
	    p[0].PPlus()/rpa->gen.PBeam(0).PPlus());
  double x1(swap?p[1].PPlus()/rpa->gen.PBeam(0).PPlus():
	    p[1].PMinus()/rpa->gen.PBeam(1).PMinus());
  double z0(p_fsmc->ERan("zeta_1")), z1(p_fsmc->ERan("zeta_2"));
  double lnRF=(mur==muf?0.:log(mur/muf));
  double K0(0.),K1(0.);
  K1+=Cgg1gg(x0,x1,z0,z1,jf->QTCut(),muf,Q);
  K1+=Cgg1gq(x0,x1,z0,z1,jf->QTCut(),muf,Q)+Cgg1gq(x1,x0,z1,z0,jf->QTCut(),muf,Q);
  K0=Cgg0gg(x0,x1,z0,z1,jf->QTCut(),muf,Q);
  // based on hf1tt(muR)=hf1tt(muF)
  double weight=1.0+as4pi*(K1/K0+H1+2.*beta0*lnRF);
  if ((m_kfmode&1) || (m_kfmode&2)) weight-=as4pi*(H1+4.*sqr(M_PI));
  weight*=norm;
  // add two-loop ggH top mass dependence; note that hf0tt(mur,Q,Mt)==ggH1l(Q,Mt,0.,0.)
  if (m_mtmode) weight+=as4pi*(2.*hf1tt(mur,Q,Mt)-H1)*hf0tt(mur,Q,Mt);
  if (IsBad(weight)) {
    if (params) p_proc->GetMEwgtinfo()->m_bkw.push_back(norm);
    return params?1.0:norm;
  }
  msg_Debugging()<<"K = "<<weight<<"\n";
  if (params) p_proc->GetMEwgtinfo()->m_bkw.push_back(weight);
  return params?1.0:weight;
}

double HF1_KFactor::KFactor(const int mode) 
{
  DEBUG_FUNC(p_proc->Name()<<" "<<p_proc->Generator()->Name()<<" "<<mode);
  m_weight = KFactor(NULL, mode);
  msg_Debugging()<<"Weight: "<<m_weight<<"\n";
  if (s_variations->Size()) {
    std::vector<double> &bkw(p_proc->GetMEwgtinfo()->m_bkw);
    bkw.clear();
    s_variations->ForEach(
        [this, &mode](size_t varindex,
                       QCD_Variation_Params& varparams) -> void {
          KFactor(&varparams, mode);
        });
    msg_Debugging()<<"New K factors: "<<bkw<<"\n";
    for (size_t i(0);i<bkw.size();++i) bkw[i]*=m_weight?1.0/m_weight:0.0;
    msg_Debugging()<<"Weight variations: "<<bkw<<"\n";
  }
  return m_weight;
}

double HF1_KFactor::KFactor(QCD_Variation_Params* params, const int& mode)
{
  if (params==NULL) {
    s_as=MODEL::as;
    s_pdf=p_proc->Integrator()->ISR()->PDF(0);
  }
  else {
    s_as=params->p_alphas;
    s_pdf=params->p_pdf1;
  }
  const Vec4D_Vector &p(p_proc->Integrator()->Momenta());
  Scale_Setter_Base *sc(p_proc->ScaleSetter());
  double mur(sqrt(sc->Scale(stp::ren)*(params?params->m_muR2fac:1.0))), Q(p[2].Mass());
  double norm(m_mtmode?ggH1l(Q,Mt,Mb,Mc):1.);
  double as4pi((*s_as)(mur*mur)/(4.0*M_PI));
  double K(0.0);
  if (!(m_kfmode&1)&&(m_kfmode&2))
    K+=as4pi*(2.0*hf1tt(mur,0.,Mt)+4.0*sqr(M_PI));
  if (params) p_proc->GetMEwgtinfo()->m_bkw.push_back(K*norm);
  return params?1.0:K*norm;
}

DECLARE_GETTER(HNNLO_KFactor,"HNNLO",
	       KFactor_Setter_Base,KFactor_Setter_Arguments);

KFactor_Setter_Base *ATOOLS::Getter
<KFactor_Setter_Base,KFactor_Setter_Arguments,HNNLO_KFactor>::
operator()(const KFactor_Setter_Arguments &args) const
{
  const Flavour_Vector &fl(args.p_proc->Flavours());
  if (fl[2].Kfcode()==kf_h0) return new HNNLO_KFactor(args);
  return NULL;
}

void ATOOLS::Getter<KFactor_Setter_Base,KFactor_Setter_Arguments,
		    HNNLO_KFactor>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"H NNLO K factor\n";
}

DECLARE_GETTER(HHF1_KFactor,"HHF1",
	       KFactor_Setter_Base,KFactor_Setter_Arguments);

KFactor_Setter_Base *ATOOLS::Getter
<KFactor_Setter_Base,KFactor_Setter_Arguments,HHF1_KFactor>::
operator()(const KFactor_Setter_Arguments &args) const
{
  const Flavour_Vector &fl(args.p_proc->Flavours());
  if (fl[2].Kfcode()==kf_h0) return new HHF1_KFactor(args);
  return NULL;
}

void ATOOLS::Getter<KFactor_Setter_Base,KFactor_Setter_Arguments,
		    HHF1_KFactor>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"H HF1 K factor\n";
}

DECLARE_GETTER(HHF2_KFactor,"HHF2",
	       KFactor_Setter_Base,KFactor_Setter_Arguments);

KFactor_Setter_Base *ATOOLS::Getter
<KFactor_Setter_Base,KFactor_Setter_Arguments,HHF2_KFactor>::
operator()(const KFactor_Setter_Arguments &args) const
{
  const Flavour_Vector &fl(args.p_proc->Flavours());
  if (fl[2].Kfcode()==kf_h0) return new HHF2_KFactor(args);
  return NULL;
}

void ATOOLS::Getter<KFactor_Setter_Base,KFactor_Setter_Arguments,
		    HHF2_KFactor>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"H HF2 K factor\n";
}

DECLARE_GETTER(HNLO_KFactor,"HNLO",
	       KFactor_Setter_Base,KFactor_Setter_Arguments);

KFactor_Setter_Base *ATOOLS::Getter
<KFactor_Setter_Base,KFactor_Setter_Arguments,HNLO_KFactor>::
operator()(const KFactor_Setter_Arguments &args) const
{
  const Flavour_Vector &fl(args.p_proc->Flavours());
  if (fl[2].Kfcode()==kf_h0) return new HNLO_KFactor(args);
  return NULL;
}

void ATOOLS::Getter<KFactor_Setter_Base,KFactor_Setter_Arguments,
		    HNLO_KFactor>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"H NLO K factor\n";
}

DECLARE_GETTER(HF1_KFactor,"HF1",
	       KFactor_Setter_Base,KFactor_Setter_Arguments);

KFactor_Setter_Base *ATOOLS::Getter
<KFactor_Setter_Base,KFactor_Setter_Arguments,HF1_KFactor>::
operator()(const KFactor_Setter_Arguments &args) const
{
  const Flavour_Vector &fl(args.p_proc->Flavours());
  if (fl[2].Kfcode()==kf_h0) return new HF1_KFactor(args);
  return NULL;
}

void ATOOLS::Getter<KFactor_Setter_Base,KFactor_Setter_Arguments,
		    HF1_KFactor>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"H HF1 K factor (NLO)\n";
}
