#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Phys/Variations.H"

namespace PHASIC {

  class DY_KFactor: public KFactor_Setter_Base {
  protected:

    double m_k0sq[2];
    int    m_fomode, m_newfs;

    PHASIC::Multi_Channel *p_fsmc;

  public:

    DY_KFactor(const KFactor_Setter_Arguments &args);

  };// end of class DY_KFactor

  class DYNNLO_KFactor: public DY_KFactor {
  public:

    inline DYNNLO_KFactor(const KFactor_Setter_Arguments &args):
      DY_KFactor(args) {}

    double KFactor(ATOOLS::Variation_Parameters *params,
		   ATOOLS::Variation_Weights *weights,
		   const int &mode);
    double KFactor(const int mode=0);
    double KFactor(const ATOOLS::NLO_subevt &evt);

  };// end of class DYNNLO_KFactor

  class DYNLO_KFactor: public DY_KFactor {
  public:

    inline DYNLO_KFactor(const KFactor_Setter_Arguments &args):
      DY_KFactor(args) {}

    double KFactor(ATOOLS::Variation_Parameters *params,
		   ATOOLS::Variation_Weights *weights,
		   const int &mode);
    double KFactor(const int mode=0);
    double KFactor(const ATOOLS::NLO_subevt &evt);

  };// end of class DYNLO_KFactor

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
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"

#include "QT_Selector.H"
#include "coeffqt.H"
#include "Tools.H"

using namespace SHNNLO;
using namespace PHASIC;
using namespace ATOOLS;

DY_KFactor::DY_KFactor
(const KFactor_Setter_Arguments &args):
  KFactor_Setter_Base(args), p_fsmc(NULL)
{
  Data_Reader read(" ",";","#","=");
  if (s_pdf==NULL) {
    s_as=MODEL::as;
    s_pdf=p_proc->Integrator()->ISR()->PDF(0);
    s_pdfmin[0]=read.GetValue<double>("CSS_PDF_MIN",1.0e-4);
    s_pdfmin[1]=read.GetValue<double>("CSS_PDF_MIN_X",1.0e-2);
  }
  m_k0sq[0]=ToType<double>(rpa->gen.Variable("CSS_FS_PT2MIN"));
  m_k0sq[1]=ToType<double>(rpa->gen.Variable("CSS_IS_PT2MIN"));
  m_fomode=read.GetValue<int>("NNLOqT_FOMODE",0);
  m_newfs=0;
  for (size_t i(0);i<p_proc->NOut();++i)
    if (!Flavour(kf_jet).Includes
	(p_proc->Flavours()[p_proc->NIn()+i])) ++m_newfs;
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

double DYNNLO_KFactor::KFactor(const ATOOLS::NLO_subevt &evt)
{
  return m_weight=1.0;
}

double DYNNLO_KFactor::KFactor(const int mode) 
{
  DEBUG_FUNC(p_proc->Name()<<" "<<p_proc->Generator()->Name());
  const int lmode(mode&~2);
  m_weight=KFactor(NULL,NULL,lmode);
  msg_Debugging()<<"Weight: "<<m_weight<<"\n";
  if (p_proc->VariationWeights()) {
    Variation_Weights vw(p_proc->VariationWeights()->GetVariations());
    std::vector<double> &bkw(p_proc->Caller()->GetMEwgtinfo()->m_bkw);
    if (mode&2) bkw.clear();
    size_t oldsize(bkw.size());
    vw.UpdateOrInitialiseWeights(&DYNNLO_KFactor::KFactor,*this,lmode);
    msg_Debugging()<<"New K factors: "<<std::vector<double>
      (&bkw[oldsize],&bkw.back()+1)<<"\n";
    for (size_t i(oldsize);i<bkw.size();++i) bkw[i]*=m_weight?1.0/m_weight:0.0;
    msg_Debugging()<<"Weight variations: "<<bkw<<"\n";
  }
  return m_weight;
}

double DYNNLO_KFactor::KFactor
(Variation_Parameters *params,Variation_Weights *weights,
 const int &mode) 
{
  if (weights==NULL) {
    s_as=MODEL::as;
    s_pdf=p_proc->Integrator()->ISR()->PDF(0);
  }
  else {
    s_as=params->p_alphas;
    s_pdf=params->p_pdf1;
  }
  Scale_Setter_Base *sc(p_proc->ScaleSetter());
  double muf(sqrt(sc->Scale(stp::fac)*(weights?params->m_muF2fac:1.0)));
  double mur(sqrt(sc->Scale(stp::ren)*(weights?params->m_muR2fac:1.0)));
  if (p_proc->NOut()>(p_proc->Info().Has(nlo_type::real)?m_newfs+1:m_newfs)) {
    double weight(1.0);
    weight=NNLODiffWeight(p_proc,weight,mur*mur,muf*muf,m_k0sq,mode,m_fomode,0,
			  weights?params->m_name:"");
    if (weights) p_proc->Caller()->GetMEwgtinfo()->m_bkw.push_back(weight);
    return weights?1.0:weight;
  }
  if (p_proc->Info().Has(nlo_type::real)) {
    double weight(NNLODeltaWeight(p_proc,1.0,m_fomode));
    if (weights) p_proc->Caller()->GetMEwgtinfo()->m_bkw.push_back(weight);
    return weights?1.0:weight;
  }
  Cluster_Amplitude *ampl(NULL);
  if (sc->Amplitudes().size()) ampl=sc->Amplitudes().front();
  if (ampl) ampl->SetNLO(4);
  if (mode!=1) {
    if (weights) p_proc->Caller()->GetMEwgtinfo()->m_bkw.push_back(1.0);
    return 1.0;
  }
  QT_Selector *jf=(QT_Selector*)
    p_proc->Selector()->GetSelector("NNLOqT_Selector");
  if (jf==NULL) THROW(fatal_error,"Must use selector \"NNLOqT\"");
  const Vec4D_Vector &p(p_proc->Integrator()->Momenta());
  Vec4D hfs;
  for (size_t i(0);i<m_newfs;++i) hfs+=p[p_proc->NIn()+i];
  double Q(hfs.Mass());
  DEBUG_VAR(muf);
  DEBUG_VAR(Q);
  Flavour fl0(p_proc->Flavours()[0]);
  Flavour fl1(p_proc->Flavours()[1]);
  int i0(fl0.IsGluon()?0:(long int)fl0);
  int i1(fl1.IsGluon()?0:(long int)fl1);
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
  DEBUG_VAR(fl0<<" "<<x0);
  DEBUG_VAR(fl1<<" "<<x1);
  // Note that when mur!=muf, need to convert expansion
  // in aS(muf) to aS(mur)
  // C = C0+C1*aS(muf)+C2*aS(muf)^2+C3*aS(muf)^3+...
  //   = C0+C1*aS(mur)
  //    +(C2+2*b0*log(mur/muf)*C1)*aS(mur)^2
  //    +(C3+4*b0*log(mur/muf)*C2+(2*b1*log(mur/muf)+4*b0^2*log(mur/muf)^2)*C1)*aS(mur)^3
  double K=0.;
  K+=Cqq2qiqj(i0,i1,x0,x1,z0,z1,jf->QTCut(),muf,Q)+Cqq2qiqj(i1,i0,x1,x0,z1,z0,jf->QTCut(),muf,Q);
  K+=Cqq2qiqbi(i0,i1,x0,x1,z0,z1,jf->QTCut(),muf,Q)+Cqq2qiqbi(i1,i0,x1,x0,z1,z0,jf->QTCut(),muf,Q);
  K+=Cqq2qg(i0,i1,x0,x1,z0,z1,jf->QTCut(),muf,Q)+Cqq2qg(i1,i0,x1,x0,z1,z0,jf->QTCut(),muf,Q);
  K+=Cqq2gg(i0,i1,x0,x1,z0,z1,jf->QTCut(),muf,Q)+Cqq2qiqi(i0,i1,x0,x1,z0,z1,jf->QTCut(),muf,Q);
  if (mur!=muf) {
    double K2=Cqq1qiqi(i0,i1,x0,x1,z0,z1,jf->QTCut(),muf,Q);
    K2+=Cqq1qg(i0,i1,x0,x1,z0,z1,jf->QTCut(),muf,Q)+Cqq1qg(i1,i0,x1,x0,z0,z1,jf->QTCut(),muf,Q);
    K+=2*beta0*log(mur/muf)*K2;
  }
  K=K/Cqq0qiqi(i0,i1,x0,x1,z0,z1,jf->QTCut(),muf,Q);
  K*=sqr((*s_as)(mur*mur)/(4.0*M_PI));
  DEBUG_VAR(K);
  if (IsBad(K)) {
    if (weights) p_proc->Caller()->GetMEwgtinfo()->m_bkw.push_back(1.0);
    return 1.0;
  }
  DEBUG_VAR(1.0+K);
  if (weights) p_proc->Caller()->GetMEwgtinfo()->m_bkw.push_back(1.0+K);
  return weights?1.0:1.0+K;
}

double DYNLO_KFactor::KFactor(const ATOOLS::NLO_subevt &evt)
{
  return m_weight=1.0;
}

double DYNLO_KFactor::KFactor(const int mode) 
{
  DEBUG_FUNC(p_proc->Name()<<" "<<p_proc->Generator()->Name());
  m_weight=KFactor(NULL,NULL,mode);
  msg_Debugging()<<"Weight: "<<m_weight<<"\n";
  if (p_proc->VariationWeights()) {
    std::vector<double> &bkw(p_proc->Caller()->GetMEwgtinfo()->m_bkw);
    bkw.clear();
    p_proc->VariationWeights()->UpdateOrInitialiseWeights
      (&DYNLO_KFactor::KFactor,*this,mode);
    for (size_t i(0);i<bkw.size();++i) bkw[i]*=m_weight?1.0/m_weight:0.0;
    msg_Debugging()<<"Weight variations: "<<bkw<<"\n";
  }
  return m_weight;
}

double DYNLO_KFactor::KFactor
(Variation_Parameters *params,Variation_Weights *weights,
 const int &mode) 
{
  if (weights==NULL) {
    s_as=MODEL::as;
    s_pdf=p_proc->Integrator()->ISR()->PDF(0);
  }
  else {
    s_as=params->p_alphas;
    s_pdf=params->p_pdf1;
  }
  Scale_Setter_Base *sc(p_proc->ScaleSetter());
  double muf(sqrt(sc->Scale(stp::fac)*(weights?params->m_muF2fac:1.0)));
  double mur(sqrt(sc->Scale(stp::ren)*(weights?params->m_muR2fac:1.0)));
  if (p_proc->NOut()>m_newfs) {
    double weight(1.0);
    weight=NLODiffWeight(p_proc,weight,mur*mur,muf*muf,m_k0sq,m_fomode,0,
			 weights?params->m_name:"");
    if (weights) p_proc->Caller()->GetMEwgtinfo()->m_bkw.push_back(weight);
    return weights?1.0:weight;
  }
  QT_Selector *jf=(QT_Selector*)
    p_proc->Selector()->GetSelector("NNLOqT_Selector");
  if (jf==NULL) THROW(fatal_error,"Must use selector \"NNLOqT\"");
  const Vec4D_Vector &p(p_proc->Integrator()->Momenta());
  Cluster_Amplitude *ampl(NULL);
  if (sc->Amplitudes().size()) ampl=sc->Amplitudes().front();
  if (ampl) ampl->SetNLO(4);
  Vec4D hfs;
  for (size_t i(0);i<m_newfs;++i) hfs+=p[p_proc->NIn()+i];
  double Q(hfs.Mass());
  msg_Debugging()<<"\\mu_F = "<<muf<<", \\mu_R = "<<mur<<"\n";
  Flavour fl0(p_proc->Flavours()[0]);
  Flavour fl1(p_proc->Flavours()[1]);
  int i0(fl0.IsGluon()?0:(long int)fl0);
  int i1(fl1.IsGluon()?0:(long int)fl1);
  int swap(p[0][3]<p[1][3]);
  double x0(swap?p[0].PMinus()/rpa->gen.PBeam(1).PMinus():
	    p[0].PPlus()/rpa->gen.PBeam(0).PPlus());
  double x1(swap?p[1].PPlus()/rpa->gen.PBeam(0).PPlus():
	    p[1].PMinus()/rpa->gen.PBeam(1).PMinus());
  double z0(p_fsmc->ERan("zeta_1")), z1(p_fsmc->ERan("zeta_2"));
  double K=0.;
  K=Cqq1qiqi(i0,i1,x0,x1,z0,z1,jf->QTCut(),muf,Q);
  K+=Cqq1qg(i0,i1,x0,x1,z0,z1,jf->QTCut(),muf,Q)+Cqq1qg(i1,i0,x1,x0,z0,z1,jf->QTCut(),muf,Q);
  K=K/Cqq0qiqi(i0,i1,x0,x1,z0,z1,jf->QTCut(),muf,Q);
  K*=(*s_as)(mur*mur)/(4.0*M_PI);
  msg_Debugging()<<"K = "<<1.0+K<<"\n";
  if (IsBad(K)) K=0.0;
  if (weights) p_proc->Caller()->GetMEwgtinfo()->m_bkw.push_back(1.0+K);
  return weights?1.0:1.0+K;
}

DECLARE_GETTER(DYNNLO_KFactor,"DYNNLO",
	       KFactor_Setter_Base,KFactor_Setter_Arguments);

KFactor_Setter_Base *ATOOLS::Getter
<KFactor_Setter_Base,KFactor_Setter_Arguments,DYNNLO_KFactor>::
operator()(const KFactor_Setter_Arguments &args) const
{
  return new DYNNLO_KFactor(args);
}

void ATOOLS::Getter<KFactor_Setter_Base,KFactor_Setter_Arguments,
		    DYNNLO_KFactor>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"DY NNLO K factor\n";
}

DECLARE_GETTER(DYNLO_KFactor,"DYNLO",
	       KFactor_Setter_Base,KFactor_Setter_Arguments);

KFactor_Setter_Base *ATOOLS::Getter
<KFactor_Setter_Base,KFactor_Setter_Arguments,DYNLO_KFactor>::
operator()(const KFactor_Setter_Arguments &args) const
{
  return new DYNLO_KFactor(args);
}

void ATOOLS::Getter<KFactor_Setter_Base,KFactor_Setter_Arguments,
		    DYNLO_KFactor>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"DY NLO K factor\n";
}
