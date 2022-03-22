#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "PDF/Main/PDF_Base.H"

namespace PHASIC {

  class DIS_KFactor: public KFactor_Setter_Base {
  protected:

    double m_k0sq[2];
    int    m_fomode;

    PHASIC::Multi_Channel *p_fsmc;

    double F1q(double x, int hel_q);
    double F2q(double x, int hel_q);
    double F3q(double x, int hel_q);

    // 2*Q^2/(x*y^2) times W_munu*L^munu
    double WL(double x,
	      double y,
	      double Q2,
	      int hel_e,
	      int hel_q);

    double m_M2;
    double m_flpre;
    
  public:

    DIS_KFactor(const KFactor_Setter_Arguments &args);

    double KFactor(ATOOLS::QCD_Variation_Params* params,
                   const int mode,
                   const int order);

  };// end of class DIS_KFactor

  class DISNNLO_KFactor: public DIS_KFactor {
  public:

    inline DISNNLO_KFactor(const KFactor_Setter_Arguments &args):
      DIS_KFactor(args) {}

    double KFactor(const int mode=0);
    double KFactor(ATOOLS::QCD_Variation_Params* params, const int& mode);

  };// end of class DISNNLO_KFactor

  class DISNLO_KFactor: public DIS_KFactor {
  public:

    inline DISNLO_KFactor(const KFactor_Setter_Arguments &args):
      DIS_KFactor(args) {}

    double KFactor(const int mode=0);
    double KFactor(ATOOLS::QCD_Variation_Params* params, const int& mode);

  };// end of class DISNLO_KFactor

}// end of namespace PHASIC

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#include "Tools.H"
#include "hard.H"
#include "param.H"
#include "DISinclu.H"

using namespace SHNNLO;
using namespace PHASIC;
using namespace ATOOLS;

DIS_KFactor::DIS_KFactor
(const KFactor_Setter_Arguments &args):
  KFactor_Setter_Base(args), p_fsmc(NULL)
{
  Settings& s = Settings::GetMainSettings();
  if (s_pdf==NULL) {
    s_as=MODEL::as;
    s_pdf=p_proc->Integrator()->ISR()->PDF(0);
    s_pdfmin[0] = s["CSS_PDF_MIN"].Get<double>();
    s_pdfmin[1] = s["CSS_PDF_MIN_X"].Get<double>();
  }
  m_k0sq[0] = s["CSS_FS_PT2MIN"].Get<double>();
  m_k0sq[1] = s["CSS_IS_PT2MIN"].Get<double>();
  m_fomode = s["DISNNLO_FOMODE"].SetDefault(0).Get<int>();
  m_M2=sqr(Flavour(2212).Mass());
  int beam1(p_proc->Flavours()[0].Kfcode());
  if(abs(beam1)!=11) THROW(fatal_error, "Electron beam must be first beam");
  m_flpre=(beam1==11)?1.0:-1.0;
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

double DIS_KFactor::F1q(double x, int hel){
  return F2q(x,hel)/(2.0*x);
}

double DIS_KFactor::F2q(double x, int hel){
  return x;
}

double DIS_KFactor::F3q(double x, int hel){
  if(hel==0) return +1.0;
  if(hel==1) return -1.0;
  THROW(fatal_error, "Internal error");
}

double DIS_KFactor::WL(double x, double y, double Q2, int hel_e, int hel_q){
  double y2=sqr(y); double x2=sqr(x); double sgn=(hel_e==0)?m_flpre:-m_flpre;
  return (1.0-y-x2*y2*m_M2/Q2)*F2q(x,hel_q) + y2*x*F1q(x,hel_q) + sgn*(y-y2/2.0)*x*F3q(x,hel_q); 
}

double DISNNLO_KFactor::KFactor(const int mode) 
{
  DEBUG_FUNC(p_proc->Name()<<" "<<p_proc->Generator()->Name());
  const int lmode(mode&~2);
  if (p_fsmc) {
    s_z[0]=p_fsmc->ERan("zeta_1'");
    s_z[1]=p_fsmc->ERan("zeta_2'");
    s_z[2]=p_fsmc->ERan("zeta_1''");
    s_z[3]=p_fsmc->ERan("zeta_2''");
  }
  m_weight = KFactor(NULL, lmode);
  msg_Debugging()<<"Weight: "<<m_weight<<"\n";
  if (s_variations->Size()) {
    std::vector<double> &bkw(p_proc->Caller()->GetMEwgtinfo()->m_bkw);
    bkw.clear();
    s_variations->ForEach(
        [this, &lmode](size_t varindex,
                       QCD_Variation_Params& varparams) -> void {
          KFactor(&varparams, lmode);
        });
    msg_Debugging()<<"New K factors: "<<bkw<<"\n";
    for (size_t i(0);i<bkw.size();++i) bkw[i]*=m_weight?1.0/m_weight:0.0;
    msg_Debugging()<<"Weight variations: "<<bkw<<"\n";
  }
  if (p_proc->NOut()>2 && m_fomode && rpa->gen.NumberOfTrials()) {
    bool gen=rpa->gen.NumberOfTrials()>s_ntrials;
    s_ntrials=rpa->gen.NumberOfTrials();
    if (gen) s_disc=ran->Get();
    else msg_Debugging()<<"keep random point\n";
    Scale_Setter_Base *sc(p_proc->ScaleSetter());
    if (sc->Amplitudes().size()) {
      Cluster_Amplitude *ampl(sc->Amplitudes().front()->Last());
      if (ampl->Legs().size()>4) return m_weight=0.0;
      m_weight*=2.0;
      if (s_disc>0.5) {
	ampl->SetNLO(256);
	m_weight*=-1.0;
	msg_Debugging()<<"project to Born\n";
      }
    }
  }
  return m_weight;
}

double DISNNLO_KFactor::KFactor(QCD_Variation_Params* params, const int& mode)
{
  if (params==NULL) {
    s_as=s_as;
    s_pdf=p_proc->Integrator()->ISR()->PDF(1);
  }
  else {
    s_as=params->p_alphas;
    s_pdf=params->p_pdf2;
  }
  if (p_proc->NOut()>2) {
    Scale_Setter_Base *sc(p_proc->ScaleSetter());
    double mur2(sc->Scale(stp::ren)*(params?params->m_muR2fac:1.0));
    double muf2(sc->Scale(stp::fac)*(params?params->m_muF2fac:1.0));
    double weight(1.0);
    weight=NNLODiffWeight(p_proc,weight,mur2,muf2,m_k0sq,
			  mode,m_fomode,1,params?params->Name():"");
    if (params) p_proc->Caller()->GetMEwgtinfo()->m_bkw.push_back(weight);
    return params?1.0:weight;
  }
  return DIS_KFactor::KFactor(params, mode, 1);
}

double DISNLO_KFactor::KFactor(const int mode) 
{
  DEBUG_FUNC(p_proc->Name()<<" "<<p_proc->Generator()->Name());
  const int lmode(mode&~2);
  if (p_fsmc) {
    s_z[0]=p_fsmc->ERan("zeta_1'");
    s_z[1]=p_fsmc->ERan("zeta_2'");
    s_z[2]=p_fsmc->ERan("zeta_1''");
    s_z[3]=p_fsmc->ERan("zeta_2''");
  }
  m_weight = KFactor(NULL, lmode);
  msg_Debugging()<<"Weight: "<<m_weight<<"\n";
  if (s_variations->Size()) {
    std::vector<double> &bkw(p_proc->Caller()->GetMEwgtinfo()->m_bkw);
    bkw.clear();
    s_variations->ForEach(
        [this, &lmode](size_t varindex,
                       QCD_Variation_Params& varparams) -> void {
          KFactor(&varparams, lmode);
        });
    for (size_t i(0);i<bkw.size();++i) bkw[i]*=m_weight?1.0/m_weight:0.0;
    msg_Debugging()<<"Weight variations: "<<bkw<<"\n";
  }
  if (p_proc->NOut()>2 && m_fomode && rpa->gen.NumberOfTrials()) {
    Scale_Setter_Base *sc(p_proc->ScaleSetter());
    if (sc->Amplitudes().size()) {
      Cluster_Amplitude *ampl(sc->Amplitudes().front()->Last());
      if (ampl->Legs().size()>4) return m_weight=0.0;
      m_weight*=2.0;
      if (ran->Get()>0.5) {
	ampl->SetNLO(256);
	m_weight*=-1.0;
	msg_Debugging()<<"project to Born\n";
      }
    }
  }
  return m_weight;
}

double DISNLO_KFactor::KFactor(QCD_Variation_Params* params, const int& mode)
{
  if (params==NULL) {
    s_as=s_as;
    s_pdf=p_proc->Integrator()->ISR()->PDF(1);
  }
  else {
    s_as=params->p_alphas;
    s_pdf=params->p_pdf2;
  }
  Scale_Setter_Base *sc(p_proc->ScaleSetter());
  double muf(sqrt(sc->Scale(stp::fac)*(params?params->m_muF2fac:1.0)));
  double mur(sqrt(sc->Scale(stp::ren)*(params?params->m_muR2fac:1.0)));
  if (p_proc->NOut()>2) {
    double weight(1.0);
    weight=NLODiffWeight(p_proc,weight,mur*mur,muf*muf,m_k0sq,m_fomode,1,
			 params?params->Name():"");
    if (params) p_proc->Caller()->GetMEwgtinfo()->m_bkw.push_back(weight);
    return params?1.0:weight;
  }
  return DIS_KFactor::KFactor(params, mode, 0);
}

double DIS_KFactor::KFactor(QCD_Variation_Params* params,
                            const int mode,
                            const int order)
{
  DEBUG_FUNC(p_proc->Name());
  const Vec4D_Vector &moms(p_proc->Integrator()->Momenta());
  Vec4D k(moms[0]),p(moms[1]),kp(moms[2]),pp(moms[3]),P(rpa->gen.PBeam(1));
  Vec4D q=k-kp; double Q2(-q.Abs2());
  double x=Q2/(2.0*q*P); double y=(q*P)/(k*P);
  if(x>1.0){
    msg_Out() << "Reject kinematics, x=" << x << ">1.0" << std::endl;
    return 0.0;
  }
  std::vector<METOOLS::Spin_Amplitudes> amps;
  std::vector<std::vector<Complex> >    cols;
  p_proc->FillAmplitudes(amps,cols);
  METOOLS::Spin_Amplitudes amp=amps[0];

  // 1: lambda=-1
  // 0: lambda=+1
  // for NLO: knlo = alphaS(mur)/4/Pi * DISinclu_1(muf)/DISinclu_0(muf)
  // for NNLO: knnlo = (alphaS(mur)/4/Pi)^2 * { DISinclu_2(muf)/DISinclu_0(muf)+2*beta0*log(mur/muf)*DISinclu_1(muf)/DISinclu_0(muf) };
  // static double beta0 = 11./3.*CA-4./3.*TF*nf;
  
  Scale_Setter_Base *sc(p_proc->ScaleSetter());
  double muf(sqrt(sc->Scale(stp::fac)*(params?params->m_muF2fac:1.0)));
  double mur(sqrt(sc->Scale(stp::ren)*(params?params->m_muR2fac:1.0)));
  double mu2(sc->Scale(stp::ren)*(params?params->m_muR2fac:1.0));
  Cluster_Amplitude *ampl(NULL);
  if (sc->Amplitudes().size()) ampl=sc->Amplitudes().front();
  if (ampl) ampl->SetNLO(4);
  double z(p_fsmc->ERan("zeta_1"));
  double anti=(p_proc->Flavours()[0]).IsAnti()?1.0:-1.0;
  anti*=(p_proc->Flavours()[1]).IsAnti()?1.0:-1.0;
  double asmur((*s_as)(sqr(mur)));
  double nom  =0.0;
  double denom=0.0;
  int    flav=(long int)p_proc->Flavours()[1];

  for(size_t pol_e(0); pol_e<2; pol_e++){
    for(int pol_q(0); pol_q<2; pol_q++){
      std::vector<int> hels;
      hels.push_back(pol_e); hels.push_back(pol_q);
      hels.push_back(pol_e); hels.push_back(pol_q);
      int h=-((pol_q==pol_e)?1:-1)*anti;
      double thishel = std::norm(amp.Get(hels));
      double fac     = 1.0;
      double order_0 = DISinclu_0(flav, h, x, y, z, sqrt(Q2), muf);
      double order_1 = DISinclu_1(flav, h, x, y, z, sqrt(Q2), muf);
      fac           += asmur/(4.0*M_PI)*order_1/order_0;
      if (order) {
	double order_2 = DISinclu_2(flav, h, x, y, z, sqrt(Q2), muf);
	fac += sqr(asmur/(4.0*M_PI))*(order_2/order_0 + 2.0*Beta0(5.0)*log(mur/muf)*order_1/order_0);
      }
      nom           += thishel*fac;
      denom         += thishel;
    }
  }
  double weight=nom/denom;
  if (IsBad(weight)) weight=1.0;
  if (params) p_proc->Caller()->GetMEwgtinfo()->m_bkw.push_back(weight);
  return params?1.0:weight;
}

DECLARE_GETTER(DISNNLO_KFactor,"DISNNLO",
	       KFactor_Setter_Base,KFactor_Setter_Arguments);

KFactor_Setter_Base *ATOOLS::Getter
<KFactor_Setter_Base,KFactor_Setter_Arguments,DISNNLO_KFactor>::
operator()(const KFactor_Setter_Arguments &args) const
{
  return new DISNNLO_KFactor(args);
}

void ATOOLS::Getter<KFactor_Setter_Base,KFactor_Setter_Arguments,
		    DISNNLO_KFactor>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"DIS NNLO K factor\n";
}

DECLARE_GETTER(DISNLO_KFactor,"DISNLO",
	       KFactor_Setter_Base,KFactor_Setter_Arguments);

KFactor_Setter_Base *ATOOLS::Getter
<KFactor_Setter_Base,KFactor_Setter_Arguments,DISNLO_KFactor>::
operator()(const KFactor_Setter_Arguments &args) const
{
  return new DISNLO_KFactor(args);
}

void ATOOLS::Getter<KFactor_Setter_Base,KFactor_Setter_Arguments,
		    DISNLO_KFactor>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"DIS NLO K factor\n";
}
