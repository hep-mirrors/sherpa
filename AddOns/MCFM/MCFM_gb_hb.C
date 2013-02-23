#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_Fermion_Mass.H"
#include "AddOns/MCFM/MCFM_Wrapper.H"

namespace MCFM {
  // README:
  // For Higgs production, choose model: MODEL = THDM
  // It is important for the Higgs production to have all five flavours 
  // in the initial state, and only the Yukawa coupling of the b must be
  // switched on.
  
  
  class MCFM_gb_hb: public PHASIC::Virtual_ME2_Base {
  private:
    int                     m_pID;
    double                * p_p, *p_msqv;
    MODEL::Running_AlphaS * p_as;
    MODEL::Running_Fermion_Mass * p_mb;
    double                  m_mh2,m_Gh2;
    double                  m_mtau2,m_Yukb,m_Yuktau;
    double                  m_cplcorr,m_normcorr;
    
    
    double CallMCFM(const int & i,const int & j);
  public:
    MCFM_gb_hb(const int & pID,
	       const PHASIC::Process_Info& pi,
	       const ATOOLS::Flavour_Vector& flavs);
    ~MCFM_gb_hb();
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
  };
  
}// end of namespace MCFM

extern "C" { 
  void   qqb_hg_v_(double *p,double *msqv); 
  double massfrun_(const double & mb,const double & scale,
		   const double & asmz,const int & nloop);
}

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MCFM;
using namespace MODEL;
using namespace PHASIC;
using namespace ATOOLS;

MCFM_gb_hb::MCFM_gb_hb(const int & pID,const Process_Info& pi,
		     const Flavour_Vector& flavs) :
  Virtual_ME2_Base(pi,flavs), m_pID(pID),
  p_as((Running_AlphaS *)
       s_model->GetScalarFunction(std::string("alpha_S"))),
  p_mb((Running_Fermion_Mass *)
       s_model->GetScalarFunction(std::string("m")+
				  std::string(Flavour(kf_b).IDName()))),
  m_mh2(sqr(Flavour(kf_h0).Mass())),
  m_Gh2(sqr(Flavour(kf_h0).Width())),
  m_mtau2(sqr(Flavour(kf_tau).Mass())), 
  m_Yukb((*p_mb)(m_mh2)), m_Yuktau(Flavour(kf_tau).Yuk()),
  m_normcorr(4.*24./qcdcouple_.ason2pi)
{
  masses_.mb = Flavour(kf_b).HadMass();
  
  rpa->gen.AddCitation
    (1,"The NLO matrix elements have been taken from MCFM \\cite{Campbell:2002zm}.");
  switch (m_pID) {
  case 141:
    double mix(s_model->ComplexMatrixElement(std::string("Z_R"),0,0).real());
    double vev(s_model->ScalarConstant(std::string("vev")));
    double tanb2(sqr(s_model->ScalarConstant(std::string("tan(beta)"))));
    susycoup_.susycoup = mix*sqrt(1.+tanb2);
    masses_.hmass = sqrt(m_mh2); 
    masses_.hwidth = sqrt(m_Gh2); 
    m_cplcorr  = pow(2.*masses_.wmass/vev,4.)*sqr(ewcouple_.xw/ewcouple_.esq);
    m_cplcorr *= sqr(m_Yukb * m_Yuktau/masses_.mb);
    break;
  }

  p_p = new double[4*MCFM_NMX];
  p_msqv = new double[sqr(2*MCFM_NF+1)];
  m_drmode=m_mode=1;

  msg_Out()<<METHOD<<"(corr="<<m_cplcorr
	   <<" & susy = "<<susycoup_.susycoup<<").\n";
}

MCFM_gb_hb::~MCFM_gb_hb()
{
  delete [] p_p;
  delete [] p_msqv;
}


double MCFM_gb_hb::CallMCFM(const int & i,const int & j) {
  switch (m_pID) {
  case 141: qqb_hg_v_(p_p,p_msqv); break;
  }
  return p_msqv[mr(i,j)];
}

void MCFM_gb_hb::Calc(const Vec4D_Vector &p)
{
  double corrfactor(m_cplcorr*m_normcorr);

  scale_.musq        = m_mur2;
  scale_.scale       = sqrt(scale_.musq);
  msbarmasses_.mb_msbar = m_Yukb;

  double mbeff_MCFM(massfrun_(msbarmasses_.mb_msbar,sqrt(m_mur2),
			      couple_.amz,2));

  double sh((p[2]+p[3]).Abs2());
  double hprop2((sqr(sh-sqr(masses_.hmass))+
		 sqr(masses_.hmass*masses_.hwidth))/
		(sqr(sh-m_mh2)+m_mh2*m_Gh2));
  double hdec((sh-4.*m_mtau2)/(3.*(sh-4.*sqr(masses_.mb))));
  corrfactor *= hprop2 * hdec / sqr(mbeff_MCFM);
  corrfactor *= 4.*M_PI*(*p_as)(m_mur2)/qcdcouple_.gsq;

  for (int n(0);n<2;++n)                GetMom(p_p,n,-p[n]);
  for (size_t n(2);n<p.size();++n) GetMom(p_p,n,p[n]);
  long int i(m_flavs[0]), j(m_flavs[1]);
  if (i==21) i=0;
  if (j==21) j=0;

  epinv_.epinv=epinv2_.epinv2=0.0;
  double res(CallMCFM(i,j)  * corrfactor);
  epinv_.epinv=1.0;
  double res1(CallMCFM(i,j) * corrfactor);
  epinv2_.epinv2=1.0;
  double res2(CallMCFM(i,j) * corrfactor);
  m_res.Finite() = res;
  m_res.IR()     = (res1-res);
  m_res.IR2()    = (res2-res1);

  msg_Out()<<METHOD<<"("<<i<<" "<<j<<", corr="<<corrfactor<<")--> "
	   <<res<<", "<<res1<<".\n";
}

double MCFM_gb_hb::Eps_Scheme_Factor(const Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(MCFM_gb_hb,"MCFM_gb_hb")
Virtual_ME2_Base *ATOOLS::Getter
<Virtual_ME2_Base,Process_Info,MCFM_gb_hb>::
operator()(const Process_Info &pi) const
{
  return NULL;
  DEBUG_FUNC("");
  if (pi.m_loopgenerator!="MCFM")                       return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo)                return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    // check for right model and absence of b Yukawa couplings
    Flavour_Vector fl(pi.ExtractFlavours());
    // two incoming strongly interacting particles.
    if (!fl[0].Strong() || !fl[1].Strong())             return NULL;
    int pID(0);
    if (pi.m_fi.m_ps.size()<1 || pi.m_fi.m_ps.size()>2) return NULL;
    Flavour flh(pi.m_fi.m_ps[0].m_fl[0]);
    // higgs propagator
    if (flh!=Flavour(kf_h0))                            return NULL;
    if (pi.m_fi.m_ps.size()==2 && 
	!pi.m_fi.m_ps[1].m_fl[0].Strong())              return NULL;

    if (Flavour(kf_b).Yuk()<=0. ||
	s_model->Name()!=std::string("THDM") ||
	!Flavour(kf_h0).IsOn()) {
      msg_Error()<<"Warning in "<<METHOD<<":"<<std::endl
		 <<"   Try to initialise process gb->Hb in MCFM.\n"
		 <<"   Inconsistent setting with Sherpa: \n"
		 <<"YUKAWA[5] = "<<Flavour(kf_b).Yuk()
		 <<" (should be >0), "
		 <<"MODEL = "<<s_model->Name()<<"(should be 'THDM', and "
		 <<"ACTIVE[25] = "<<Flavour(kf_h0).IsOn()<<" (should be 1).\n"
		 <<"   Will exit the run."<<std::endl;
      THROW(not_implemented,"Incompatible setting with MCFM.");
      return NULL;
    }
    // tau tau final state
    if (fl.size()==5 && 
	(fl[2]==fl[3].Bar() && fl[2].Kfcode()==15)) {
      if (Flavour(kf_tau).Yuk()<=0.) {
	msg_Error()<<"Error in "<<METHOD<<":\n"
		   <<"   Setup for gb->[h->tau tau] b, but tau Yukawa = 0.\n";
	THROW(not_implemented,"Incompatible setting with MCFM.");
      }
      pID = 141;
    }
    if (pID>0) {
      zerowidth_.zerowidth=true;
      if (nproc_.nproc>=0) {
	if (nproc_.nproc!=pID)
	  THROW(not_implemented,
		"Only one process class allowed when using MCFM");
      }
      nproc_.nproc=pID;
      chooser_();
      msg_Info()<<"Initialise MCFM with nproc = "<<nproc_.nproc<<"\n";
      return new MCFM_gb_hb(pID,pi,fl);
    }
  }
  return NULL;
}
