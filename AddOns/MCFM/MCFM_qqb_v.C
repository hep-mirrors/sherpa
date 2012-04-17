#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "AddOns/MCFM/MCFM_Wrapper.H"

namespace MCFM {

  class MCFM_qqb_v: public PHASIC::Virtual_ME2_Base {
  private:
    int      m_pID, m_njets;         // process ID from MCFM (this is taken 
                                     //  from the "Getter" function, hardwired
                                     //  there) and # of jets
    double * p_p, * p_msqv;          // momenta as arrays of doubles for Fortran
                                     // compliance and |M|^2
    MODEL::Running_AlphaS * p_as;    // SHERPA's strong coupling constant
    double   m_normcorr;             // part of the correction to make them
                                     // talk with each other in a meaningful
                                     // way.

    double CallMCFM(const int & i,const int & j); 
    // produce the p_msqv in dependence on the flavour tags i and j
  public:
    MCFM_qqb_v(const int & pID,const PHASIC::Process_Info& pi,
	       const ATOOLS::Flavour_Vector& flavs);
    ~MCFM_qqb_v();
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
  };

}// end of namespace MCFM


// routines to be taken from MCFM, that's why they're declared external! 
extern "C" { 
  void qqb_z_v_(double *p,double *msqv); 
  void qqb_z1jet_v_(double *p,double *msqv); 
  void qqb_z2jet_v_(double *p,double *msqv); 
  void qqb_w_v_(double *p,double *msqv); 
  void qqb_w1jet_v_(double *p,double *msqv); 
  void qqb_w2jet_v_(double *p,double *msqv); 
}

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MCFM;
using namespace PHASIC;
using namespace ATOOLS;

MCFM_qqb_v::MCFM_qqb_v(const int & pID,const Process_Info& pi,
		       const Flavour_Vector& flavs):
  Virtual_ME2_Base(pi,flavs), m_pID(pID), m_njets(0),
  p_as((MODEL::Running_AlphaS *)
       MODEL::s_model->GetScalarFunction(std::string("alpha_S"))),
  m_normcorr(4.*9./qcdcouple_.ason2pi *
	     ((m_pID==32||m_pID==42||m_pID==46)?1./3.:1.))
{
  p_p = new double[4*MCFM_NMX];
  p_msqv = new double[sqr(2*MCFM_NF+1)];
  m_drmode=m_mode=1;
  if (m_pID==11||m_pID==16||m_pID==41||m_pID==42) {
    m_njets = 1;
    rpa->gen.AddCitation
      (1,"The NLO matrix elements have been taken from MCFM.");
  }
  else if (m_pID==22||m_pID==27||m_pID==44||m_pID==46) {
    m_njets = 2;
    rpa->gen.AddCitation
      (1,"The NLO matrix elements have been taken from MCFM \\cite{Campbell:2002tg}.");
  }
  else {
    rpa->gen.AddCitation
      (1,"The NLO matrix elements have been taken from MCFM.");
  }
}

MCFM_qqb_v::~MCFM_qqb_v()
{
  delete [] p_p;
  delete [] p_msqv;
}

double MCFM_qqb_v::CallMCFM(const int & i,const int & j) {
  switch (m_pID) {
  case 31:
  case 32:  qqb_z_v_(p_p,p_msqv); break;
  case 41: 
  case 42:  qqb_z1jet_v_(p_p,p_msqv); break;
  case 44:  
  case 46:  qqb_z2jet_v_(p_p,p_msqv); break;
  case 1:  
  case 6:   qqb_w_v_(p_p,p_msqv); break;
  case 11: 
  case 16:  qqb_w1jet_v_(p_p,p_msqv); break;
  case 22: 
  case 27:  qqb_w2jet_v_(p_p,p_msqv); break;
  }
  return p_msqv[mr(i,j)];
}

void MCFM_qqb_v::Calc(const Vec4D_Vector &p)
{
  scale_.musq=m_mur2;
  scale_.scale=sqrt(scale_.musq);
  double corrfactor(m_normcorr);
  if (m_njets>0) {
    //msg_Out()<<"   QCD factor: "<<pow((*p_as)(m_mur2)/qcdcouple_.as,m_njets)
    //	     <<"."<<std::endl;
    corrfactor *= pow((*p_as)(m_mur2)/qcdcouple_.as,m_njets);
  }
  for (int n(0);n<2;++n)           GetMom(p_p,n,-p[n]);
  for (size_t n(2);n<p.size();++n) GetMom(p_p,n,p[n]);
  long int i(m_flavs[0]), j(m_flavs[1]);
  if (i==21) { i=0; corrfactor *= 8./3.; }
  if (j==21) { j=0; corrfactor *= 8./3.; }


  epinv_.epinv=epinv2_.epinv2=0.0;
  double res(CallMCFM(i,j)  * corrfactor);
  epinv_.epinv=1.0;
  double res1(CallMCFM(i,j) * corrfactor);
  epinv2_.epinv2=1.0;
  double res2(CallMCFM(i,j) * corrfactor);
  m_res.Finite() = res;
  m_res.IR()     = (res1-res);
  m_res.IR2()    = (res2-res1);
}

double MCFM_qqb_v::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(MCFM_qqb_v_Getter,"MCFM_qqb_v")
Virtual_ME2_Base *MCFM_qqb_v_Getter::operator()(const Process_Info &pi) const
{
  msg_Out()<<"\n In "<<METHOD<<" with "<<pi.m_fi.m_ps.size()<<"... .\n";
  DEBUG_FUNC("");
  if (pi.m_loopgenerator!="MCFM")                    return NULL;
  if (MODEL::s_model->Name()!=std::string("SM") &&
      MODEL::s_model->Name()!=std::string("THDM") &&
      MODEL::s_model->Name()!=std::string("SM+EHC")) return NULL;
  if (pi.m_oew>2)                                    return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo)             return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    int pID(0);
    Flavour_Vector fl(pi.ExtractFlavours());
    if (fl[0].Strong() && fl[1].Strong() &&
	fl[2].IsLepton() && fl[3].IsLepton()) {
      if ((MODEL::s_model->Name()==std::string("THDM") ||
	   MODEL::s_model->Name()==std::string("SM+EHC")) &&
	  pi.m_fi.m_ps.size()>0 && 
	  pi.m_fi.m_ps[0].m_fl[0]==ATOOLS::Flavour(kf_h0)) {
	msg_Error()<<"Warning in "<<METHOD<<":\n"
		   <<"   Try to initialise process with intermediate "
		   <<"Higgs boson.\n"
		   <<"   Will return 0 and hope for the best.\n";
	return NULL;
      }
      if ((fl[2]==fl[3].Bar() && 
	   MODEL::s_model->ScalarConstant("Yukawa_b")>0. &&
	   ATOOLS::Flavour(kf_b).IsMassive()==0) ||
	  MODEL::s_model->Name()!=std::string("SM")) {
	msg_Error()<<"Warning in "<<METHOD<<":"<<std::endl
		   <<"   Try to initialise process qqb->llbar in MCFM.\n"
		   <<"   Inconsistent setting with Sherpa: "<<std::endl
		   <<"YUKAWA_B = "<<MODEL::s_model->ScalarConstant("Yukawa_b")
		   <<" (should be 0 for llbar, to play it safe), and "
		   <<"MODEL = "<<MODEL::s_model->Name()<<" (should be 'SM')."
		   <<std::endl<<"   Will exit the run."<<std::endl;
	exit(1);
	return NULL;
      }
      bool flag = false;
      if (fl[2]==fl[3].Bar()) {
	if (fl[2].IsUptype()) flag = true;
	if (fl.size()==4)                                     pID = flag?32:31;
	if (fl.size()==5 && fl[4].Strong())                   pID = flag?42:41;
	if (fl.size()==6 && fl[4].Strong() && fl[5].Strong()) pID = flag?46:44;
      }
      else if ((fl[2].IsUptype() && fl[3].IsDowntype()) ||
	       (fl[3].IsUptype() && fl[2].IsDowntype())) {
	if (fl[2].IsUptype()) flag = true;
	if (fl.size()==4)                                     pID = flag?1:6;
	if (fl.size()==5 && fl[4].Strong())                   pID = flag?11:16;
	if (fl.size()==6 && fl[4].Strong() && fl[5].Strong()) pID = flag?22:27;
      }
    }
    if (pID>0) {
      if (nproc_.nproc>=0) {
	if (nproc_.nproc!=pID)
	  THROW(not_implemented,
		"Only one process class allowed when using MCFM");
      }
      else {
	msg_Info()<<"\n";
	nproc_.nproc=pID;
	chooser_();
      }
      msg_Info()<<"Initialise MCFM with nproc = "<<nproc_.nproc<<".\n";
      return new MCFM_qqb_v(pID,pi,fl);
    }
  }
  return NULL;
}
