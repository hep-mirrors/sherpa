#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "AddOns/MCFM/MCFM_Wrapper.H"

namespace MCFM {
  class MCFM_qqb_ZZ: public PHASIC::Virtual_ME2_Base {
  private:
    int     m_pID;
    double  m_aqed,m_cplcorr, m_normcorr;
    double *p_p, *p_msqv;

    double CallMCFM(const int & i,const int & j);
  public:
    MCFM_qqb_ZZ(const int & pID,const PHASIC::Process_Info& pi,
		const ATOOLS::Flavour_Vector& flavs);
    ~MCFM_qqb_ZZ();
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
  };

}// end of namespace MCFM

extern "C" { 
  void qqb_zz_v_(double *p,double *msqv); 
}

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MCFM;
using namespace PHASIC;
using namespace ATOOLS;

MCFM_qqb_ZZ::MCFM_qqb_ZZ(const int & pID,const PHASIC::Process_Info& pi,
			 const Flavour_Vector& flavs):
  Virtual_ME2_Base(pi,flavs), m_pID(pID),
  m_aqed(MODEL::s_model->ScalarFunction(std::string("alpha_QED"))),
  m_cplcorr(ATOOLS::sqr(4.*M_PI*m_aqed/ewcouple_.esq)),
  m_normcorr(4.*9./qcdcouple_.ason2pi)
{
  rpa->gen.AddCitation
    (1,"The NLO matrix elements have been taken from MCFM \\cite{Campbell:1999ah}.");
  p_p = new double[4*MCFM_NMX];
  p_msqv = new double[sqr(2*MCFM_NF+1)];
  m_drmode=m_mode=1;
  if (m_pID==82) m_normcorr /= 3.;
}

MCFM_qqb_ZZ::~MCFM_qqb_ZZ()
{
  delete [] p_p;
  delete [] p_msqv;
}

double MCFM_qqb_ZZ::CallMCFM(const int & i,const int & j) {
  qqb_zz_v_(p_p,p_msqv); 
  return p_msqv[mr(i,j)];
}

void MCFM_qqb_ZZ::Calc(const Vec4D_Vector &p)
{
  double corrfactor(m_cplcorr*m_normcorr);

  msg_Out()<<METHOD<<": "
	   <<"m34 = "<<sqrt((p[2]+p[3]).Abs2())<<", "
	   <<"m56 = "<<sqrt((p[4]+p[5]).Abs2())<<".\n";

  // u( q2)+dbar( q1)-->nu(q3)+e^+(q4)+mu^-(q6)+mu^+(q5)

  for (int n(0);n<2;++n) GetMom(p_p,n,p[n]);
  GetMom(p_p,2,p[2]); 
  GetMom(p_p,3,p[4]); 
  GetMom(p_p,4,p[3]); 
  GetMom(p_p,5,p[5]); 

  long int i(m_flavs[0]), j(m_flavs[1]);
  if (i==21) { i=0; }
  if (j==21) { j=0; }
  scale_.musq=m_mur2;
  scale_.scale=sqrt(scale_.musq);

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

double MCFM_qqb_ZZ::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(MCFM_qqb_ZZ_Getter,"MCFM_qqb_ZZ")
Virtual_ME2_Base *MCFM_qqb_ZZ_Getter::operator()(const Process_Info &pi) const
{
  DEBUG_FUNC("");
  if (MODEL::s_model->Name()!=std::string("SM"))        return NULL;
  if (pi.m_loopgenerator!="MCFM")                       return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo)                return NULL;
  if (!pi.m_fi.m_nloqcdtype&nlo_type::loop)             return NULL;
  Flavour_Vector fl(pi.ExtractFlavours());
  // check for right flavour structure:
  // - two incoming quarks + 4 particle FS
  // - assume two propagators - check for alternatives (i.e. 4 lepton FS) later
  // - have W/Z
  if (!(fl[0].IsQuark() && fl[1].IsQuark()))            return NULL;
  if (fl.size()!=6)                                     return NULL;
  msg_Out()<<METHOD<<":";
  for (size_t i=0;i<6;i++) msg_Out()<<" "<<fl[i];
  msg_Out()<<"\n";
  if (!(fl[2].IsLepton() && fl[3].IsLepton() &&
	fl[4].IsLepton() && fl[5].IsLepton())) {
    msg_Error()<<"Error in "<<METHOD<<":\n"
	       <<"   ZZ production for lepton/neutrino final states only.\n";
    THROW(fatal_error,"Not yet working."); 
  }
  int pID(0);
  if (pi.m_fi.m_ps.size()==2 &&
      pi.m_fi.m_ps[0].m_fl[0].Kfcode()==23 && pi.m_fi.m_ps[0].m_fl[1].Kfcode()==23) {
    if (MODEL::s_model->ScalarConstant("Yukawa_b")>0.) {
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   Must switch off Yukawa couplings of b quarks.\n";
      THROW(fatal_error,"Wrong model assupmtions."); 
    }

    if (fl[2].IsDowntype() && fl[4].IsUptype())   pID = 87;
    if (fl[2].IsDowntype() && fl[4].IsDowntype()) pID = 86;
    if (fl[2].IsUptype() && fl[4].IsDowntype()){
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   Order of decays in runcard is: "
		 <<"Z -> nu nub , Z -> l+ l- should be Z -> l+ l- , Z-> nu nub.\n";
      THROW(fatal_error,"Not working."); 
    }
  }
  else if (pi.m_fi.m_ps.size()==4) {
    msg_Error()<<"Error in "<<METHOD<<":\n"
	       <<"   4-lepton final states not yet unfeasible.\n";
    THROW(fatal_error,"Not working."); 
  }  
  


  if (pID!=0) {
    msg_Error()<<"Error in "<<METHOD<<":\n"
	       <<"   this class of processes is not yet working.\n";
    THROW(fatal_error,"Not yet working."); 

    if (nproc_.nproc>=0) {
      if (nproc_.nproc!=pID)
	THROW(not_implemented,
	      "Only one process class allowed when using MCFM");
    }
    nproc_.nproc=pID;
    if (pID==86 || pID==87) zerowidth_.zerowidth=true;
    if (pID==81 || pID==82) zerowidth_.zerowidth=false;
    chooser_();
    msg_Info()<<"Initialise MCFM with nproc = "<<nproc_.nproc<<"\n";
    return new MCFM_qqb_ZZ(pID,pi,fl);
  }
  return NULL;
}

