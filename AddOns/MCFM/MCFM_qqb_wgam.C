#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "AddOns/MCFM/MCFM_Wrapper.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Phys/Color.H"
#include "MODEL/Main/Standard_Model.H"

namespace MCFM { 

  class MCFM_qqb_wgam: public PHASIC::Virtual_ME2_Base {
  private:
    int                       m_pID; 
    bool                      m_anom;
    double                  * p_p, *p_msqv;
    double                    m_sin2_thetaW;
    double                    m_aqed;
    double                    m_unitarization_n;
    double                    m_normcorr;
    double                    m_cplcorr;
    int                       m_count;
     
  public: 
    MCFM_qqb_wgam(int & pID, bool & swapped,const PHASIC::Process_Info& pi,
		  const Flavour_Vector& flavs, bool anom);
    ~MCFM_qqb_wgam();
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
  };

}

extern "C" { 
  void spinoru_(const int & N, double *p,
		Complex * za, Complex * zb);
  void qqb_wgam_v_(double * p, double & n, double * msqv);
}

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MCFM;
using namespace PHASIC;
using namespace ATOOLS;
using namespace MODEL;
using namespace std;

MCFM_qqb_wgam::MCFM_qqb_wgam(int & pID, bool & swapped,
			     const PHASIC::Process_Info& pi,
			     const Flavour_Vector& flavs,bool anom) :
  Virtual_ME2_Base(pi,flavs), m_pID(pID),
  m_anom(anom),
  m_sin2_thetaW(s_model->ScalarConstant(string("sin2_thetaW"))),
  m_aqed(s_model->ScalarFunction(string("alpha_QED"))),
  m_unitarization_n(0),
  m_normcorr(36./qcdcouple_.ason2pi),
  m_cplcorr(pow(4*M_PI*m_aqed/(m_sin2_thetaW*ewcouple_.gwsq),2)
	    *4*M_PI*m_aqed/ewcouple_.esq),
  m_count(0)
{
  rpa->gen.AddCitation
    (1,"The NLO matrix elements have been taken from MCFM \\cite{}.");
   
  p_p      = new double[4*MCFM_NMX];
  p_msqv   = new double[sqr(2*MCFM_NF+1)];
  m_drmode = m_mode=1;
}

MCFM_qqb_wgam::~MCFM_qqb_wgam()
{
  delete [] p_p;
  delete [] p_msqv;
}

void MCFM_qqb_wgam::Calc(const Vec4D_Vector &p)
{
  double corrfactor(m_cplcorr*m_normcorr);
  // All momenta outgoing
  for (int n(0);n<2;++n)           GetMom(p_p,n,-p[n]);
  // reorder particles
  GetMom(p_p,2,p[3]);
  GetMom(p_p,3,p[4]);
  GetMom(p_p,4,p[2]);

  // set up helicity amplitudes
  spinoru_(5,p_p,zprods_.za,zprods_.zb);

  // Set the unitarisation power in the case of anomalous couplings
  if (m_anom) m_unitarization_n=s_model->ScalarConstant(string("UNITARIZATION_N"));

  long int i(m_flavs[0]), j(m_flavs[1]); 
  if (i==21) { i=0; } // correct for gluons
  if (j==21) { j=0; }
  scale_.musq=m_mur2; // set scale
  scale_.scale=sqrt(scale_.musq);
  if (scale_.scale<1.3) {
    int prevcount = m_count;
    m_count++;
    std::cout <<"count:  " << m_count << std::endl;
  }

  i+=MCFM_NF; j+=MCFM_NF; // so the index runs from 0 over the array

  // set residues
  epinv_.epinv=epinv2_.epinv2=0.0;
  qqb_wgam_v_(p_p,m_unitarization_n,p_msqv);
  double res(p_msqv[j*(2*MCFM_NF+1)+i]*corrfactor);

  epinv_.epinv=1.0;
  qqb_wgam_v_(p_p,m_unitarization_n,p_msqv);
  double res1(p_msqv[j*(2*MCFM_NF+1)+i]*corrfactor);

  epinv2_.epinv2=1.0;
  qqb_wgam_v_(p_p,m_unitarization_n,p_msqv);
  double res2(p_msqv[j*(2*MCFM_NF+1)+i]*corrfactor);
     
  m_res.Finite() = res;
  m_res.IR()     = res1-res;
  m_res.IR2()    = res2-res1;
}

double MCFM_qqb_wgam::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(MCFM_qqb_wgam_Getter,"MCFM_qqb_wgam")
Virtual_ME2_Base *MCFM_qqb_wgam_Getter::operator()(const Process_Info &pi) const
{
  DEBUG_FUNC("");
  if (pi.m_loopgenerator!="MCFM")                       return NULL;
  if ((MODEL::s_model->Name()!=std::string("SM")) &&
      (MODEL::s_model->Name()!=std::string("SM+AGC")))  return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo)                return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    Flavour_Vector fl(pi.ExtractFlavours());
    if (!fl[0].Strong() || !fl[1].Strong())             return NULL;
    if (fl.size()!=5)                                   return NULL;
    // check for fully leptonic FS
    if (!(fl[0].IsQuark() && fl[1].IsQuark()))          return NULL;
    // check for outgoing photon
    if (!(fl[4].Kfcode()==22 || fl[2].Kfcode()==22
	  || fl[3].Kfcode()==22))                       return NULL;
    // Check it is not a Z
    if (fl[3].IsLepton() && fl[4]==fl[3].Bar() ||
	fl[2].IsLepton() && fl[3]==fl[2].Bar())         return NULL;
  
    if (s_model->Name()==std::string("SM+AGC")){
      // make sure there are no CP violating anomalous couplings
      if ( (s_model->ScalarConstant(string("kappat_gamma"))!=0)  || 
    	   (s_model->ScalarConstant(string("lambdat_gamma"))!=0) ||
    	   (s_model->ScalarConstant(string("g1_gamma"))!=1)      ||
    	   (s_model->ScalarConstant(string("g4_gamma"))!=0)     ||
    	   (s_model->ScalarConstant(string("g5_gamma"))!=0)){
    	msg_Error()<<"Warning in "<<METHOD<<":"<<std::endl
    		   <<"   Try to set :"<<std::endl
		   <<"kappat_gamma = "
		   <<s_model->ScalarConstant(string("kappat_gamma"))
		   <<", lambdat_gamma = "
		   <<s_model->ScalarConstant(string("lambdat_gamma"))<<", "
		   <<std::endl
		   <<"g1_gamma = "<<s_model->ScalarConstant(string("g1_gamma"))
		   <<", g4_gamma = "<<s_model->ScalarConstant(string("g4_gamma"))
		   <<", "<<std::endl
		   <<"g5_gamma = "
		   <<s_model->ScalarConstant(string("g5_gamma"))<<std::endl
		   <<"Should be : "<<std::endl
		   <<"kappat_gamma = 0, lambdat_gamma = 0, "<<std::endl
		   <<"g1_gamma = 1, g4_gamma = 0, "<<std::endl
		   <<" and g5_gamma = 0"<<std::endl
		   <<"for MCFM."<<std::endl
		   <<"   Will exit the run."<<std::endl;
	THROW(not_implemented,"Incompatible setting with MCFM.");
	return NULL;
      }
    }

    int pID(0);
    bool swapped(false);
    bool anom(false);

    // Check that the model is SM or SM with anomalous gauge couplings    
    if (pi.m_fi.m_ps.size()==2 || pi.m_fi.m_ps.size()==3) {
      ATOOLS::Flavour fl1(pi.m_fi.m_ps[0].m_fl[0]);
      ATOOLS::Flavour fl2(pi.m_fi.m_ps[1].m_fl[0]);
      if (fl[3].IsLepton() && fl[4]==fl[3].Bar()) {
	if (MODEL::s_model->Name()!=std::string("SM")&&
	    MODEL::s_model->Name()!=std::string("SM+AGC")) {
	  msg_Error()<<"Warning in "<<METHOD<<":"<<std::endl
		     <<"   Try to initialise process qqb->Vgamma in MCFM."
		     <<std::endl
		     <<"   Inconsistent setting with Sherpa: "<<std::endl
		     <<"model = "<<MODEL::s_model->Name()
		     <<"(should be 'SM'or 'SM+AGC')."
		     <<std::endl<<"   Will exit the run."<<std::endl;
	  THROW(not_implemented,"Incompatible setting with MCFM.");
	  return NULL;
	}
      }
    
      // W + gamma final state
      if ((fl[2].IsLepton() && fl[3].IsLepton())
	  || (fl[3].IsLepton() && fl[4].IsLepton())) {
	if (MODEL::s_model->Name()!=std::string("SM") &&
	    MODEL::s_model->Name()!=std::string("SM+AGC")) {
	  msg_Error()<<"Warning in "<<METHOD<<":"<<std::endl
		     <<"   Try to initialise process qqb->Vgamma in MCFM."
		     <<std::endl
		     <<"   Inconsistent setting with Sherpa: "<<std::endl
		     <<"model = "<<MODEL::s_model->Name()<<"(should be 'SM')."
		     <<std::endl<<"   Will exit the run."<<std::endl;
	  THROW(not_implemented,"Incompatible setting with MCFM.");
	  return NULL;
	}
	if ((fl[2].IsUptype() && fl[3].IsDowntype()) ||
	    (fl[4].IsDowntype() && fl[3].IsUptype())){
	  if (fl[2].IsUptype()) {
	    zerowidth_.zerowidth=true;
	  }
	  pID = 290; // W+ + gamma
	  if (MODEL::s_model->Name()=="SM") zerowidth_.zerowidth=false;
	  else if (MODEL::s_model->Name()=="SM+AGC") {
	    zerowidth_.zerowidth=false;
	    anom = true;
	  }
	  swapped=true;
	}
	else if (fl[3].IsDowntype() && fl[4].IsUptype()) {
	  pID = 295; // W- + gamma
	  if (MODEL::s_model->Name()=="SM") zerowidth_.zerowidth=false;
	  else if (MODEL::s_model->Name()=="SM+AGC"){
	    zerowidth_.zerowidth=false;
	    anom = true;
	  }
	  swapped=true;
	}
      } 
    }

    // check that the unitarisation powers match
    if ( (MODEL::s_model->Name()=="SM+AGC") && 
	 (MODEL::s_model->ScalarConstant(string("UNITARIZATION_N"))!=
	  MODEL::s_model->ScalarConstant(string("UNITARIZATION_N3"))) )
      {
	msg_Error()<<"Warning in "<<METHOD<<":"<<std::endl
		   <<"   Should have UNITARIZATION_N = UNITARIZATION_N3"
		   <<std::endl
		   <<"   but have: " << std::endl
		   <<"UNITARIZATION_N = "
		   <<MODEL::s_model->ScalarConstant(string("UNITARIZATION_N"))
		   <<" and UNITARIZATION_N3 = "
		   <<MODEL::s_model->ScalarConstant(string("UNITARIZATION_N3"))
		   <<std::endl<<"   Will exit the run."<<std::endl;
	THROW(not_implemented,"Incompatible setting with MCFM.");
	return NULL;
      }
    
    if (pID!=0) {
      if (nproc_.nproc>=0) {
	if (nproc_.nproc!=pID)
	  THROW(not_implemented,
		"Only one process class allowed when using MCFM");
      }
      nproc_.nproc=pID;
      chooser_();
      msg_Info()<<"Initialise MCFM with nproc = "<<nproc_.nproc<<"\n";
      return new MCFM_qqb_wgam(pID,swapped,pi,fl,anom);
    }
  }
  return NULL;
}



