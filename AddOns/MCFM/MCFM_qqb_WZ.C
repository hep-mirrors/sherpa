#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "AddOns/MCFM/MCFM_Wrapper.H"

namespace MCFM {
  class MCFM_qqb_WZ: public PHASIC::Virtual_ME2_Base {
  private:
    int     m_pID; 
    double  m_aqed,m_cplcorr, m_normcorr;
    double *p_p, *p_msqv;

    double CallMCFM(const int & i,const int & j);
  public:
    MCFM_qqb_WZ(const int & pID, const PHASIC::Process_Info& pi,
		const ATOOLS::Flavour_Vector& flavs);
    ~MCFM_qqb_WZ();
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
  };

}// end of namespace MCFM

extern "C" { 
  void qqb_wz_v_(double *p,double *msqv); 
  void qqb_wz_(double *p,double *a);
}

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MCFM;
using namespace PHASIC;
using namespace ATOOLS;

MCFM_qqb_WZ::MCFM_qqb_WZ(const int & pID, const PHASIC::Process_Info& pi,
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
  if (m_pID==72 || m_pID==77) m_normcorr /= 3.;
}

MCFM_qqb_WZ::~MCFM_qqb_WZ()
{
  delete [] p_p;
  delete [] p_msqv;
}

double MCFM_qqb_WZ::CallMCFM(const int & i,const int & j) {
  qqb_wz_v_(p_p,p_msqv); 
  return p_msqv[mr(i,j)];
}

void MCFM_qqb_WZ::Calc(const Vec4D_Vector &p)
{
  double corrfactor(m_cplcorr*m_normcorr);

  for (int n(0);n<2;++n) GetMom(p_p,n,-p[n]);
  GetMom(p_p,4,p[2]); 
  GetMom(p_p,5,p[3]); 
  GetMom(p_p,2,p[4]); 
  GetMom(p_p,3,p[5]); 

  long int i(m_flavs[0]), j(m_flavs[1]);
  if (i==21) { i=0; } 
  if (j==21) { j=0; }
  scale_.musq=m_mur2; 
  scale_.scale=sqrt(scale_.musq);
  double a[11][11];

  qqb_wz_(p_p,*a);
 
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

double MCFM_qqb_WZ::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(MCFM_qqb_WZ_Getter,"MCFM_qqb_WZ")
Virtual_ME2_Base *MCFM_qqb_WZ_Getter::operator()(const Process_Info &pi) const
{
  DEBUG_FUNC("");
  if (MODEL::s_model->Name()!=std::string("SM")
      && MODEL::s_model->Name()!=std::string("SM+AGC")) return NULL;
  // checking for NLO with MCFM
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
  if ((fl[2].IsUptype() && fl[3].IsDowntype())
      || fl[2].IsDowntype() && fl[3].IsUptype())        return NULL;
  
  bool four_lepton(true);

  if (pi.m_fi.m_ps.size()==2) {
    if (fl[2].IsDowntype()){
      msg_Error()<<"Error in "<<METHOD<<": "
		 <<"   Cannot have Z->l+l- final state due to photon interference.\n";
      THROW(fatal_error,"Not working.");  
    }
    four_lepton=false;
  }

  if (!(fl[2].IsLepton() && fl[3].IsLepton() &&
	fl[4].IsLepton() && fl[5].IsLepton())) {
    msg_Error()<<"Error in "<<METHOD<<":\n"
	       <<"   WZ production for lepton/neutrino final states only.\n";
    THROW(fatal_error,"Not yet working."); 
  }
  if (((fl[2]!=fl[3].Bar()) || abs(fl[4].Kfcode()-fl[5].Kfcode())!=1) 
      && (fl[2]!=fl[4].Bar() || abs(fl[3].Kfcode()-fl[5].Kfcode())!=1)){
    msg_Error()<<"Error in "<<METHOD<<":\n"
	       <<"   Lepton final states unfeasible.\n";
    return NULL;
  } 
  if (pi.m_fi.m_ps.size()==4) { 
    if (MODEL::s_model->Name()==std::string("SM+AGC")) {
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   Must put in explicit boson decays for anomalous couplings.\n";
      THROW(fatal_error,"Not working.");
    }
    if ((fl[2]==fl[4]) || (fl[3]==fl[4]) || (fl[2]==fl[3])) {
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   No interferences in this process allowed.\n";
      THROW(fatal_error,"Not working.");
    }  
    if (fl[2].IsUptype()){
      msg_Error()<<"Error in "<<METHOD<<":\n"
    		 <<"   Cannot decay Z to neutrinos in this way.\n";
      THROW(fatal_error,"Not yet working.");
    }
  } 

  int pID(0); 

  if (fl[2].IsUptype()   && fl[4].IsDowntype()) pID = 77;
  if (fl[2].IsUptype()   && fl[4].IsUptype())   pID = 72;
  if (fl[2].IsDowntype() && fl[4].IsDowntype()) pID = 76;
  if (fl[2].IsDowntype() && fl[4].IsUptype())   pID = 71;
  
  if (pID!=0) {
    if (nproc_.nproc>=0) {
      if (nproc_.nproc!=pID)
	THROW(not_implemented,
	      "Only one process class allowed when using MCFM");
    }
    nproc_.nproc=pID;
    if (four_lepton) zerowidth_.zerowidth=false;
    else             zerowidth_.zerowidth=true;
    chooser_();
    msg_Info()<<"Initialise MCFM with nproc = "<<nproc_.nproc<<"\n";
    return new MCFM_qqb_WZ(pID,pi,fl);
  } 
  return NULL; 
}
 
