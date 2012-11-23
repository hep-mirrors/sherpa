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
    bool                      m_anomcoup;
    double                  * p_p, *p_msqv;
    MODEL::Running_AlphaS   * p_as;
    double                    m_sin2_thetaW;
    double                    m_MW;
    double                    m_GW;
    double                    m_Qu;
    double                    m_Qd;
    double                    m_CF;
    double                    m_NC;
    double                    m_aqed;
    double                    m_kappa_gamma;
    double                    m_lambda_gamma;
    double                    m_unitarization_scale;
    double                    m_unitarization_n;
     
  public:
    MCFM_qqb_wgam(int & pID, bool & swapped,const PHASIC::Process_Info& pi,
		  const Flavour_Vector& flavs, bool anomcoup);
    ~MCFM_qqb_wgam();
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
    Complex agamtree(int p1, int p2, int p3, int p4, int p5, Complex * za,
		     Complex * zb, int hgamma);
    Complex agamvirt(int p1, int p2, int p3, int p4, int p5, Complex * za,
		     Complex * zb, int hgamma);
    void DoMCFM();
  };

}

extern "C" { 
  Complex vpole_(double * s);
  void spinoru_(const int & N, double *p,
		Complex * za, Complex * zb);
  Complex fagamma_(const int &i1, const int &i2,
		const int &i3, const int &i4,
		const int &i5, Complex * za,
		Complex * zb);
  Complex fbgamma_(const int & i1, const int & i2,
		const int & i3, const int & i4,
		const int & i5, Complex * za,
		Complex * zb);
  void qqb_wgam_v_(double * p, double * msqv);
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
			     const Flavour_Vector& flavs,
			     bool anomcoup) :
  Virtual_ME2_Base(pi,flavs), m_pID(pID),
  m_anomcoup(anomcoup),
  p_as((Running_AlphaS *)s_model->GetScalarFunction(string("alpha_S"))),
  m_sin2_thetaW(s_model->ScalarConstant(string("sin2_thetaW"))),
  m_MW(Flavour(kf_Wplus).Mass()),
  m_GW(Flavour(kf_Wplus).Width()),
  m_Qu(Flavour(kf_u).Charge()),
  m_Qd(Flavour(kf_d).Charge()),
  m_CF(CF),
  m_NC(NC),
  m_aqed(s_model->ScalarFunction(string("alpha_QED"))),
  m_kappa_gamma(1),
  m_lambda_gamma(0),
  m_unitarization_scale(0),
  m_unitarization_n(0)
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

Complex MCFM_qqb_wgam::agamtree(int p1, int p2, int p3, int p4, int p5, 
				Complex * za, Complex * zb, 
				int hgamma)
{
  double xfac;
  Complex prp12,prp34,Agamtree;

  int p1mx = (p1-1)*MCFM_NMX-1; // convert indices into usable form.
  int p2mx = (p2-1)*MCFM_NMX-1; // Converting from fortran array to
  int p3mx = (p3-1)*MCFM_NMX-1; // c++ 1d array.
  int p4mx = (p4-1)*MCFM_NMX-1;
  int p5mx = (p5-1)*MCFM_NMX-1;

  prp34=sprods_.s[p3+p4mx]/Complex((sprods_.s[p3+p4mx]-pow(m_MW,2)),m_MW*m_GW);
  if (m_anomcoup){
    //std::cout << "anomcoup" << std::endl;
    m_kappa_gamma         =s_model->ScalarConstant(string("kappa_gamma"));
    m_lambda_gamma        =s_model->ScalarConstant(string("lambda_gamma"));
    m_unitarization_scale =s_model->ScalarConstant(string("UNITARIZATION_SCALE"))/1000;
    m_unitarization_n     =s_model->ScalarConstant(string("UNITARIZATION_N"));
    anomcoup_.tevscale    = m_unitarization_scale;
    anomcoup_.delk_g      = m_kappa_gamma-1;
    anomcoup_.lambda_g    = m_lambda_gamma;
  }
  //std::cout << "tevscale:  " << m_unitarization_scale << std::endl;
  //std::cout << "delk_g:  " << m_kappa_gamma-1 << std::endl;
  //std::cout << "lambda_g:  " << m_lambda_gamma << std::endl;
  xfac=1.0/pow(1.0+sprods_.s[p1+p2mx]/
  		       pow(anomcoup_.tevscale*1000.,2),m_unitarization_n);
  //std::cout << "unitarization MCFM:  " << xfac << std::endl;
  anomcoup_.xdelk_g=xfac*anomcoup_.delk_g;
  anomcoup_.xlambda_g=xfac*anomcoup_.lambda_g;
  
  if (zerowidth_.zerowidth){
    // zerowidth: no final state radiation, so we can set prp12 to zero
    prp12=Complex(0.,0.);}
  else{
    // otherwise, usual Breit-Wigner form
    prp12=sprods_.s[p1+p2mx]/Complex((sprods_.s[p1+p2mx]-pow(m_MW,2)),m_MW*m_GW);
  }
  //  c.f. Eqs.(4.4),(4.5) of hep-ph/9803250 (multiplied by -i)
  //       for the terms proportional to prp34
 
  if (hgamma == -1){  
     Agamtree=-pow(zb[p2+p4mx],2)
      /(sprods_.s[p1+p2mx]-sprods_.s[p3+p4mx]) 
      *(m_Qu*(za[p2+p5mx]/(zb[p4+p3mx]*zb[p1+p5mx])*prp34
	      -za[p4+p5mx]/(zb[p1+p2mx]*zb[p3+p5mx])*prp12)
      +m_Qd*(za[p1+p5mx]/(zb[p4+p3mx]*zb[p2+p5mx])*prp34
	     +za[p4+p5mx]/(zb[p1+p2mx]*zb[p3+p5mx])*prp12)); 
   }
  else if (hgamma == +1){
     Agamtree=-pow(za[p1+p3mx],2)
      /(sprods_.s[p1+p2mx]-sprods_.s[p3+p4mx]) 
      *(m_Qd*(zb[p1+p5mx]/(za[p3+p4mx]*za[p2+p5mx])*prp34
	      +zb[p4+p5mx]/(za[p2+p1mx]*za[p3+p5mx])*prp12)
	+m_Qu*(zb[p2+p5mx]/(za[p3+p4mx]*za[p1+p5mx])*prp34
	       -zb[p4+p5mx]/(za[p2+p1mx]*za[p3+p5mx])*prp12));  
  }
  // additional anomalous coupling component, Eqs. (7)-(9) of hep-ph/0002138
  if  (hgamma == -1) {
    Agamtree=Agamtree+prp34*(m_Qu-m_Qd)*za[p5+p3mx]/
      (2.0*sprods_.s[p3+p4mx]*(sprods_.s[p1+p2mx]-sprods_.s[p3+p4mx])
       *za[p4+p3mx])*((anomcoup_.xdelk_g+anomcoup_.xlambda_g)*zb[p2+p4mx]
    		      *za[p1+p5mx]*za[p4+p3mx]
		      +anomcoup_.xlambda_g*zb[p2+p5mx]*za[p5+p1mx]*za[p3+p5mx]);
  }
  else if (hgamma == +1) {
    Agamtree=Agamtree+prp34*(m_Qd-m_Qu)*zb[p5+p4mx]/
      (2.0*sprods_.s[p3+p4mx]*(sprods_.s[p1+p2mx]-sprods_.s[p3+p4mx])
       *zb[p3+p4mx])*((anomcoup_.xdelk_g+anomcoup_.xlambda_g)*za[p1+p3mx]
    		      *zb[p2+p5mx]*zb[p3+p4mx]
    		      +anomcoup_.xlambda_g*za[p1+p5mx]*zb[p5+p2mx]*zb[p4+p5mx]);
  }
  return Agamtree;
}

Complex MCFM_qqb_wgam::agamvirt(int p1, int p2, int p3, int p4, int p5, 
				Complex * za, Complex * zb, int hgamma)
{
  Complex prop,vpl,Agamvirt;
  double s34,s12;
  s34=real(za[(p3-1) + (p4-1)*MCFM_NMX]*zb[(p4-1) + (p3-1)*MCFM_NMX]);
  s12=real(za[(p1-1) + (p2-1)*MCFM_NMX]*zb[(p2-1) + (p1-1)*MCFM_NMX]);
  prop=s34/Complex(s34-pow(m_MW,2),m_MW * m_GW);
  vpl=vpole_(&s12); // set poles

  if    (hgamma == +1) {
    Agamvirt=prop*(m_Qd*fagamma_(p1,p2,p3,p4,p5,za,zb)
		   +m_Qu*fbgamma_(p1,p2,p3,p4,p5,za,zb))
      +vpl*agamtree(p1,p2,p3,p4,p5,za,zb,+1);
  }
  else if (hgamma == -1) {
    Agamvirt=prop*(m_Qu*fagamma_(p2,p1,p4,p3,p5,zb,za)
		  + m_Qd*fbgamma_(p2,p1,p4,p3,p5,zb,za))
      +vpl*agamtree(p1,p2,p3,p4,p5,za,zb,-1);
  }
  std::cout << "Born MCFM:  " << 
    (abs(pow(agamtree(p1,p2,p3,p4,p5,za,zb,+1),2))
     + abs(pow(agamtree(p1,p2,p3,p4,p5,za,zb,-1),2)))*
    2.*3.*pow(4.*M_PI*m_aqed/m_sin2_thetaW,2)*4.*M_PI*m_aqed
	    << std::endl;
   
  return Agamvirt;
} 

void MCFM_qqb_wgam::DoMCFM(){

  double fac;
  double qbq;
  double qqb;
  // factor in the colour factors plus the symmetrisation
  spinoru_(5,p_p,zprods_.za,zprods_.zb); // set inner products
                                         // za and zb.
  fac=m_CF*2.*m_NC*pow(4*M_PI*m_aqed/m_sin2_thetaW,2)*4*M_PI*m_aqed;
     
  if (nwz_.nwz == -1){
    // ie d-ub
    qqb=fac*2.0*real(conj(agamtree(2,1,3,4,5,zprods_.za,zprods_.zb,+1))
    		     *agamvirt(2,1,3,4,5,zprods_.za,zprods_.zb,+1))
      +fac*2.0*real(conj(agamtree(2,1,3,4,5,zprods_.za,zprods_.zb,-1))
    		    *agamvirt(2,1,3,4,5,zprods_.za,zprods_.zb,-1));
  }

  else if (nwz_.nwz == +1){ 
    // ie u-db
    qqb=fac*2.0*real(conj(agamtree(1,2,4,3,5,zprods_.zb,zprods_.za,+1))
		     *agamvirt(1,2,4,3,5,zprods_.zb,zprods_.za,+1))
      +fac*2.0*real(conj(agamtree(1,2,4,3,5,zprods_.zb,zprods_.za,-1))
		    *agamvirt(1,2,4,3,5,zprods_.zb,zprods_.za,-1));
  }
  int k,l;

  double Vsq[2*MCFM_NF+1][2*MCFM_NF+1]={ 0 }; // set matrix with elements 
                                              // of CKM matrix squared.
  Vsq[MCFM_NF+2][MCFM_NF-1]=pow(cabib_.Vud,2);
  Vsq[MCFM_NF+2][MCFM_NF-3]=pow(cabib_.Vus,2);
  Vsq[MCFM_NF+2][MCFM_NF-5]=pow(cabib_.Vub,2);
  Vsq[MCFM_NF+4][MCFM_NF-1]=pow(cabib_.Vcd,2);
  Vsq[MCFM_NF+4][MCFM_NF-3]=pow(cabib_.Vcs,2);
  Vsq[MCFM_NF+4][MCFM_NF-5]=pow(cabib_.Vcb,2);
  
  Vsq[MCFM_NF+1][MCFM_NF-2]=pow(cabib_.Vud,2);
  Vsq[MCFM_NF+3][MCFM_NF-2]=pow(cabib_.Vus,2);
  Vsq[MCFM_NF+5][MCFM_NF-2]=pow(cabib_.Vub,2);
  Vsq[MCFM_NF+1][MCFM_NF-4]=pow(cabib_.Vcd,2);
  Vsq[MCFM_NF+3][MCFM_NF-4]=pow(cabib_.Vcs,2);
  Vsq[MCFM_NF+5][MCFM_NF-4]=pow(cabib_.Vcb,2);
 
  for  (k=0;k<=2*MCFM_NF;k++){
    for (l=0;l<=2*MCFM_NF;l++){
      // set msqv=0 to initalize
      p_msqv[k*(2*MCFM_NF+1) + l]=0.;
      p_msqv[(k*(2*MCFM_NF+1))+l]=Vsq[k][l]*qqb;
    }
  }
}

void MCFM_qqb_wgam::Calc(const Vec4D_Vector &p)
{
  for (int n(0);n<2;++n)           GetMom(p_p,n,-p[n]);
  GetMom(p_p,2,p[3]);
  GetMom(p_p,3,p[4]);
  GetMom(p_p,4,p[2]);
  // All momenta outgoing

  long int i(m_flavs[0]), j(m_flavs[1]); // incoming partons
  if (i==21) { i=0; } // correct for gluons
  if (j==21) { j=0; }
  scale_.musq=m_mur2; // set scale
  scale_.scale=sqrt(scale_.musq);
  double corr = 36./qcdcouple_.ason2pi*
    pow(4*M_PI*m_aqed/(m_sin2_thetaW*ewcouple_.gwsq),2)*4*M_PI*m_aqed/ewcouple_.esq;
   //double corr = 1.;
  double a[11][11];
  i+=MCFM_NF; j+=MCFM_NF; // so the index runs from 0 over the array

  epinv_.epinv=epinv2_.epinv2=0.0;// set single and double poles to zero
  DoMCFM();
  qqb_wgam_v_(p_p,*a);
  double resMCFM(a[j][i]*corr);
  double res(p_msqv[(i*(2*MCFM_NF+1))+j]);
  //  std::cout << "res MCFM  " << resMCFM << std::endl;
  //  std::cout << "res       " << res     << std::endl;

  epinv_.epinv=1.0; // set single poles to one 
  DoMCFM();
  qqb_wgam_v_(p_p,*a);
  double res1MCFM(a[j][i]*corr);
  double res1(p_msqv[(i*(2*MCFM_NF+1))+j]);
  //  std::cout << "res1 MCFM " << res1MCFM << std::endl;
  //  std::cout << "res1      " << res1     << std::endl;

  epinv2_.epinv2=1.0; // set double poles to one
  DoMCFM(); 
  qqb_wgam_v_(p_p,*a);
  double res2MCFM(a[j][i]*corr);
  double res2(p_msqv[(i*(2*MCFM_NF+1))+j]);
  //  std::cout << "res2 MCFM " << res2MCFM << std::endl;
  //  std::cout << "res2      " << res2     << std::endl;

  
   // std::cout << "From MCFM: " << a[j][i]*36.
  //  *pow(4*M_PI*m_aqed/m_sin2_thetaW/ewcouple_.gwsq,2)
  //  *4.*M_PI*m_aqed/ewcouple_.esq << std::endl;
    
  // get poles
  m_res.Finite() = res;
  m_res.IR()     = (res1-res);
  m_res.IR2()    = (res2-res1);
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
      if ( (s_model->ScalarConstant(string("kappat_gamma"))!=0)  || 
	   (s_model->ScalarConstant(string("lambdat_gamma"))!=0) ||
	   (s_model->ScalarConstant(string("g1_gamma"))!=1)      ||
	   (s_model->ScalarConstant(string("g4_gamma"))!=0)     ||
	   (s_model->ScalarConstant(string("g5_gamma"))!=0)){
	msg_Error()<<"Warning in "<<METHOD<<":"<<std::endl
		   <<"   Try to set :"<<std::endl
		   <<"kappat_gamma = "<<s_model->ScalarConstant(string("kappat_gamma"))
		   <<", lambdat_gamma = "<<s_model->ScalarConstant(string("lambdat_gamma"))<<", "
		   <<std::endl
		   <<"g1_gamma = "<<s_model->ScalarConstant(string("g1_gamma"))
		   <<", g4_gamma = "<<s_model->ScalarConstant(string("g4_gamma"))<<", "<<std::endl
		   <<"g5_gamma = "<<s_model->ScalarConstant(string("g5_gamma"))<<std::endl
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
    bool anomcoup(false);
    
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
		     <<"model = "<<MODEL::s_model->Name()<<
	    "(should be 'SM'or 'SM+AGC')."
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
	    anomcoup = true;
	  }
	  swapped=true;
	}
	else if (fl[3].IsDowntype() && fl[4].IsUptype()) {
	  pID = 295; // W- + gamma
	  if (MODEL::s_model->Name()=="SM") zerowidth_.zerowidth=false;
	  else if (MODEL::s_model->Name()=="SM+AGC"){
	    zerowidth_.zerowidth=false;
	    anomcoup = true;
	  }
	  swapped=true;
	}
      } 
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
      return new MCFM_qqb_wgam(pID,swapped,pi,fl,anomcoup);
    }
  }
  return NULL;
}



