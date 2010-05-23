#include "AMEGIC++/DipoleSubtraction/Single_Virtual_Correction.H"
#include "AMEGIC++/DipoleSubtraction/DipoleSplitting_Base.H"
#include "AMEGIC++/DipoleSubtraction/Flavour_KernelsA.H"
#include "AMEGIC++/DipoleSubtraction/Massive_Kernels.H"
#include "AMEGIC++/DipoleSubtraction/Single_LOProcess_MHV.H"
#include "AMEGIC++/Phasespace/Phase_Space_Generator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "PDF/Main/ISR_Handler.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"

#include "PHASIC++/Process/Virtual_ME2_Base.H"

using namespace AMEGIC;
using namespace PHASIC;
using namespace MODEL;
using namespace PDF;
using namespace BEAM;
using namespace ATOOLS;
using namespace std;

/*-------------------------------------------------------------------------------

  Constructors

  ------------------------------------------------------------------------------- */

Single_Virtual_Correction::Single_Virtual_Correction() :   
  p_psgen(0), p_dsij(0), p_partner(this), p_LO_process(NULL), 
  p_dipole(0), p_flkern(0), p_masskern(0), p_loopme(0), m_massive(0)
{ 
  m_dalpha = 1.;
  double helpd;
  Data_Reader reader(" ",";","!","=");
  reader.AddComment("#");
  reader.SetInputPath(rpa.GetPath());
  reader.SetInputFile(rpa.gen.Variable("ME_DATA_FILE"));
  if (reader.ReadFromFile(helpd,"DIPOLE_ALPHA")) {
    m_dalpha = helpd;
    msg_Tracking()<<"Set dipole cut alpha="<<m_dalpha<<" . "<<std::endl;
  }
  for (int i=0;i<8;i++) m_kpca[i]=0.;
  for (int i=0;i<8;i++) m_kpcb[i]=0.;
  m_cmur[0]=0.;
  m_cmur[1]=0.;

  rpa.gen.AddCitation(1,"The automated generation of Catani-Seymour Dipole\
 Terms is published under \\cite{Gleisberg:2007md}.");
}



Single_Virtual_Correction::~Single_Virtual_Correction()
{
  if (p_psgen)      {delete p_psgen; p_psgen=0;}
  if (p_dsij) {
    for (size_t i=0;i<p_LO_process->PartonList().size();i++) delete[] p_dsij[i];
    delete[] p_dsij;
  }
  if (p_LO_process) {delete p_LO_process; p_LO_process=0;}
  if (p_dipole)     {delete p_dipole; p_dipole=0;}
  if (p_flkern)     {delete p_flkern; p_flkern=0;}
  if (p_masskern)   {delete p_masskern; p_masskern=0;}
  if (p_loopme)     {delete p_loopme; p_loopme=0;}
}

/*------------------------------------------------------------------------------
  
  Generic stuff for initialization of Single_Virtual_Correctiones
      
  ------------------------------------------------------------------------------*/

void Single_Virtual_Correction::PolarizationNorm() {
  m_Norm = SymmetryFactors() * p_LO_process->GetPolarisation()->Spin_Average(m_nin,&m_flavs.front());
}



/*------------------------------------------------------------------------------

  Initializing libraries, amplitudes, etc.

  ------------------------------------------------------------------------------*/
void Single_Virtual_Correction::SelectLoopProcess()
{
  if (m_pinfo.m_fi.m_nloqcdtype&nlo_type::loop) {
    Process_Info loop_pi(m_pinfo);
    loop_pi.m_fi.m_nloqcdtype=nlo_type::loop;
   
    p_loopme=PHASIC::Virtual_ME2_Base::GetME2(loop_pi);  
  }
  else {
    if (m_pinfo.m_loopgenerator!="Internal")
      p_loopme = PHASIC::Virtual_ME2_Base::GetME2(m_pinfo);
    else p_loopme = new PHASIC::Trivial_Virtual(m_pinfo, Flavour_Vector());
  }

  if (!p_loopme) {
    msg_Error()<<"No Loop ME found for "<<Name()<<endl;
    p_loopme=new PHASIC::Trivial_Virtual(m_pinfo, Flavour_Vector());
  }
}



int Single_Virtual_Correction::InitAmplitude(Model_Base * model,Topology* top,
				      vector<Process_Base *> & links,
				      vector<Process_Base *> & errs)
{
  Init();
  if (!model->CheckFlavours(m_nin,m_nout,&m_flavs.front())) return 0;

  m_pslibname = ToString(m_nin)+"_"+ToString(m_nout);
  m_ptypename = "P"+m_pslibname;
//   m_name+= "_VIRT";

  if (m_pinfo.m_amegicmhv>0) {
    if (CF.MHVCalculable(m_pinfo)) p_LO_process = new Single_LOProcess_MHV(m_pinfo);
    if (m_pinfo.m_amegicmhv==2) return 0;
  }
  if (!p_LO_process) p_LO_process = new Single_LOProcess(m_pinfo);
  p_LO_process->SetTestMoms(p_testmoms);

  p_dipole = new DipoleSplitting_Base();
  p_dipole->SetAlpha(m_dalpha);
  p_flkern = new Flavour_KernelsA();
  p_flkern->SetAlpha(m_dalpha);

  size_t fsgluons(0);
  for (size_t i=m_nin;i<m_flavs.size();i++) {
    if (m_flavs[i].IsMassive()) m_massive=1;
    if (m_flavs[i].IsGluon()) fsgluons++; 
  }
  if (m_massive||fsgluons) {
    p_masskern = new Massive_Kernels();
    if (p_masskern->Nmf()>0) m_massive=1;
    if (!m_massive) {
      delete p_masskern;
      p_masskern=0;
    }
    else {
      m_xpa.resize(p_masskern->Nmf()*fsgluons);
      m_xpb.resize(p_masskern->Nmf()*fsgluons);
    }
    if (p_masskern) p_masskern->SetAlpha(m_dalpha);
  }

  PolarizationNorm();

  if (!p_LO_process->InitAmplitude(model,top,links,errs)) return 0;
  m_iresult = p_LO_process->Result();
  m_oqcd = p_LO_process->OrderQCD()+1;
  m_oew = p_LO_process->OrderEW();
  m_pinfo.m_oqcd=m_oqcd;
  m_pinfo.m_oew=m_oew;

  p_dipole->SetCoupling(p_LO_process->CouplingMap());
  p_flkern->SetCoupling(p_LO_process->CouplingMap());
  if (p_masskern) p_masskern->SetCoupling(p_LO_process->CouplingMap());

  if (p_LO_process!=p_LO_process->Partner()) {
    string partnerID=p_LO_process->Partner()->Name();
    for (size_t j=0;j<links.size();j++) if (Type()==links[j]->Type()) {
      if (partnerID==links[j]->Name()) {
	msg_Tracking()<<"Can map full virtual process: "<<Name()<<" -> "<<partnerID<<" Factor: "<<p_LO_process->GetSFactor()<<endl;
 	p_mapproc = p_partner = (Single_Virtual_Correction*)links[j];
	m_sfactor = p_LO_process->GetSFactor();
 	break;
      }
    }
  }

  if (p_partner==this) {
    SelectLoopProcess();
    p_loopme->SetCouplings((Coupling_Map*)p_LO_process->CouplingMap());

    links.push_back(this);
    p_dsij = new double*[p_LO_process->PartonList().size()];
    for (size_t i=0;i<p_LO_process->PartonList().size();i++) {
      p_dsij[i] = new double[p_LO_process->PartonList().size()];
      for (size_t j=0;j<p_LO_process->PartonList().size();j++) p_dsij[i][j]=0.;
    }
  }
  if (m_pinfo.m_fi.m_nloqcdtype&&nlo_type::vsub) m_wgtinfo.AddMEweights(18);
  else if (m_pinfo.m_fi.m_nloqcdtype&&nlo_type::loop) m_wgtinfo.AddMEweights(2);
  Minimize();
  if (p_partner==this && Result()>0.) SetUpIntegrator();
  return 1;
}

bool AMEGIC::Single_Virtual_Correction::NewLibs() 
{
  bool newloops(0);
  if (p_loopme) newloops=p_loopme->NewLibs();
  return newloops||(p_partner->GetLOProcess()->NewLibs());
}
/*------------------------------------------------------------------------------
  
  Phase space initialization
  
  ------------------------------------------------------------------------------*/

bool AMEGIC::Single_Virtual_Correction::FillIntegrator
(PHASIC::Phase_Space_Handler *const psh)
{
  if (p_partner!=this) return true;
  if (p_LO_process!=p_LO_process->Partner()) return 1;
  if (!SetUpIntegrator()) THROW(fatal_error,"No integrator");
  return Process_Base::FillIntegrator(psh);
}


bool Single_Virtual_Correction::SetUpIntegrator() 
{  
  if (m_nin==2) {
    if ( (m_flavs[0].Mass() != p_int->ISR()->Flav(0).Mass()) ||
	 (m_flavs[1].Mass() != p_int->ISR()->Flav(1).Mass()) ) p_int->ISR()->SetPartonMasses(&m_flavs.front());
    if (CreateChannelLibrary()) return 1;
  }
  if (m_nin==1) if (CreateChannelLibrary()) return 1;
  return 0;
}

bool Single_Virtual_Correction::CreateChannelLibrary()
{
  if (!p_LO_process) return 1;
  p_psgen     = new Phase_Space_Generator(m_nin,m_nout);
  bool newch  = 0;
  if (m_nin>=1)  newch = p_psgen->Construct(p_channellibnames,m_ptypename,m_pslibname,&m_flavs.front(),p_LO_process); 
  if (newch>0) return 0;
  return 1;
}

/*------------------------------------------------------------------------------
  
  Process management
  
  ------------------------------------------------------------------------------*/
void Single_Virtual_Correction::SetLookUp(const bool lookup)
{
  m_lookup=lookup; 
  if (p_LO_process) p_LO_process->SetLookUp(lookup);
}

void Single_Virtual_Correction::Minimize()
{
  if (p_partner==this) return;
  if (p_psgen)      { delete p_psgen; p_psgen=0; }
  if (p_dipole)     { delete p_dipole; p_dipole=0;}
  if (p_flkern)     { delete p_flkern; p_flkern=0;}
  if (p_loopme)     { delete p_loopme; p_loopme=0;}

  m_oqcd      = p_partner->OrderQCD();
  m_oew       = p_partner->OrderEW();
}

/*------------------------------------------------------------------------------

  Calculating total cross sections
  
  ------------------------------------------------------------------------------*/


double Single_Virtual_Correction::Partonic(const ATOOLS::Vec4D_Vector &_moms)
{
  Vec4D_Vector moms(_moms);
  if (m_nin==2 && p_int->ISR() && p_int->ISR()->On()) {
    Poincare cms(moms[0]+moms[1]);
    for (size_t i(0);i<moms.size();++i) cms.Boost(moms[i]);
  }
  return DSigma(moms,m_lookup);
}

double Single_Virtual_Correction::Partonic2()
{ 
  if (p_int->ISR()->On()==0) return 0.;
  return DSigma2(); 
}


double Single_Virtual_Correction::DSigma(const ATOOLS::Vec4D_Vector &_moms,bool lookup)
{
  m_lastxs = 0.;
  if (m_nin==2) {
    for (size_t i=0;i<m_nin+m_nout;i++) {
      if (_moms[i][0]<m_flavs[i].Mass()) return 0.;
    }
  }
  if (m_nin==1) {
    for (size_t i=m_nin;i<m_nin+m_nout;i++) {
      if (_moms[i][0]<m_flavs[i].Mass()) return 0.;
    }
  }
  if (p_partner == this) {
    m_lastdxs = operator()(_moms);
  }
  else {
    if (lookup) {
      m_lastdxs = p_partner->LastDXS()*m_sfactor;
    }
    else m_lastdxs = p_partner->operator()(_moms)*m_sfactor;
  }
  double kpterm = p_partner->Get_KPterms(p_int->ISR()->PDF(0),p_int->ISR()->PDF(1),m_flavs);
  if (p_partner != this) kpterm*=m_sfactor;

  m_wgtinfo.m_w0 = m_lastdxs/m_sfactor;
  p_partner->FillMEwgts(m_wgtinfo); 
  m_wgtinfo*=m_Norm*m_sfactor;

  return m_lastxs = m_Norm * (m_lastdxs+kpterm);
}

double Single_Virtual_Correction::DSigma2() 
{
  if (m_lastxs==0.) return 0.;
  double kpterm = p_partner->Get_KPterms(p_int->ISR()->PDF(1),p_int->ISR()->PDF(0),m_flavs);
  if (p_partner != this) kpterm*=m_sfactor;
  return m_Norm * (m_lastdxs+kpterm);
}


double Single_Virtual_Correction::Calc_Imassive(const ATOOLS::Vec4D *mom) 
{
  double res=0.;
  double mur = p_scale->Scale(stp::ren);
  for (size_t i=0;i<p_LO_process->PartonList().size();i++) {
    for (size_t k=i+1;k<p_LO_process->PartonList().size();k++) {
      int typei = 2*m_flavs[p_LO_process->PartonList()[i]].IntSpin();
      int typek = 2*m_flavs[p_LO_process->PartonList()[k]].IntSpin();
      double sik=2.*mom[p_LO_process->PartonList()[i]]*mom[p_LO_process->PartonList()[k]];
      double mi=m_flavs[p_LO_process->PartonList()[i]].Mass();
      double mk=m_flavs[p_LO_process->PartonList()[k]].Mass();

      p_masskern->Calculate(typei,mur,sik,mi,mk,p_LO_process->PartonList()[i]<m_nin);
      double splf  = p_masskern->I_Fin();
      double splf1 = p_masskern->I_E1();
      double splf2 = p_masskern->I_E2();
      p_masskern->Calculate(typek,mur,sik,mk,mi,p_LO_process->PartonList()[k]<m_nin);
      splf  += p_masskern->I_Fin();
      splf1 += p_masskern->I_E1();
      splf2 += p_masskern->I_E2();

      Vec4D_Vector momv(mom, &mom[m_nin+m_nout]);
      double lsc = log(4.*M_PI*mur/dabs(sik)/p_loopme->Eps_Scheme_Factor(momv));

      splf+=splf1*lsc+splf2*0.5*sqr(lsc);
      res+=p_dsij[i][k]*splf;
    }
  }
  return -res*p_flkern->Coupling();
}

double Single_Virtual_Correction::Calc_I(const ATOOLS::Vec4D *mom) 
{
  if (m_massive) return Calc_Imassive(mom);

  double res=0.;
  double mur = p_scale->Scale(stp::ren);
  for (size_t i=0;i<p_LO_process->PartonList().size();i++) {
    for (size_t k=i+1;k<p_LO_process->PartonList().size();k++) {
      int typei = 2*m_flavs[p_LO_process->PartonList()[i]].IntSpin();
      int typek = 2*m_flavs[p_LO_process->PartonList()[k]].IntSpin();
      Vec4D_Vector momv(mom, &mom[m_nin+m_nout]);
      double lsc = log(4.*M_PI*mur/dabs(2.*mom[p_LO_process->PartonList()[i]]*mom[p_LO_process->PartonList()[k]])/p_loopme->Eps_Scheme_Factor(momv));
      double splf = p_dipole->Vif(typei)+p_dipole->Vif(typek);
      double splf1 = p_dipole->Vie1(typei)+p_dipole->Vie1(typek);
      double splf2 = p_dipole->Vie2(typei)+p_dipole->Vie2(typek);

      splf+=splf1*lsc+splf2*0.5*sqr(lsc);
      res+=p_dsij[i][k]*splf;
      m_cmur[0]+=p_dsij[i][k]*(splf1+splf2*lsc);
      m_cmur[1]+=p_dsij[i][k]*splf2;
    }
  }
  m_cmur[0]*=-p_flkern->Coupling();
  m_cmur[1]*=-p_flkern->Coupling();

  return -res*p_flkern->Coupling();
}

void Single_Virtual_Correction::Calc_KP(const ATOOLS::Vec4D *mom, double x0, double x1, double eta0, double eta1, double weight) 
{
  bool sa=m_flavs[0].Strong();
  bool sb=m_flavs[1].Strong();
  if (!sa && !sb) return;
  if (x0<eta0 || x1<eta1) return; 
  size_t pls=1;
  if (sa&&sb) pls++;
  double muf = p_scale->Scale(stp::fac);
  for (int i=0;i<8;i++) m_kpca[i]=0.;
  for (int i=0;i<8;i++) m_kpcb[i]=0.;

  if (sa) {
    double w=1.-eta0;
    int type=m_flavs[0].IntSpin();
    m_kpca[0]=-w*p_flkern->Kb1(type,x0)+p_flkern->Kb2(type)-p_flkern->Kb4(type,eta0);
    m_kpca[1]=w*(p_flkern->Kb1(type,x0)+p_flkern->Kb3(type,x0));
    m_kpca[2]=-w*p_flkern->Kb1(type+2,x0)+p_flkern->Kb2(type+2)-p_flkern->Kb4(type+2,eta0);
    m_kpca[3]=w*(p_flkern->Kb1(type+2,x0)+p_flkern->Kb3(type+2,x0));
    for (int i=0;i<4;i++) m_kpca[i]*=p_dsij[0][0];

    double t=0.;
    if (!m_massive) {
      for (size_t i=pls;i<p_LO_process->PartonList().size();i++) {
	int itype=m_flavs[p_LO_process->PartonList()[i]].IntSpin();
	t+= p_dsij[0][i]*p_flkern->ft(itype);
      }
      m_kpca[type*2-2]+=t*(-w*p_flkern->t1(x0)+p_flkern->t2()-p_flkern->t4(eta0));
      m_kpca[type*2-1]+=t*w*p_flkern->t1(x0);
    }
    else {
      size_t xpcnt(0);
      for (size_t i=pls;i<p_LO_process->PartonList().size();i++) {
	int spin=m_flavs[p_LO_process->PartonList()[i]].IntSpin();
	double saj=dabs(2.*mom[p_LO_process->PartonList()[0]]*mom[p_LO_process->PartonList()[i]]);
	double muq2=saj;
	if (spin!=2) muq2=sqr(m_flavs[p_LO_process->PartonList()[i]].Mass())/saj;
	m_kpca[0]+=p_dsij[0][i]*
	  (-w*p_masskern->t1(type,spin,muq2,x0)+p_masskern->t2(type,spin,muq2)-p_masskern->t4(type,spin,muq2,eta0));
	m_kpca[1]+=p_dsij[0][i]*
	  (w*(p_masskern->t1(type,spin,muq2,x0)+p_masskern->t3(type,spin,muq2,x0)));
	m_kpca[2]+=p_dsij[0][i]*
	  (-w*p_masskern->t1(type+2,spin,muq2,x0)+p_masskern->t2(type+2,spin,muq2)-p_masskern->t4(type+2,spin,muq2,eta0));
	m_kpca[3]+=p_dsij[0][i]*
	  (w*(p_masskern->t1(type+2,spin,muq2,x0)+p_masskern->t3(type+2,spin,muq2,x0)));	
	if (spin==2) {
	  for (int j=0;j<p_masskern->Nmf();j++) {
	    m_xpa[xpcnt].xp=1.-4.*sqr(p_masskern->FMass(j))/saj;
	    if (m_xpa[xpcnt].xp>eta0) {
	      m_kpca[0]+=p_dsij[0][i]*p_masskern->t6(type,m_xpa[xpcnt].xp);
	      m_kpca[1]+=p_dsij[0][i]*w*p_masskern->t5(type,x0,m_xpa[xpcnt].xp);
	      m_kpca[2]+=p_dsij[0][i]*p_masskern->t6(type+2,m_xpa[xpcnt].xp);
	      m_kpca[3]+=p_dsij[0][i]*w*p_masskern->t5(type+2,x0,m_xpa[xpcnt].xp);
	      
	      m_xpa[xpcnt].kpc=p_dsij[0][i]*
		(-w*p_masskern->t5(type,x0,m_xpa[xpcnt].xp)-p_masskern->t6(type,m_xpa[xpcnt].xp)-p_masskern->t7(type,eta0,m_xpa[xpcnt].xp));
	    }
	    xpcnt++;
	  }
	}
      }
    }

    if (sb) {
      m_kpca[0]-=p_dsij[0][1]*(-w*p_flkern->Kt1(type,x0)+p_flkern->Kt2(type)-p_flkern->Kt4(type,eta0));
      m_kpca[1]-=p_dsij[0][1]*w*(p_flkern->Kt1(type,x0)+p_flkern->Kt3(type,x0));
      m_kpca[2]-=p_dsij[0][1]*(-w*p_flkern->Kt1(type+2,x0)+p_flkern->Kt2(type+2)-p_flkern->Kt4(type+2,eta0));
      m_kpca[3]-=p_dsij[0][1]*w*(p_flkern->Kt1(type+2,x0)+p_flkern->Kt3(type+2,x0));
    }
    
    double asum=0.,fsum=0.;
    for (size_t i=1;i<p_LO_process->PartonList().size();i++) {
      fsum+=p_dsij[0][i];
      asum+=p_dsij[0][i]*log(muf/dabs(2.*mom[p_LO_process->PartonList()[0]]*mom[p_LO_process->PartonList()[i]]));
    }
    asum/=fsum;
    m_kpca[4]=fsum*(-w*p_flkern->P1(type,x0)+p_flkern->P2(type)-p_flkern->P4(type,eta0));
    m_kpca[5]=fsum*w*(p_flkern->P1(type,x0)+p_flkern->P3(type,x0));
    m_kpca[6]=fsum*(-w*p_flkern->P1(type+2,x0)+p_flkern->P2(type+2)-p_flkern->P4(type+2,eta0));
    m_kpca[7]=fsum*w*(p_flkern->P1(type+2,x0)+p_flkern->P3(type+2,x0));
    m_kpca[0]+=asum*m_kpca[4];
    m_kpca[1]+=asum*m_kpca[5];
    m_kpca[2]+=asum*m_kpca[6];
    m_kpca[3]+=asum*m_kpca[7];
  }
  
  if (sb) {
    double w=1.-eta1;
    int type=m_flavs[1].IntSpin();
    m_kpcb[0]=-w*p_flkern->Kb1(type,x1)+p_flkern->Kb2(type)-p_flkern->Kb4(type,eta1);
    m_kpcb[1]=w*(p_flkern->Kb1(type,x1)+p_flkern->Kb3(type,x1));
    m_kpcb[2]=-w*p_flkern->Kb1(type+2,x1)+p_flkern->Kb2(type+2)-p_flkern->Kb4(type+2,eta1);
    m_kpcb[3]=w*(p_flkern->Kb1(type+2,x1)+p_flkern->Kb3(type+2,x1));
    for (int i=0;i<4;i++) m_kpcb[i]*=p_dsij[0][0];

    double t=0.;
    if (!m_massive) {
      for (size_t i=pls;i<p_LO_process->PartonList().size();i++) {
	int itype=m_flavs[p_LO_process->PartonList()[i]].IntSpin();
	t+= p_dsij[pls-1][i]*p_flkern->ft(itype);
      }
      m_kpcb[type*2-2]+=t*(-w*p_flkern->t1(x1)+p_flkern->t2()-p_flkern->t4(eta1));
      m_kpcb[type*2-1]+=t*w*p_flkern->t1(x1);
    }
    else {
      size_t xpcnt(0);
      for (size_t i=pls;i<p_LO_process->PartonList().size();i++) {
	int spin=m_flavs[p_LO_process->PartonList()[i]].IntSpin();
	double saj=dabs(2.*mom[p_LO_process->PartonList()[pls-1]]*mom[p_LO_process->PartonList()[i]]);
	double muq2=saj;
	if (spin!=2) muq2=sqr(m_flavs[p_LO_process->PartonList()[i]].Mass())/saj;
	m_kpcb[0]+=p_dsij[pls-1][i]*
	  (-w*p_masskern->t1(type,spin,muq2,x1)+p_masskern->t2(type,spin,muq2)-p_masskern->t4(type,spin,muq2,eta1));
	m_kpcb[1]+=p_dsij[pls-1][i]*
	  (w*(p_masskern->t1(type,spin,muq2,x1)+p_masskern->t3(type,spin,muq2,x1)));
	m_kpcb[2]+=p_dsij[pls-1][i]*
	  (-w*p_masskern->t1(type+2,spin,muq2,x1)+p_masskern->t2(type+2,spin,muq2)-p_masskern->t4(type+2,spin,muq2,eta1));
	m_kpcb[3]+=p_dsij[pls-1][i]*
	  (w*(p_masskern->t1(type+2,spin,muq2,x1)+p_masskern->t3(type+2,spin,muq2,x1)));	
	if (spin==2) {
	  for (int j=0;j<p_masskern->Nmf();j++) {
	    m_xpb[xpcnt].xp=1.-4.*sqr(p_masskern->FMass(j))/saj;
	    if (m_xpb[xpcnt].xp>eta1) {
	      m_kpcb[0]+=p_dsij[pls-1][i]*p_masskern->t6(type,m_xpb[xpcnt].xp);
	      m_kpcb[1]+=p_dsij[pls-1][i]*w*p_masskern->t5(type,x1,m_xpb[xpcnt].xp);
	      m_kpcb[2]+=p_dsij[pls-1][i]*p_masskern->t6(type+2,m_xpb[xpcnt].xp);
	      m_kpcb[3]+=p_dsij[pls-1][i]*w*p_masskern->t5(type+2,x1,m_xpb[xpcnt].xp);

	      m_xpb[xpcnt].kpc=p_dsij[pls-1][i]*
		(-w*p_masskern->t5(type,x1,m_xpb[xpcnt].xp)-p_masskern->t6(type,m_xpb[xpcnt].xp)-p_masskern->t7(type,eta1,m_xpb[xpcnt].xp));
	    }
	    xpcnt++;
	  }
	}
      }
    }

    if (sa) {
      m_kpcb[0]-=p_dsij[0][1]*(-w*p_flkern->Kt1(type,x1)+p_flkern->Kt2(type)-p_flkern->Kt4(type,eta1));
      m_kpcb[1]-=p_dsij[0][1]*w*(p_flkern->Kt1(type,x1)+p_flkern->Kt3(type,x1));
      m_kpcb[2]-=p_dsij[0][1]*(-w*p_flkern->Kt1(type+2,x1)+p_flkern->Kt2(type+2)-p_flkern->Kt4(type+2,eta1));
      m_kpcb[3]-=p_dsij[0][1]*w*(p_flkern->Kt1(type+2,x1)+p_flkern->Kt3(type+2,x1));
    }

    double asum=0.,fsum=0.;
    for (size_t i=0;i<p_LO_process->PartonList().size();i++) if (i!=pls-1) {
      fsum+=p_dsij[pls-1][i];
      asum+=p_dsij[pls-1][i]*log(muf/dabs(2.*mom[p_LO_process->PartonList()[pls-1]]*mom[p_LO_process->PartonList()[i]]));
    }
    asum/=fsum;
    m_kpcb[4]=fsum*(-w*p_flkern->P1(type,x1)+p_flkern->P2(type)-p_flkern->P4(type,eta1));
    m_kpcb[5]=fsum*w*(p_flkern->P1(type,x1)+p_flkern->P3(type,x1));
    m_kpcb[6]=fsum*(-w*p_flkern->P1(type+2,x1)+p_flkern->P2(type+2)-p_flkern->P4(type+2,eta1));
    m_kpcb[7]=fsum*w*(p_flkern->P1(type+2,x1)+p_flkern->P3(type+2,x1));
    m_kpcb[0]+=asum*m_kpcb[4];
    m_kpcb[1]+=asum*m_kpcb[5];
    m_kpcb[2]+=asum*m_kpcb[6];
    m_kpcb[3]+=asum*m_kpcb[7];
  }

  double gfac=weight;
  if (sa) gfac/=(1.-eta0);
  if (sb) gfac/=(1.-eta1);
  gfac*=p_flkern->Coupling();

  if (sa) {
    for (int i=0;i<8;i++) m_kpca[i]*=gfac;
    for (size_t i=0;i<m_xpa.size();i++) m_xpa[i].kpc*=gfac;
  }
  if (sb) {
    for (int i=0;i<8;i++) m_kpcb[i]*=gfac;
    for (size_t i=0;i<m_xpb.size();i++) m_xpb[i].kpc*=gfac;
  }
}

double Single_Virtual_Correction::Get_KPterms(PDF_Base *pdfa, PDF_Base *pdfb, ATOOLS::Flavour_Vector& flav) 
{
  if ((m_pinfo.m_fi.m_nloqcdtype&nlo_type::vsub)==0) return 0.;
  bool sa=flav[0].Strong();
  bool sb=flav[1].Strong();
  if (!sa && !sb) return 0.;
  double eta0 = 1.; if (sa) eta0 = p_int->ISR()->X1();
  double eta1 = 1.; if (sb) eta1 = p_int->ISR()->X2();
  if (m_x0<eta0 || m_x1<eta1) return 0.; 
  size_t pls=1;
  if (sa&&sb) pls++;
  Flavour gluon(kf_gluon);
  Flavour quark(kf_quark);
  double fa=0.,faq=0.,fag=0.,faqx=0.,fagx=0.;
  double fb=0.,fbq=0.,fbg=0.,fbqx=0.,fbgx=0.;
  double g2massq(0.);
  double muf = p_scale->Scale(stp::fac);

  if (sa) {
    pdfa->Calculate(eta0/m_x0,muf);
    fagx = pdfa->GetXPDF(gluon)/eta0;
    if (flav[0].IsQuark()) faqx = pdfa->GetXPDF(flav[0])/eta0;
    else {
      for (size_t i=0;i<quark.Size();i++) faqx+= pdfa->GetXPDF(quark[i]);
      faqx/=eta0;
    }
    pdfa->Calculate(eta0,muf);
    fa  = pdfa->GetXPDF(flav[0])/eta0;
    if (!fa>0.) return 0.;
    fag = pdfa->GetXPDF(gluon)/eta0;
    if (flav[0].IsQuark()) faq = fa;
    else {
      for (size_t i=0;i<quark.Size();i++) faq+= pdfa->GetXPDF(quark[i]);
      faq/=eta0;
    }

    if (m_massive) {
      for (size_t i=0;i<m_xpa.size();i++) if (m_xpa[i].xp>eta0) {
	pdfa->Calculate(eta0/m_xpa[i].xp,muf);
	g2massq+=m_xpa[i].kpc*pdfa->GetXPDF(m_flavs[0])/eta0/fa;
      }
    }    
  }
  if (sb) {
    pdfb->Calculate(eta1/m_x1,muf);
    fbgx = pdfb->GetXPDF(gluon)/eta1;
    if (flav[1].IsQuark()) fbqx = pdfb->GetXPDF(flav[1])/eta1;
    else {
      for (size_t i=0;i<quark.Size();i++) fbqx+= pdfb->GetXPDF(quark[i]);
      fbqx/=eta1;
    }
    pdfb->Calculate(eta1,muf);
    fb = pdfb->GetXPDF(flav[1])/eta1;
    if (!fb>0.) return 0.;
    fbg = pdfb->GetXPDF(gluon)/eta1;
    if (flav[1].IsQuark()) fbq = fb;
    else {
      for (size_t i=0;i<quark.Size();i++) fbq+= pdfb->GetXPDF(quark[i]);
      fbq/=eta1;
    }

    if (m_massive) {
      for (size_t i=0;i<m_xpb.size();i++) if (m_xpb[i].xp>eta1) {
	pdfb->Calculate(eta1/m_xpb[i].xp,muf);
	g2massq+=m_xpb[i].kpc*pdfb->GetXPDF(m_flavs[1])/eta1/fb;
      }
    }    
  }

  double res=g2massq;
  if (sa) {
    res+= (m_kpca[0]*faq+m_kpca[1]*faqx+m_kpca[2]*fag+m_kpca[3]*fagx)/fa;
  }
  
  if (sb) {
    res+= (m_kpcb[0]*fbq+m_kpcb[1]*fbqx+m_kpcb[2]*fbg+m_kpcb[3]*fbgx)/fb;
  }
  return res * KFactor();
}

void Single_Virtual_Correction::CheckPoleCancelation(const ATOOLS::Vec4D *mom)
{
  cout.precision(15);
    for (size_t i=0;i<m_nin+m_nout;i++) cout<<i<<": "<<mom[i]<<endl;
    cout<<"Process: "<<Name()<<endl;
  cout.precision(6);
  double doublepole=0.;
  double singlepole=0.;
  double mur = p_scale->Scale(stp::ren);
  for (size_t i=0;i<p_LO_process->PartonList().size();i++) {
    for (size_t k=i+1;k<p_LO_process->PartonList().size();k++) {
      int typei = 2*m_flavs[p_LO_process->PartonList()[i]].IntSpin();
      int typek = 2*m_flavs[p_LO_process->PartonList()[k]].IntSpin();
      Vec4D_Vector momv(mom, &mom[m_nin+m_nout]);
      double lsc = log(4.*M_PI*mur/dabs(2.*mom[p_LO_process->PartonList()[i]]*mom[p_LO_process->PartonList()[k]])/p_loopme->Eps_Scheme_Factor(momv));

      doublepole+=p_dsij[i][k]*(p_dipole->Vie2(typei)+p_dipole->Vie2(typek));
      singlepole+=p_dsij[i][k]*(p_dipole->Vie1(typei)+p_dipole->Vie1(typek)+
				(p_dipole->Vie2(typei)+p_dipole->Vie2(typek))*lsc);
    }
  }
  doublepole*=-p_flkern->Coupling();
  singlepole*=-p_flkern->Coupling();
  if (!ATOOLS::IsEqual(doublepole,-p_loopme->ME_E2())) {
    cout<<"Double poles do not cancel: "<<doublepole<<" vs. "<<p_loopme->ME_E2()<<endl;
  }
  if (!ATOOLS::IsEqual(singlepole,-p_loopme->ME_E1())) {
    cout<<"Single poles do not cancel: "<<singlepole<<" vs. "<<p_loopme->ME_E1()<<endl;
  }
}

double Single_Virtual_Correction::operator()(const ATOOLS::Vec4D_Vector &mom)
{
  if (p_partner!=this) {
    p_partner->Integrator()->SetMomenta(p_int->Momenta());
    return p_partner->operator()(mom)*m_sfactor;
  }
  Integrator()->RestoreInOrder();

  double M2(0.);
  m_cmur[0]=0.;
  m_cmur[1]=0.;

  p_LO_process->Calc_AllXS(p_int->Momenta(),&mom.front(),p_dsij);
  if (p_loopme->NeedsBorn()) p_loopme->SetBorn(p_dsij[0][0]);
  p_loopme->SetRenScale(p_scale->Scale(stp::ren));
  p_loopme->Calc(mom);
  double I=0.;
  if (m_pinfo.m_fi.m_nloqcdtype&nlo_type::vsub) I=Calc_I(&mom.front());
  m_x0=1.,m_x1=1.;
  double eta0=1.,eta1=1.;
  double w=1.;
  if (m_flavs[0].Strong()) {
    eta0 = p_int->ISR()->X1();
    m_x0 = eta0+ran.Get()*(1.-eta0);
    w *= (1.-eta0);
//       m_x0 = eta0*std::exp(-ran.Get()*log(eta0));
//       w *= -m_x0*log(eta0);
  }
  if (m_flavs[1].Strong()) {
    eta1 = p_int->ISR()->X2();
    m_x1 = eta1+ran.Get()*(1.-eta1);
    w *= (1.-eta1);
//        m_x1 = eta1*std::exp(-ran.Get()*log(eta1));
//        w *= -m_x1*log(eta1);
  }

  if (m_pinfo.m_fi.m_nloqcdtype&nlo_type::vsub) Calc_KP(&mom.front(),m_x0,m_x1,eta0,eta1,w);

  double lme = 0.;
  if (m_pinfo.m_fi.m_nloqcdtype&nlo_type::loop) {
    lme= p_loopme->ME_Finite();
    m_cmur[0]+=p_loopme->ME_E1()+
      double(m_nin+m_nout-4)*p_dipole->G2()*p_dsij[0][0]*p_flkern->Coupling();
    m_cmur[1]+=p_loopme->ME_E2();
  }
  if (m_pinfo.m_fi.m_nloqcdtype&nlo_type::polecheck) CheckPoleCancelation(&mom.front());
  M2+=I+lme;
  if (m_pinfo.m_fi.m_nloqcdtype&nlo_type::born) M2+=p_dsij[0][0]; 
  if (!(M2>0.) && !(M2<0.)
      && !(M2==0.)) {
    std::cout<<"NAN!!! "<<eta0<<" "<<m_x0<<" "<<eta1<<" "<<m_x1<<" "<<p_dsij[0][0]<<" "<<I<<std::endl;
    for (size_t i=0;i<m_nin+m_nout;i++) cout<<i<<": "<<mom[i]<<endl;
  }

  if (!(M2>0.) && !(M2<0.)
      && !(M2==0.)) {
    std::cout<<"NAN!!! "<<M2<<std::endl;
  }
  return M2 * KFactor();
}



int Single_Virtual_Correction::NumberOfDiagrams() { 
  if (p_partner==this) return p_LO_process->NumberOfDiagrams(); 
  return p_partner->NumberOfDiagrams();
}

Point * Single_Virtual_Correction::Diagram(int i) { 
  if (p_partner==this) return p_LO_process->Diagram(i); 
  return p_partner->Diagram(i);
} 

void Single_Virtual_Correction::FillAmplitudes(METOOLS::Amplitude_Tensor* atensor,double sfactor)
{ }

void Single_Virtual_Correction::AddChannels(std::list<std::string>* tlist) 
{ 
  if (p_partner==this) {    
    list<string>* clist = p_channellibnames;
    for (list<string>::iterator it=clist->begin();it!=clist->end();++it) {
      bool hit = 0;
      for (list<string>::iterator jt=tlist->begin();jt!=tlist->end();++jt) {
	if ((*it)==(*jt)) {
	  hit = 1;
	  break;
	}
      }
      if (!hit) tlist->push_back((*it));
    }
  }
}

void Single_Virtual_Correction::SetScale(const Scale_Setter_Arguments &args)
{
  if (!p_LO_process->IsMapped()) p_LO_process->SetScale(args);
  SetScaleSetter(p_LO_process->Partner()->ScaleSetter());
}

void Single_Virtual_Correction::FillMEwgts(PHASIC::ME_wgtinfo& wgtinfo) {
  wgtinfo.m_y1=m_x0;
  wgtinfo.m_y2=m_x1;
  if (wgtinfo.m_nx<2) return;
  for (int i=0;i<2;i++) wgtinfo.p_wx[i]=m_cmur[i];
  if (wgtinfo.m_nx<18) return;
  for (int i=0;i<4;i++) wgtinfo.p_wx[i+2]=m_kpca[i];
  for (int i=0;i<4;i++) wgtinfo.p_wx[i+6]=m_kpcb[i];
  for (int i=0;i<4;i++) wgtinfo.p_wx[i+10]=m_kpca[i+4];
  for (int i=0;i<4;i++) wgtinfo.p_wx[i+14]=m_kpcb[i+4];
}
