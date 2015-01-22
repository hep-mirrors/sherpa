#include "AMEGIC++/DipoleSubtraction/Single_Virtual_Correction.H"
#include "AMEGIC++/DipoleSubtraction/DipoleSplitting_Base.H"
#include "AMEGIC++/DipoleSubtraction/Single_LOProcess_MHV.H"
#include "AMEGIC++/DipoleSubtraction/Single_LOProcess_External.H"
#include "AMEGIC++/Phasespace/Phase_Space_Generator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "PDF/Main/ISR_Handler.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"

#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/My_MPI.H"

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
  p_psgen(0), p_partner(this), p_LO_process(NULL), 
  p_dipole(0), p_kpterms(0), p_loopme(0),
  m_bsum(0.0), m_vsum(0.0), m_isum(0.0), m_n(0.0),
  m_mbsum(0.0), m_mvsum(0.0), m_misum(0.0), m_mn(0.0)
{ 
  m_dalpha = 1.;
  m_kpcemode = 0;
  int helpi;
  double helpd;
  Data_Reader reader(" ",";","!","=");
  reader.AddComment("#");
  reader.SetInputPath(rpa->GetPath());
  reader.SetInputFile(rpa->gen.Variable("ME_DATA_FILE"));
  if (reader.ReadFromFile(helpi,"KP_CHECK_ENERGY")) {
    m_kpcemode = helpi;
    msg_Tracking()<<"Set KP-term energy check mode "<<m_kpcemode<<" . "<<std::endl;
  }
  m_checkborn = false;
  if (reader.ReadFromFile(helpi,"CHECK_BORN")) {
    m_checkborn = helpi;
    msg_Tracking()<<"Set Born check mode "<<m_checkborn<<" . "<<std::endl;
  }
  m_checkpoles = false;
  if (reader.ReadFromFile(helpi,"CHECK_POLES")) {
    m_checkpoles = helpi;
    msg_Tracking()<<"Set infrared poles check mode "<<m_checkpoles<<" . "<<std::endl;
  }
  m_checkpolesthreshold = 0.;
  if (reader.ReadFromFile(helpd,"CHECK_POLES_THRESHOLD")) {
    m_checkpolesthreshold = helpd;
    msg_Tracking()<<"Set infrared poles check threshold to "<<m_checkpolesthreshold<<" . "<<std::endl;
  }
  m_checkfinite = false;
  if (reader.ReadFromFile(helpi,"CHECK_FINITE")) {
    m_checkfinite = helpi;
    msg_Tracking()<<"Set infrared finite check mode "<<m_checkfinite<<" . "<<std::endl;
  }
  if (reader.ReadFromFile(helpd,"DIPOLE_ALPHA")) {
    m_dalpha = helpd;
    msg_Tracking()<<"Set dipole cut alpha="<<m_dalpha<<" . "<<std::endl;
  }
  if (reader.ReadFromFile(helpd,"DIPOLE_ALPHA_FF")) {
    m_dalpha_ff = helpd;
    msg_Tracking()<<"Set ff dipole cut alpha="<<m_dalpha_ff<<" . "<<std::endl;
  }
  else m_dalpha_ff=m_dalpha;
  if (reader.ReadFromFile(helpd,"DIPOLE_ALPHA_FI")) {
    m_dalpha_fi = helpd;
    msg_Tracking()<<"Set fi dipole cut alpha="<<m_dalpha_fi<<" . "<<std::endl;
  }
  else m_dalpha_fi=m_dalpha;
  if (reader.ReadFromFile(helpd,"DIPOLE_ALPHA_IF")) {
    m_dalpha_if = helpd;
    msg_Tracking()<<"Set if dipole cut alpha="<<m_dalpha_if<<" . "<<std::endl;
  }
  else m_dalpha_if=m_dalpha;
  if (reader.ReadFromFile(helpd,"DIPOLE_ALPHA_II")) {
    m_dalpha_ii = helpd;
    msg_Tracking()<<"Set ii dipole cut alpha="<<m_dalpha_ii<<" . "<<std::endl;
  }
  else m_dalpha_ii=m_dalpha;
  m_force_init=reader.GetValue("LOOP_ME_INIT",0);
  m_checkloopmap=reader.GetValue("CHECK_LOOP_MAP",0);
  m_sccmur=reader.GetValue("USR_WGT_MODE",1);
  m_user_bvimode=reader.GetValue("NLO_BVI_MODE",0);
  m_imode=reader.GetValue<int>("NLO_IMODE",0);
  m_cmur[0]=0.;
  m_cmur[1]=0.;

  static bool addcite(false);
  if (!addcite) {
    addcite=true;
  rpa->gen.AddCitation(1,"The automated generation of Catani-Seymour dipole\
 terms in Amegic is published under \\cite{Gleisberg:2007md}.");
  }
}



Single_Virtual_Correction::~Single_Virtual_Correction()
{
  m_cpls.clear();
  p_selector=NULL;
  p_scale=NULL;
  if (p_psgen)      {delete p_psgen; p_psgen=0;}
  if (p_LO_process) {delete p_LO_process; p_LO_process=0;}
  if (p_dipole)     {delete p_dipole; p_dipole=0;}
  if (p_kpterms)    {delete p_kpterms; p_kpterms=0;}
  if (p_loopme)     {delete p_loopme; p_loopme=0;}
}

/*------------------------------------------------------------------------------
  
  Generic stuff for initialization of Single_Virtual_Correctiones
      
  ------------------------------------------------------------------------------*/

void Single_Virtual_Correction::PolarizationNorm() {
  m_Norm = SymmetryFactors() * p_LO_process->GetPolarisation()->Spin_Average(m_nin,&m_flavs.front());
}

double Single_Virtual_Correction::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom) {
  if (p_loopme) {
    return p_loopme->Eps_Scheme_Factor(mom);
  }
  else {
    // MSbar
    return 4.*M_PI;
    // return 2.*M_PI*p_scale->Scale(stp::ren,1)/(mom[0]*mom[1]);
  }
}

bool AMEGIC::Single_Virtual_Correction::Combinable
(const size_t &idi,const size_t &idj)
{
  return p_LO_process->Combinable(idi, idj);
}

  
const Flavour_Vector &AMEGIC::Single_Virtual_Correction::CombinedFlavour
(const size_t &idij)
{
  return p_LO_process->CombinedFlavour(idij);
}


/*------------------------------------------------------------------------------

  Initializing libraries, amplitudes, etc.

  ------------------------------------------------------------------------------*/
void Single_Virtual_Correction::SelectLoopProcess()
{
  p_loopme=NULL;
  if (m_pinfo.m_fi.m_nloqcdtype&nlo_type::loop || m_force_init) {
    Process_Info loop_pi(m_pinfo);
    loop_pi.m_fi.m_nloqcdtype=nlo_type::loop;
    p_loopme=PHASIC::Virtual_ME2_Base::GetME2(loop_pi);  
    if (!p_loopme) {
      PRINT_VAR(loop_pi);
      THROW(not_implemented, "Couldn't find virtual ME for this process.");
    }
    p_loopme->SetCouplings(*p_LO_process->CouplingMap());
    p_loopme->SetNorm(m_Norm);
  }
}



int Single_Virtual_Correction::InitAmplitude(Model_Base * model,Topology* top,
				      vector<Process_Base *> & links,
				      vector<Process_Base *> & errs)
{
  Init();
  if (m_user_bvimode!=0) m_bvimode=m_user_bvimode;
  else m_bvimode=7;
  m_eoreset = (m_bvimode!=7);
  if (!model->CheckFlavours(m_nin,m_nout,&m_flavs.front())) return 0;

  m_pslibname = ToString(m_nin)+"_"+ToString(m_nout);
  m_ptypename = "P"+m_pslibname;
//   m_name+= "_VIRT";

  if (m_pinfo.m_amegicmhv>0) {
    if (m_pinfo.m_amegicmhv==12) {
      p_LO_process = new Single_LOProcess_External(m_pinfo, p_int->Beam(), p_int->ISR());
    }
    else if (CF.MHVCalculable(m_pinfo))
      p_LO_process = new Single_LOProcess_MHV(m_pinfo, p_int->Beam(), p_int->ISR());
    if (m_pinfo.m_amegicmhv==2) return 0;
  }
  if (!p_LO_process)
    p_LO_process = new Single_LOProcess(m_pinfo, p_int->Beam(), p_int->ISR());
  p_LO_process->SetTestMoms(p_testmoms);
  p_proc=p_LO_process;

  int massive(0);
  if (m_dalpha_ii!=m_dalpha ||
      m_dalpha_if!=m_dalpha ||
      m_dalpha_fi!=m_dalpha ||
      m_dalpha_ff!=m_dalpha) massive=1;
  p_kpterms = new KP_Terms(this,massive|(m_kpcemode?2:0));
  p_kpterms->SetIMode(m_imode);
  p_flkern=p_kpterms->FlavKern();
  p_masskern=p_kpterms->MassKern();
  if (!p_masskern) p_kpterms->SetAlpha(m_dalpha);
  else p_kpterms->SetAlpha(m_dalpha_ff, m_dalpha_fi,m_dalpha_if,m_dalpha_ii);
  p_dipole = new DipoleSplitting_Base();
  if (!p_masskern) p_dipole->SetAlpha(m_dalpha);

  PolarizationNorm();

  string sfname = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename+"/"+Name()+".sym";
  if (FileExists(sfname)) {
    My_In_File in(sfname);
    if (!in.Open()) THROW(fatal_error,"Cannot open "+sfname);
    double symf;
    *in>>symf;
    m_Norm*=symf;
    in.Close();
  }

  if (!p_LO_process->InitAmplitude(model,top,links,errs,m_checkloopmap)) return 0;
  m_iresult = p_LO_process->Result();
  nlo_type::code nlot(nlo_type::loop|nlo_type::vsub);
  m_oqcd = p_LO_process->OrderQCD()+
    ((m_pinfo.m_fi.m_nloqcdtype&nlot)?1:0);
  m_oew = p_LO_process->OrderEW()+
    ((m_pinfo.m_fi.m_nloewtype&nlot)?1:0);
  m_pinfo.m_oqcd=m_oqcd;
  m_pinfo.m_oew=m_oew;

  p_dipole->SetCoupling(p_LO_process->CouplingMap());
  p_kpterms->SetCoupling(p_LO_process->CouplingMap());
  m_cpls=*p_LO_process->CouplingMap();

  if (p_LO_process!=p_LO_process->Partner()) {
    string partnerID=p_LO_process->Partner()->Name();
    for (size_t j=0;j<links.size();j++) if (Type()==links[j]->Type()) {
      if (partnerID==links[j]->Name()) {
	msg_Tracking()<<"Can map full virtual process: "<<Name()<<" -> "<<partnerID<<" Factor: "<<p_LO_process->GetSFactor()<<endl;
 	p_mapproc = p_partner = (Single_Virtual_Correction*)links[j];
	InitFlavmap(p_partner);
	m_sfactor = p_LO_process->GetSFactor();
	m_cpls=p_partner->m_cpls;
 	break;
      }
    }
  }

  if (p_partner==this) {
    SelectLoopProcess();

    links.push_back(this);
    p_dsij.resize(p_LO_process->PartonList().size());
    for (size_t i=0;i<p_LO_process->PartonList().size();i++) {
      p_dsij[i].resize(p_LO_process->PartonList().size());
      for (size_t j=0;j<p_LO_process->PartonList().size();j++) p_dsij[i][j]=0.;
    }
  }
  m_mewgtinfo.m_type=mewgttype::muR;
  if (m_pinfo.m_fi.m_nloqcdtype&nlo_type::vsub)
    m_mewgtinfo.m_type|=mewgttype::muF;
  Minimize();
  if (p_partner==this && Result()>0.) SetUpIntegrator();
  return 1;
}

void AMEGIC::Single_Virtual_Correction::AddPoint(const double &value)
{
  if (m_bvimode!=7) return;
  double last(m_lastb+m_lastv+m_lasti+m_lastkp);
  if (value!=0.0 && last==0.0) {
    msg_Error()<<METHOD<<"(): Zero result in '"<<m_name<<"'."<<std::endl;
    return;
  }
#ifdef USING__MPI
  ++m_mn;
  if (value==0.0) return;
  m_mbsum+=sqr(value*m_lastb/last);
  m_mvsum+=sqr(value*m_lastv/last);
  m_misum+=sqr(value*(m_lasti+m_lastkp)/last);
#else
  ++m_n;
  if (value==0.0) return;
  m_bsum+=sqr(value*m_lastb/last);
  m_vsum+=sqr(value*m_lastv/last);
  m_isum+=sqr(value*(m_lasti+m_lastkp)/last);
#endif
}

bool AMEGIC::Single_Virtual_Correction::ReadIn(const std::string &pid)
{
  std::string name;
  My_In_File from(pid+"/"+m_name+".bvi");
  if (!from.Open()) return false;
  from->precision(16);
  *from>>name>>m_n>>m_bsum>>m_vsum>>m_isum;
  if (name!=m_name) THROW(fatal_error,"Corrupted results file");
  return true;
}

void AMEGIC::Single_Virtual_Correction::WriteOut(const std::string &pid)
{
  My_Out_File outfile(pid+"/"+m_name+".bvi");
  outfile.Open();
  outfile->precision(16);
  *outfile<<m_name<<"  "<<m_n<<" "<<m_bsum<<" "<<m_vsum<<" "<<m_isum<<"\n";
}

void AMEGIC::Single_Virtual_Correction::EndOptimize()
{
  m_bvimode=7;
  if (m_eoreset) p_int->Reset();
}

bool AMEGIC::Single_Virtual_Correction::NewLibs() 
{
  return p_partner->GetLOProcess()->NewLibs();
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
	 (m_flavs[1].Mass() != p_int->ISR()->Flav(1).Mass()) ) p_int->ISR()->SetPartonMasses(m_flavs);
    if (CreateChannelLibrary()) return 1;
  }
  if (m_nin==1) if (CreateChannelLibrary()) return 1;
  return 0;
}

bool Single_Virtual_Correction::CreateChannelLibrary()
{
  if (!p_LO_process || p_LO_process->NumberOfDiagrams()==0) return 1;
  if (p_LO_process->Partner()!=p_LO_process || p_psgen) return true;
  p_psgen     = new Phase_Space_Generator(m_nin,m_nout);
  bool newch  = 0;
  if (m_nin>=1)  newch = p_psgen->Construct(p_channellibnames,m_ptypename,p_LO_process->PSLibName(),&m_flavs.front(),p_LO_process); 
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
  if (p_kpterms)    { delete p_kpterms; p_kpterms=0;}
  if (p_loopme)     { delete p_loopme; p_loopme=0;}

  m_oqcd      = p_partner->OrderQCD();
  m_oew       = p_partner->OrderEW();
}

/*------------------------------------------------------------------------------

  Calculating total cross sections
  
  ------------------------------------------------------------------------------*/


double Single_Virtual_Correction::Partonic(const ATOOLS::Vec4D_Vector &moms,const int mode)
{
  if (mode==1) THROW(fatal_error,"Invalid call");
  if (!Selector()->Result()) return m_lastxs = m_lastdxs = m_lastbxs = 0.0;
  return DSigma(moms,m_lookup,mode);
}

double Single_Virtual_Correction::DSigma(const ATOOLS::Vec4D_Vector &_moms,bool lookup,const int mode)
{
  m_lastxs = m_lastdxs = m_lastbxs = 0.;
  double wgt(1.0);
  int bvimode(p_partner->m_bvimode);
  if (!lookup && m_user_bvimode!=0) {
    double sum(((m_user_bvimode&1)?dabs(m_bsum):0.0)+
	       ((m_user_bvimode&2)?dabs(m_isum):0.0)+
	       ((m_user_bvimode&4)?dabs(m_vsum):0.0)), disc(sum*ran->Get());
    if (disc>dabs(m_bsum)+dabs(m_isum)) {
      p_partner->m_bvimode=4;
      wgt=sum/dabs(m_vsum);
    }
    else if (disc>dabs(m_bsum)) {
      p_partner->m_bvimode=2;
      wgt=sum/dabs(m_isum);
    }
    else {
      p_partner->m_bvimode=1;
      wgt=sum/dabs(m_bsum);
    }
  }
  if (p_partner == this) {
    m_lastdxs = operator()(_moms,mode);
  }
  else {
    if (lookup) {
      m_lastdxs = p_partner->LastDXS()*m_sfactor;
    }
    else {
      p_LO_process->Integrator()->SetMomenta(p_int->Momenta());
      m_lastdxs = p_partner->operator()(_moms,mode)*m_sfactor;
    }
    m_lastbxs = p_partner->m_lastbxs*m_sfactor;
    m_lastb=p_partner->m_lastb*m_sfactor;
    m_lastv=p_partner->m_lastv*m_sfactor;
    m_lasti=p_partner->m_lasti*m_sfactor;
  }
  double eta0(0.0), eta1(0.0);
  if (mode==0) {
    eta0=p_int->Momenta()[0].PPlus()/rpa->gen.PBeam(0).PPlus();
    eta1=p_int->Momenta()[1].PMinus()/rpa->gen.PBeam(1).PMinus();
  }
  else {
    eta0=p_int->Momenta()[0].PPlus()/rpa->gen.PBeam(1).PMinus();
    eta1=p_int->Momenta()[1].PMinus()/rpa->gen.PBeam(0).PPlus();
  }
  double kpterm = (p_partner->m_bvimode&2)?p_partner->Get_KPterms
    (p_int->ISR()->PDF(mode),p_int->ISR()->PDF(1-mode),eta0,eta1,m_flavs):0.0;
  if (p_partner != this) kpterm*=m_sfactor;
  m_lastkp=kpterm;

  m_mewgtinfo.m_w0 = m_lastdxs/m_sfactor;
  p_partner->FillMEwgts(m_mewgtinfo);
  m_mewgtinfo*=m_Norm*m_sfactor;

  p_partner->m_bvimode=bvimode;

  m_lastbxs*=m_Norm;
  return m_lastxs = wgt * m_Norm * (m_lastdxs+kpterm);
}

double Single_Virtual_Correction::Calc_Imassive(const ATOOLS::Vec4D *mom) 
{
  double res=0.;
  double mur = p_scale->Scale(stp::ren,1);
  int lm(p_loopme?p_loopme->DRMode():0);
  for (size_t i=0;i<p_LO_process->PartonList().size();i++) {
    for (size_t k=i+1;k<p_LO_process->PartonList().size();k++) {
      int typei = m_flavs[p_LO_process->PartonList()[i]].IntSpin();
      int typek = m_flavs[p_LO_process->PartonList()[k]].IntSpin();
      double sik=2.*mom[p_LO_process->PartonList()[i]]
                   *mom[p_LO_process->PartonList()[k]];
      double mi=m_flavs[p_LO_process->PartonList()[i]].Mass();
      double mk=m_flavs[p_LO_process->PartonList()[k]].Mass();
      bool susyi = m_flavs[p_LO_process->PartonList()[i]].IsSusy();
      bool susyk = m_flavs[p_LO_process->PartonList()[k]].IsSusy();

      p_masskern->Calculate(typei,mur,sik,mi,mk,
                            p_LO_process->PartonList()[i]<m_nin,
                            p_LO_process->PartonList()[k]<m_nin,susyi,lm);
      double splf  = p_masskern->I_Fin();
      double splf1 = p_masskern->I_E1();
      double splf2 = p_masskern->I_E2();
      p_masskern->Calculate(typek,mur,sik,mk,mi,
                            p_LO_process->PartonList()[k]<m_nin,
                            p_LO_process->PartonList()[i]<m_nin,susyk,lm);
      splf  += p_masskern->I_Fin();
      splf1 += p_masskern->I_E1();
      splf2 += p_masskern->I_E2();

      Vec4D_Vector momv(mom, &mom[m_nin+m_nout]);
      double lsc = log(4.*M_PI*mur/dabs(sik)/Eps_Scheme_Factor(momv));

      splf+=splf1*lsc+splf2*0.5*sqr(lsc);
      res+=p_dsij[i][k]*splf;
      m_cmur[0]+=p_dsij[i][k]*(splf1+splf2*lsc);
      m_cmur[1]+=p_dsij[i][k]*splf2;
    }
  }
  m_cmur[0]*=-p_kpterms->Coupling();
  m_cmur[1]*=-p_kpterms->Coupling();

  return -res*p_kpterms->Coupling();
}

double Single_Virtual_Correction::Calc_I(const ATOOLS::Vec4D *mom) 
{
  if (m_imode && !(m_imode&1)) return 0.;
  if (p_masskern) return Calc_Imassive(mom);

  double res=0.;
  double mur = p_scale->Scale(stp::ren,1);
  int lm(p_loopme?p_loopme->DRMode():0);
  for (size_t i=0;i<p_LO_process->PartonList().size();i++) {
    for (size_t k=i+1;k<p_LO_process->PartonList().size();k++) {
      int typei = 2*m_flavs[p_LO_process->PartonList()[i]].IntSpin();
      int typek = 2*m_flavs[p_LO_process->PartonList()[k]].IntSpin();
      Vec4D_Vector momv(mom, &mom[m_nin+m_nout]);
      double lsc = log(4.*M_PI*mur/dabs(2.*mom[p_LO_process->PartonList()[i]]*mom[p_LO_process->PartonList()[k]])/Eps_Scheme_Factor(momv));
      double splf = p_dipole->Vif(typei,lm)+p_dipole->Vif(typek,lm);
      double splf1 = p_dipole->Vie1(typei)+p_dipole->Vie1(typek);
      double splf2 = p_dipole->Vie2(typei)+p_dipole->Vie2(typek);
      splf+=splf1*lsc+splf2*0.5*sqr(lsc);
      res+=p_dsij[i][k]*splf;
      m_cmur[0]+=p_dsij[i][k]*(splf1+splf2*lsc);
      m_cmur[1]+=p_dsij[i][k]*splf2;
    }
  }
  m_cmur[0]*=-p_kpterms->Coupling();
  m_cmur[1]*=-p_kpterms->Coupling();

  return -res*p_kpterms->Coupling();
}

void Single_Virtual_Correction::Calc_KP(const ATOOLS::Vec4D *mom, double x0, double x1, double eta0, double eta1, double weight) 
{
  if (m_imode && !(m_imode&2 || m_imode&4)) return;
  p_kpterms->SetDSij(p_dsij);
  p_kpterms->Calculate(p_int->Momenta(),x0,x1,eta0,eta1,-weight*p_dipole->SPFac()/(16.0*sqr(M_PI)));
}

double Single_Virtual_Correction::Get_KPterms(PDF_Base *pdfa, PDF_Base *pdfb,
					      const double &eta0,const double &eta1,
					      ATOOLS::Flavour_Vector& flav) 
{
  if (m_imode && !(m_imode&2 || m_imode&4)) return 0.;
  if ((m_pinfo.m_fi.m_nloqcdtype&nlo_type::vsub)==0) return 0.;
  int mode(pdfa==p_int->ISR()->PDF(0)?0:1);
  return p_kpterms->Get(m_x0,m_x1,eta0,eta1,flav,mode) * KFactor();
}

void Single_Virtual_Correction::CheckPoleCancelation(const ATOOLS::Vec4D *mom)
{
  if (!p_loopme) {
    cout<<"Didn't initialise virtual ME. Ignoring pole check."<<endl;
    return;
  }
  double doublepole=0.;
  double singlepole=0.;
  double mur = p_scale->Scale(stp::ren,1);
  if (!p_masskern) {
    for (size_t i=0;i<p_LO_process->PartonList().size();i++) {
      for (size_t k=i+1;k<p_LO_process->PartonList().size();k++) {
        int typei = 2*m_flavs[p_LO_process->PartonList()[i]].IntSpin();
        int typek = 2*m_flavs[p_LO_process->PartonList()[k]].IntSpin();
        Vec4D_Vector momv(mom, &mom[m_nin+m_nout]);
        double lsc = log(4.*M_PI*mur/dabs(2.*mom[p_LO_process->PartonList()[i]]
                                            *mom[p_LO_process->PartonList()[k]])
                         /Eps_Scheme_Factor(momv));

        doublepole+=p_dsij[i][k]*(p_dipole->Vie2(typei)+p_dipole->Vie2(typek));
        singlepole+=p_dsij[i][k]*(p_dipole->Vie1(typei)+p_dipole->Vie1(typek)
                                  +(p_dipole->Vie2(typei)+p_dipole->Vie2(typek))
                                   *lsc);
      }
    }
  }
  else {
    int lm(p_loopme?p_loopme->DRMode():0);
    for (size_t i=0;i<p_LO_process->PartonList().size();i++) {
      for (size_t k=i+1;k<p_LO_process->PartonList().size();k++) {
        int typei = m_flavs[p_LO_process->PartonList()[i]].IntSpin();
        int typek = m_flavs[p_LO_process->PartonList()[k]].IntSpin();
        double sik=2.*mom[p_LO_process->PartonList()[i]]
                     *mom[p_LO_process->PartonList()[k]];
        double mi=m_flavs[p_LO_process->PartonList()[i]].Mass();
        double mk=m_flavs[p_LO_process->PartonList()[k]].Mass();
        bool susyi = m_flavs[p_LO_process->PartonList()[i]].IsSusy();
        bool susyk = m_flavs[p_LO_process->PartonList()[k]].IsSusy();

        p_masskern->Calculate(typei,mur,sik,mi,mk,
                              p_LO_process->PartonList()[i]<m_nin,
                              p_LO_process->PartonList()[k]<m_nin,
                              susyi,lm);
        double splf1 = p_masskern->I_E1();
        double splf2 = p_masskern->I_E2();
        p_masskern->Calculate(typek,mur,sik,mk,mi,
                              p_LO_process->PartonList()[k]<m_nin,
                              p_LO_process->PartonList()[i]<m_nin,
                              susyk,lm);
        splf1 += p_masskern->I_E1();
        splf2 += p_masskern->I_E2();

        Vec4D_Vector momv(mom, &mom[m_nin+m_nout]);
        double lsc = log(4.*M_PI*mur/dabs(sik)/Eps_Scheme_Factor(momv));

        splf1+=splf2*lsc;
        doublepole+=p_dsij[i][k]*splf2;
        singlepole+=p_dsij[i][k]*splf1;
      }
    }
  }
  doublepole*=m_Norm*p_kpterms->Coupling();
  singlepole*=m_Norm*p_kpterms->Coupling();
  double p1(m_Norm*p_loopme->ME_E1()*p_kpterms->Coupling()),
         p2(m_Norm*p_loopme->ME_E2()*p_kpterms->Coupling());
  if (p_loopme->Mode()==0) {
    p1*=p_dsij[0][0];
    p2*=p_dsij[0][0];
  }
  size_t precision(msg->Out().precision());
  msg->SetPrecision(16);
  if (!m_checkpolesthreshold ||
      !ATOOLS::IsEqual(doublepole,p2,m_checkpolesthreshold) ||
      !ATOOLS::IsEqual(singlepole,p1,m_checkpolesthreshold)) {
    msg_Out()<<"------------------------------------------------------------\n";
    msg_Out()<<"Process: "<<Name()<<std::endl;
    for (size_t i=0;i<m_nin+m_nout;i++) msg_Out()<<i<<": "<<mom[i]<<std::endl;
  }
  if (!m_checkpolesthreshold ||
      !ATOOLS::IsEqual(doublepole,p2,m_checkpolesthreshold)) {
    msg_Out()<<"Double poles: Sherpa = "<<doublepole<<" vs. OLP = "<<p2
             <<", rel. diff.: "<<(doublepole-p2)/(doublepole+p2)
             <<", ratio: "<<doublepole/p2<<std::endl;
  }
  if (!m_checkpolesthreshold ||
      !ATOOLS::IsEqual(singlepole,p1,m_checkpolesthreshold)) {
    msg_Out()<<"Single poles: Sherpa = "<<singlepole<<" vs. OLP = "<<p1
             <<", rel. diff.: "<<(singlepole-p1)/(singlepole+p1)
             <<", ratio: "<<singlepole/p1<<std::endl;
  }
  msg->SetPrecision(precision);
}

double Single_Virtual_Correction::operator()(const ATOOLS::Vec4D_Vector &mom,const int mode)
{
  if (p_partner!=this) {
    p_partner->Integrator()->SetMomenta(p_int->Momenta());
    return p_partner->operator()(mom,mode)*m_sfactor;
  }

  double M2(0.);
  m_cmur[0]=0.;
  m_cmur[1]=0.;

  Vec4D_Vector _mom(mom);
  Poincare cms;
  if (m_nin==2 && p_int->ISR() && p_int->ISR()->On()) {
    cms=Poincare(_mom[0]+_mom[1]);
    for (size_t i(0);i<_mom.size();++i) cms.Boost(_mom[i]);
  }
  p_LO_process->Calc_AllXS(p_int->Momenta(),&_mom.front(),p_dsij,mode);
  if (p_loopme && (m_bvimode&4)) {
    p_loopme->SetRenScale(p_scale->Scale(stp::ren,1));
    p_loopme->Calc(mom);
  }
  double I=0.;
  if ((m_pinfo.m_fi.m_nloqcdtype&nlo_type::vsub) &&
      (m_bvimode&2)) I=Calc_I(&mom.front());
  m_x0=1.,m_x1=1.;
  double eta0=1.,eta1=1.;
  double w=1.;
  if (m_flavs[0].Strong()) {
    eta0 = mom[0].PPlus()/p_int->Beam()->GetBeam(0)->OutMomentum().PPlus();
    m_x0 = eta0+ran->Get()*(1.-eta0);
    w *= (1.-eta0);
//       m_x0 = eta0*std::exp(-ran->Get()*log(eta0));
//       w *= -m_x0*log(eta0);
  }
  if (m_flavs[1].Strong()) {
    eta1 = mom[1].PMinus()/p_int->Beam()->GetBeam(1)->OutMomentum().PMinus();
    m_x1 = eta1+ran->Get()*(1.-eta1);
    w *= (1.-eta1);
//        m_x1 = eta1*std::exp(-ran->Get()*log(eta1));
//        w *= -m_x1*log(eta1);
  }

  if ((m_pinfo.m_fi.m_nloqcdtype&nlo_type::vsub) &&
      (m_bvimode&2)) Calc_KP(&mom.front(),m_x0,m_x1,eta0,eta1,w);

  double lme = 0.;
  // virtual me2 is supposed to return local nlo kfactor to born
  if ((m_pinfo.m_fi.m_nloqcdtype&nlo_type::loop) &&
      (m_bvimode&4)) {
    if (p_loopme->Mode()==0) {
      lme = p_dsij[0][0]*p_kpterms->Coupling()*p_loopme->ME_Finite();
      if (m_sccmur) {
      m_cmur[0]+=(p_loopme->ME_E1()+(OrderQCD()-1)*p_dipole->G2())*
	p_dsij[0][0]*p_kpterms->Coupling();
      m_cmur[1]+=p_loopme->ME_E2()*p_dsij[0][0]*p_kpterms->Coupling();
      }
      else {
      m_cmur[0]+=p_dsij[0][0]*p_kpterms->Coupling()*
	p_loopme->ScaleDependenceCoefficient(1);
      m_cmur[1]+=p_dsij[0][0]*p_kpterms->Coupling()*
	p_loopme->ScaleDependenceCoefficient(2);
      }
    }
    else if (p_loopme->Mode()==1) {
      lme = p_kpterms->Coupling()*p_loopme->ME_Finite();
      if (m_sccmur) {
      m_cmur[0]+=p_kpterms->Coupling()*
	(p_loopme->ME_E1()+(OrderQCD()-1)*p_dipole->G2());
      m_cmur[1]+=p_kpterms->Coupling()*p_loopme->ME_E2();
      }
      else {
      m_cmur[0]+=p_kpterms->Coupling()*
	p_loopme->ScaleDependenceCoefficient(1);
      m_cmur[1]+=p_kpterms->Coupling()*
	p_loopme->ScaleDependenceCoefficient(2);
      }
    }
    else if (p_loopme->Mode()==2) {
      // loop ME already contains I
      I=0;
      lme = p_dsij[0][0]*p_kpterms->Coupling()*p_loopme->ME_Finite();
    }
    else THROW(not_implemented,"Unknown mode");
  }
  if (m_checkpoles)
    CheckPoleCancelation(&mom.front());
  if (m_checkfinite) {
    msg->SetPrecision(16);
    msg_Out()<<"Finite_I = "<<m_Norm*I<<" vs. Finite_OLE = "<<-m_Norm*lme
             <<", rel. diff. "<<(I+lme)/(I-lme)<<std::endl;
    msg->SetPrecision(6);
  }
  if (m_checkborn &&
      (!m_checkpolesthreshold ||
       !ATOOLS::IsEqual(m_Norm*p_dsij[0][0],p_loopme->ME_Born(),m_checkpolesthreshold))) {
    msg->SetPrecision(16);
    msg_Out()<<"Born check: "
	     <<"Sherpa = "<<m_Norm*p_dsij[0][0]<<" vs. OLE = "<<p_loopme->ME_Born()
	     <<", rel. diff.: "<<(m_Norm*p_dsij[0][0]-p_loopme->ME_Born())/(m_Norm*p_dsij[0][0]+p_loopme->ME_Born())
	     <<", ratio: "<<m_Norm*p_dsij[0][0]/p_loopme->ME_Born()<<std::endl;
    msg_Out()<<"Virtual: OLE = "<<p_loopme->ME_Finite()<<std::endl;
    msg->SetPrecision(6);
  }
  double kfactor(KFactor());
  m_lastb=p_dsij[0][0] * kfactor;
  m_lasti=I            * kfactor;
  m_lastv=lme          * kfactor;
  M2+=I+lme;
  if ((m_pinfo.m_fi.m_nloqcdtype&nlo_type::born) &&
      (m_bvimode&1)) {
    M2+=p_dsij[0][0];
    m_lastbxs=p_dsij[0][0] * kfactor;
  }
  if (!(M2>0.) && !(M2<0.) && !(M2==0.)) {
    msg->SetPrecision(16);
    msg_Error()<<METHOD<<"("<<Name()<<"){\n  M2 = nan.\n  eta0 = "<<eta0
        <<" ,  m_x0 = "<<m_x0<<" ,  eta1 = "<<eta1<<" ,  m_x1 = "<<m_x1
        <<" ,  p_dsij[0][0] = "<<p_dsij[0][0]
        <<"\n  L = "<<lme<<" , I = "<<I<<"\n";
    for (size_t i=0;i<m_nin+m_nout;i++)
      msg_Error()<<"  "<<i<<": "<<mom[i]<<std::endl;
    msg_Error()<<"}\n";
    msg->SetPrecision(6);
  }
  return M2 * kfactor;
}

void Single_Virtual_Correction::FillAmplitudes(vector<METOOLS::Spin_Amplitudes>& amps,
                                               vector<vector<Complex> >& cols)
{
  p_LO_process->FillAmplitudes(amps, cols);
}



int Single_Virtual_Correction::NumberOfDiagrams() { 
  if (p_partner==this) return p_LO_process->NumberOfDiagrams(); 
  return p_partner->NumberOfDiagrams();
}

Point * Single_Virtual_Correction::Diagram(int i) { 
  if (p_partner==this) return p_LO_process->Diagram(i); 
  return p_partner->Diagram(i);
} 

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

void Single_Virtual_Correction::SetSelector(const Selector_Key &key)
{
  p_LO_process->SetSelector(key);
  p_selector=p_LO_process->Selector();
}

void Single_Virtual_Correction::SetScale(const Scale_Setter_Arguments &args)
{
  if (!p_LO_process->IsMapped()) p_LO_process->SetScale(args);
  p_scale=p_LO_process->Partner()->ScaleSetter();
}

void Single_Virtual_Correction::SetGenerator(ME_Generator_Base *const gen) 
{ 
  if (p_LO_process) p_LO_process->SetGenerator(gen);
  p_gen=gen;
}

void Single_Virtual_Correction::SetShower(PDF::Shower_Base *const ps)
{
  p_LO_process->SetShower(ps);
  p_shower=ps;
}

void Single_Virtual_Correction::SetFixedScale(const std::vector<double> &s)
{
  p_LO_process->SetFixedScale(s);
}

void Single_Virtual_Correction::SetSelectorOn(const bool on)
{
  p_LO_process->SetSelectorOn(on);
}

void Single_Virtual_Correction::FillMEwgts(ATOOLS::ME_Weight_Info& wgtinfo)
{
  wgtinfo.m_y1=m_x0;
  wgtinfo.m_y2=m_x1;
  if (wgtinfo.m_wren.size()==2)
    for (size_t i=0;i<2;i++) wgtinfo.m_wren[i]=m_cmur[i];
  if (p_kpterms) p_kpterms->FillMEwgts(wgtinfo);
}

void Single_Virtual_Correction::MPISync()
{
#ifdef USING__MPI
  int size=MPI::COMM_WORLD.Get_size();
  if (size>1) {
    double val[4];
    val[0]=m_mn;
    val[1]=m_mbsum;
    val[2]=m_mvsum;
    val[3]=m_misum;
    mpi->MPIComm()->Allreduce(MPI_IN_PLACE,val,4,MPI::DOUBLE,MPI::SUM);
    m_mn=val[0];
    m_mbsum=val[1];
    m_mvsum=val[2];
    m_misum=val[3];
  }
  m_n+=m_mn;
  m_bsum+=m_mbsum;
  m_vsum+=m_mvsum;
  m_isum+=m_misum;
  m_mn=m_mbsum=m_mvsum=m_misum=0.0;
#endif
}

Flavour Single_Virtual_Correction::ReMap(const Flavour &fl,const size_t &id) const
{
  return p_LO_process->ReMap(fl,id);
}

void Single_Virtual_Correction::FillProcessMap(NLOTypeStringProcessMap_Map *apmap)
{
  Process_Base::FillProcessMap(apmap);
  p_LO_process->SetProcMap(apmap);
}
