#include "AMEGIC++/DipoleSubtraction/Single_Virtual_Correction.H"
#include "AMEGIC++/DipoleSubtraction/DipoleSplitting_Base.H"
#include "AMEGIC++/DipoleSubtraction/Single_LOProcess_MHV.H"
#include "AMEGIC++/DipoleSubtraction/Single_LOProcess_External.H"
#include "AMEGIC++/Phasespace/Phase_Space_Generator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "PDF/Main/ISR_Handler.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PHASIC++/Main/Event_Reader.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"

#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Scoped_Settings.H"

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
  m_cmur(2,0.), m_wass(4,0.),
  m_stype(sbt::none), m_user_stype(sbt::none),
  m_user_imode(cs_itype::none), m_iresult(0.0),
  p_psgen(NULL), p_partner(this), p_LO_process(NULL),
  p_kernel_qcd(NULL), p_kernel_ew(NULL),
  p_kpterms_qcd(NULL), p_kpterms_ew(NULL), p_loopme(NULL), p_reqborn(NULL),
  m_x0(1.), m_x1(1.), m_eta0(1.), m_eta1(1.), m_z0(0.), m_z1(0.),
  m_loopmapped(false),
  m_pspisrecscheme(0), m_pspfsrecscheme(0),
  m_pspissplscheme(0), m_pspfssplscheme(0),
  m_bvimode(0),
  m_bsum(0.0), m_vsum(0.0), m_isum(0.0), m_n(0.0),
  m_mbsum(0.0), m_mvsum(0.0), m_misum(0.0), m_mn(0.0),
  m_lastb(0.0), m_lastv(0.0), m_lasti(0.0), m_lastkp(0.0),
  m_finite(0.0), m_singlepole(0.0), m_doublepole(0.0),
  p_fsmc{ NULL }
{
  m_calcv=1;
  Settings& s = Settings::GetMainSettings();
  Scoped_Settings amegicsettings{ s["AMEGIC"] };
  p_fsmc=NULL;
  m_checkborn = amegicsettings["CHECK_BORN"].Get<bool>();
  m_checkpoles = amegicsettings["CHECK_POLES"].Get<bool>();
  m_checkfinite = amegicsettings["CHECK_FINITE"].Get<bool>();
  m_checkthreshold = amegicsettings["CHECK_THRESHOLD"].Get<double>();
  m_force_init = amegicsettings["LOOP_ME_INIT"].Get<size_t>();
  m_sccmur = s["USR_WGT_MODE"].Get<bool>();
  m_murcoeffvirt = s["NLO_MUR_COEFFICIENT_FROM_VIRTUAL"].Get<bool>();
  m_user_bvimode = amegicsettings["NLO_BVI_MODE"].Get<size_t>();
  m_user_imode = s["NLO_IMODE"].Get<cs_itype::type>();
  m_user_stype = s["NLO_SUBTRACTION_MODE"].Get<sbt::subtype>();
  m_epsmode = amegicsettings["NLO_EPS_MODE"].Get<size_t>();
  m_drmode = amegicsettings["NLO_DR_MODE"].Get<size_t>();
  m_checkloopmap = amegicsettings["CHECK_LOOP_MAP"].Get<size_t>();
  m_pspisrecscheme = s["DIPOLES"]["PFF_IS_RECOIL_SCHEME"].Get<size_t>();
  m_pspfsrecscheme = s["DIPOLES"]["PFF_FS_RECOIL_SCHEME"].Get<size_t>();
  m_pspissplscheme = s["DIPOLES"]["PFF_IS_SPLIT_SCHEME"].Get<size_t>();
  m_pspfssplscheme = s["DIPOLES"]["PFF_FS_SPLIT_SCHEME"].Get<size_t>();
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
  if (p_psgen)       {delete p_psgen;       p_psgen=NULL;}
  if (p_LO_process)  {delete p_LO_process;  p_LO_process=NULL;}
  if (p_kpterms_qcd) {delete p_kpterms_qcd; p_kpterms_qcd=NULL;}
  if (p_kpterms_ew)  {delete p_kpterms_ew;  p_kpterms_ew=NULL;}
  if (p_loopme)      {delete p_loopme;      p_loopme=NULL;}
}

/*----------------------------------------------------------------------------
  
  Generic stuff for initialization of Single_Virtual_Correctiones
      
  ----------------------------------------------------------------------------*/

void Single_Virtual_Correction::PolarizationNorm()
{
  m_Norm = SymmetryFactors()
           * p_LO_process->GetPolarisation()
                         ->Spin_Average(m_nin,&m_flavs.front());
}

double Single_Virtual_Correction::Eps_Scheme_Factor
(const ATOOLS::Vec4D_Vector& mom)
{
  if (p_loopme) {
    return p_loopme->Eps_Scheme_Factor(mom);
  }
  else {
    if (m_epsmode==0) {
      // MSbar
      return 4.*M_PI;
    }
    else if (m_epsmode==1) {
      // DRED (??)
      return 2.*M_PI*p_scale->Scale(stp::ren,1)/(mom[0]+mom[1]).Abs2();
    }
    else if (m_epsmode==2) {
      // DIS
      return 2.*M_PI*p_scale->Scale(stp::ren,1)/(mom[0]*mom[1]);
    }
    else {
      THROW(fatal_error, "Unknown NLO_EPS_MODE.")
    }
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
  if (m_pinfo.m_fi.m_nlotype&nlo_type::loop || m_force_init) {
    Process_Info loop_pi(m_pinfo);
    loop_pi.m_fi.m_nlotype=nlo_type::loop;
    msg_Debugging()<<"Looking for loop\n";
    p_loopme=PHASIC::Virtual_ME2_Base::GetME2(loop_pi);
    if (!p_loopme) {
      msg_Error()<<"Could not find Loop-ME from "<<loop_pi.m_loopgenerator
                 <<" for\n"<<loop_pi<<std::endl;
      THROW(not_implemented, "Couldn't find virtual ME for this process.");
    }
    p_loopme->SetCouplings(m_cpls);
    p_loopme->SetNorm(m_Norm);
    p_loopme->SetSubType(m_stype);
    p_loopme->SetPoleCheck(m_checkpoles);
    m_drmode=p_loopme->DRMode();
  }
}



int Single_Virtual_Correction::InitAmplitude(Amegic_Model * model,Topology* top,
				      vector<Process_Base *> & links,
				      vector<Process_Base *> & errs)
{
  DEBUG_FUNC("");
  Init();
  msg_Debugging()<<m_pinfo<<std::endl;
  if (m_user_bvimode!=0) m_bvimode=m_user_bvimode;
  else m_bvimode=7;
  m_eoreset = (m_bvimode!=7);
  if (!model->p_model->CheckFlavours(m_nin,m_nout,&m_flavs.front())) return 0;

  m_pslibname = ToString(m_nin)+"_"+ToString(m_nout);
  m_ptypename = "P"+m_pslibname;

  Process_Info lopi(m_pinfo);
  // need amps with both oqcd-1 and oew-1
  std::vector<double> maxcpliqcd(m_maxcpl), mincpliqcd(m_mincpl);
  std::vector<double> maxcpliew(m_maxcpl),  mincpliew(m_mincpl);
  lopi.m_mincpl[0]=Max(0.,lopi.m_mincpl[0]-1.);
  lopi.m_mincpl[1]=Max(0.,lopi.m_mincpl[1]-1.);
  // set subtraction types generally (what may be possible)
  if (lopi.m_mincpl[0]!=m_pinfo.m_mincpl[0]) m_stype|=sbt::qcd;
  if (lopi.m_mincpl[1]!=m_pinfo.m_mincpl[1]) m_stype|=sbt::qed;
  msg_Tracking()<<"Subtraction type for "<<m_name<<": "<<m_stype<<std::endl;
  msg_Debugging()<<"Order: "<<m_mincpl<<" .. "<<m_maxcpl<<std::endl;
  if (m_stype&sbt::qcd) {
    maxcpliqcd[0]-=1; mincpliqcd[0]-=1;
  }
  else {
    maxcpliqcd[0]=0; maxcpliqcd[1]=0; mincpliqcd[0]=0; mincpliqcd[1]=0;
  }
  if (m_stype&sbt::qed) {
    maxcpliew[1]-=1; mincpliew[1]-=1;
  }
  else {
    maxcpliew[0]=0; maxcpliew[1]=0; mincpliew[0]=0; mincpliew[1]=0;
  }
  if (m_stype&sbt::qcd)
    msg_Debugging()<<"I_QCD couplings: "
                   <<mincpliqcd<<" .. "<<maxcpliqcd<<std::endl;
  if (m_stype&sbt::qed)
    msg_Debugging()<<"I_EW couplings: "
                   <<mincpliew<<" .. "<<maxcpliew<<std::endl;

  if (m_pinfo.m_amegicmhv>0) {
    if (m_pinfo.m_amegicmhv==12) {
      p_LO_process = new Single_LOProcess_External(lopi,
                                                   p_int->Beam(),
                                                   p_int->ISR(),
                                                   m_stype);
    }
    else if (CF.MHVCalculable(lopi))
      p_LO_process = new Single_LOProcess_MHV(lopi,
                                              p_int->Beam(),
                                              p_int->ISR(),
                                              m_stype);
    if (lopi.m_amegicmhv==2) return 0;
  }
  if (!p_LO_process)
    p_LO_process = new Single_LOProcess(lopi,
                                        p_int->Beam(),
                                        p_int->ISR(),
                                        m_stype);
  p_LO_process->SetTestMoms(p_testmoms);
  p_LO_process->SetPrintGraphs(lopi.m_gpath);
  p_LO_process->SetPhotonSplittingModes(m_pspissplscheme,m_pspfssplscheme);
  p_LO_process->SetMaxOrdersIQCD(maxcpliqcd);
  p_LO_process->SetMinOrdersIQCD(mincpliqcd);
  p_LO_process->SetMaxOrdersIEW(maxcpliew);
  p_LO_process->SetMinOrdersIEW(mincpliew);

  PolarizationNorm();

  if (!p_LO_process->InitAmplitude(model,top,links,errs,m_checkloopmap))
    return 0;

  m_iresult = p_LO_process->Result();
  m_cpls=*p_LO_process->CouplingMap();

  // assign subtraction types according to what is actually needed
  m_stype=p_LO_process->GetSubType();
  msg_Tracking()<<"Subtraction type for "<<m_name<<": "<<m_stype<<std::endl;

  // invalid proc if no subtraction exists
  if (m_stype==sbt::none) return 0;

  // reduce subtraction type, if only IKP_QCD or IKP_QED is requested
  m_stype&=m_user_stype;
  msg_Tracking()<<"Subtraction type for "<<m_name<<": "<<m_stype<<std::endl;

  // initialise KP-Terms and get Kernels for I-Term
  if (m_stype&sbt::qcd) {
    p_kpterms_qcd = new KP_Terms(this,
                                 sbt::qcd,
                                 p_LO_process->PartonListQCD());
    p_kpterms_qcd->SetIType(m_user_imode);
    p_kpterms_qcd->SetCoupling(p_LO_process->CouplingMap());
    p_kernel_qcd=p_kpterms_qcd->Kernel();
  }
  if (m_stype&sbt::qed) {
    p_kpterms_ew = new KP_Terms(this,
                                sbt::qed,
                                p_LO_process->PartonListQED());
    p_kpterms_ew->SetIType(m_user_imode);
    p_kpterms_ew->SetCoupling(p_LO_process->CouplingMap());
    p_kernel_ew=p_kpterms_ew->Kernel();
  }

  // re-enable once we NLO-merge
//  m_maxcpl[0] = p_LO_process->MaxOrder(0)+((m_stype&sbt::qcd)?1:0);
//  m_mincpl[0] = p_LO_process->MinOrder(0)+((m_stype&sbt::qcd)?1:0);
//  m_maxcpl[1] = p_LO_process->MaxOrder(1)+((m_stype&sbt::qed)?1:0);
//  m_mincpl[1] = p_LO_process->MinOrder(1)+((m_stype&sbt::qcd)?1:0);
//  m_pinfo.m_mincpl.resize(m_mincpl.size());
//  m_pinfo.m_maxcpl.resize(m_maxcpl.size());
//  for (size_t i(0);i<m_mincpl.size();++i) m_pinfo.m_mincpl[i]=m_mincpl[i];
//  for (size_t i(0);i<m_maxcpl.size();++i) m_pinfo.m_maxcpl[i]=m_maxcpl[i];
  msg_Debugging()<<"couplings are "<<m_mincpl<<" .. "<<m_maxcpl<<std::endl;

  SelectLoopProcess();

  // if LO-proc is mapped, see whether V can be mapped as well
  if (p_LO_process!=p_LO_process->Partner()) {
    msg_Debugging()<<Name()<<": "<<this<<std::endl;
    if (p_loopme) {
      m_loopmapped=p_loopme->IsMappableTo(p_LO_process->Partner()->Info());
      msg_Tracking()<<(m_loopmapped?"Can":"Cannot")<<" map one-loop process: "
                      <<p_loopme->Name()<<std::endl;
      if (m_loopmapped) { delete p_loopme; p_loopme=NULL;}
    }
    string partnerID=p_LO_process->Partner()->Name();
    for (size_t j=0;j<links.size();j++) if (Type()==links[j]->Type()) {
      if (partnerID==links[j]->Name()) {
        p_mapproc = p_partner = (Single_Virtual_Correction*)links[j];
        InitFlavmap(p_partner);
        m_sfactor = p_LO_process->GetSFactor();
        m_cpls=p_partner->m_cpls;
        if (p_loopme && !m_loopmapped) p_loopme->SetCouplings(m_cpls);
        break;
      }
    }
    msg_Tracking()<<"Can map LO process: "<<Name()
                  <<" -> "<<partnerID
                  <<" Factor: "<<p_LO_process->GetSFactor()<<endl;
  }
  else {
    p_LO_process->SetMinOrders(m_mincpl);
    p_LO_process->SetMaxOrders(m_maxcpl);
  }

  // if not mapped, init m_dsij
  if (p_partner==this) {
    links.push_back(this);
    size_t nijqcd(p_LO_process->PartonListQCD().size());
    m_dsijqcd.resize(nijqcd);
    for (size_t i=0;i<nijqcd;i++) {
      m_dsijqcd[i].resize(nijqcd);
      for (size_t j=0;j<nijqcd;j++) {
        m_dsijqcd[i][j]=0.;
      }
    }
    size_t nijew(p_LO_process->PartonListQED().size());
    m_dsijew.resize(nijew);
    for (size_t i=0;i<nijew;i++) {
      m_dsijew[i].resize(nijew);
      for (size_t j=0;j<nijew;j++) {
        m_dsijew[i][j]=0.;
      }
    }
    // init charge factors
    m_Q2ij.resize(nijew);
    for (size_t i=0;i<nijew;i++) {
      m_Q2ij[i].resize(nijew);
      for (size_t j=0;j<nijew;j++) {
        m_Q2ij[i][j]=0.;
      }
    }
    ComputeChargeFactors();
    // set requested Born to correct entry
    if      (m_stype==sbt::qcd)               p_reqborn=&m_dsijqcd[0][0];
    else if (m_stype==sbt::qed)               p_reqborn=&m_dsijew[0][0];
    else if (m_stype==(sbt::qcd|sbt::qed)) {
      if      (m_pinfo.m_fi.m_nlocpl[0]==1 &&
               m_pinfo.m_fi.m_nlocpl[1]==0)   p_reqborn=&m_dsijqcd[0][0];
      else if (m_pinfo.m_fi.m_nlocpl[0]==0 &&
               m_pinfo.m_fi.m_nlocpl[1]==1)   p_reqborn=&m_dsijew[0][0];
      else THROW(fatal_error,"Cannot assign correct Born for mixed correction.");
    }
    else THROW(fatal_error,"Cannot assign correct Born.");
  }
  else {
    p_reqborn=p_partner->RequestedBorn();
  }

  // init ME_Weight_Info
  // wgts needed: 2 for muR qcd, 0 for muR ew, 16 for muF qcd, 16 for muF ew
  m_mewgtinfo.m_relqcdcpl=m_pinfo.m_fi.m_nlocpl[0];
  if (m_pinfo.m_fi.m_nlotype&nlo_type::born)
    m_mewgtinfo.m_type|=mewgttype::B;
  if (m_pinfo.m_fi.m_nlotype&nlo_type::loop ||
      (m_pinfo.m_fi.m_nlotype&nlo_type::vsub && m_user_imode&cs_itype::I))
    m_mewgtinfo.m_type|=mewgttype::VI;
  if (m_pinfo.m_fi.m_nlotype&nlo_type::vsub && (m_user_imode&cs_itype::K ||
                                                m_user_imode&cs_itype::P))
    m_mewgtinfo.m_type|=mewgttype::KP;
  Minimize();
  if (p_partner==this && (Result()>0. || Result()<0.)) SetUpIntegrator();
  if (p_partner==this) msg_Info()<<"."<<std::flush;

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
  My_In_File::OpenDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/");
  if (!SetUpIntegrator()) THROW(fatal_error,"No integrator");
  My_In_File::CloseDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/");
  RequestVariables(psh);
  return Process_Base::FillIntegrator(psh);
}

void Single_Virtual_Correction::RequestVariables(Phase_Space_Handler *const psh)
{
  p_fsmc=psh->ISRIntegrator();
  if (p_fsmc==NULL) return;
  p_fsmc->AddERan("z_1");
  p_fsmc->AddERan("z_2");
} 

bool Single_Virtual_Correction::SetUpIntegrator() 
{  
  if (m_nin==2) {
    if ( (m_flavs[0].Mass() != p_int->ISR()->Flav(0).Mass()) ||
         (m_flavs[1].Mass() != p_int->ISR()->Flav(1).Mass()) )
      p_int->ISR()->SetPartonMasses(m_flavs);
    if (CreateChannelLibrary()) return 1;
  }
  if (m_nin==1) if (CreateChannelLibrary()) return 1;
  return 0;
}

bool Single_Virtual_Correction::CreateChannelLibrary()
{
  if (!p_LO_process || p_LO_process->NumberOfDiagrams()==0) return 1;
  if (p_LO_process->Partner()!=p_LO_process || p_psgen) return true;
  p_psgen     = new Phase_Space_Generator(m_nin, m_nout);
  bool newch  = 0;
  if (m_nin>=1)  newch = p_psgen->Construct(p_channellibnames,m_ptypename,
                                            p_LO_process->PSLibName(),
                                            &m_flavs.front(),p_LO_process);
  if (newch>0) return 0;
  return 1;
}

/*----------------------------------------------------------------------------
  
  Process management
  
  ----------------------------------------------------------------------------*/
void Single_Virtual_Correction::SetLookUp(const bool lookup)
{
  m_lookup=lookup; 
  if (p_LO_process) p_LO_process->SetLookUp(lookup);
  if (p_loopme && lookup==0) p_loopme->SwitchMode(lookup);
}

void Single_Virtual_Correction::Minimize()
{
  if (p_partner==this) return;
  if (p_psgen)                  { delete p_psgen;       p_psgen=NULL; }
  if (p_kpterms_qcd)            { delete p_kpterms_qcd; p_kpterms_qcd=NULL;}
  if (p_kpterms_ew)             { delete p_kpterms_ew;  p_kpterms_ew=NULL;}
  if (p_loopme && m_loopmapped) { delete p_loopme;      p_loopme=NULL;}

  m_maxcpl     = p_partner->MaxOrders();
  m_mincpl     = p_partner->MinOrders();
}

/*----------------------------------------------------------------------------

  Calculating total cross sections
  
  ----------------------------------------------------------------------------*/


double Single_Virtual_Correction::Partonic(const Vec4D_Vector &moms,
                                           Variations_Mode varmode,
                                           int mode)
{
  if (mode==1) THROW(fatal_error,"Invalid call");
  if (!Selector()->Result()) return m_lastxs = m_lastdxs = m_lastbxs = 0.0;
  return DSigma(moms,m_lookup,varmode,mode);
}

double Single_Virtual_Correction::DSigma(const Vec4D_Vector &_moms,
                                         bool lookup,
                                         Variations_Mode varmode,
                                         const int mode)
{
  DEBUG_FUNC(Name());
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
    m_lastdxs = m_Norm*operator()(_moms,varmode,mode);
  }
  else {
    if (lookup) {
      m_lastdxs = p_partner->LastDXS()*m_sfactor;
      m_lastbxs = p_partner->m_lastbxs*m_sfactor / p_partner->m_Norm;
    }
    else {
      p_LO_process->Integrator()->SetMomenta(p_int->Momenta());
      if (!m_loopmapped) p_partner->SetCalcV(0);
      m_lastdxs = m_Norm*p_partner->operator()(_moms,varmode,mode)*m_sfactor;
      p_partner->SetCalcV(1);
      m_lastbxs = p_partner->m_lastbxs*m_sfactor;
    }
    m_lastb=p_partner->m_lastb*m_sfactor;
    m_lastv=Calc_V_WhenMapped(_moms, varmode);
    m_lasti=p_partner->m_lasti*m_sfactor;
    for (size_t i(0);i<m_wass.size();++i)
      m_wass[i]=p_partner->m_wass[i]*m_sfactor;
    if (!m_loopmapped) {
      if (m_checkpoles)  {
        m_finite     = p_partner->Finite()*m_sfactor;
        m_singlepole = p_partner->SinglePole()*m_sfactor;
        m_doublepole = p_partner->DoublePole()*m_sfactor;
        CheckPoleCancelation(_moms);
      }
      // m_lasti, m_lastv already have kfactor applied
      if (m_checkfinite) CheckFinite(m_lasti,m_lastv);
      if (m_checkborn)   CheckBorn();
      // recalculate m_lastdxs
      m_lastdxs = m_lastbxs+m_lastv+m_lasti;
      m_lastdxs *= m_Norm;
    }
  }

  m_mewgtinfo.m_B = m_lastbxs/m_sfactor;
  m_mewgtinfo.m_VI = (m_lastv+m_lasti)/m_sfactor;
  for (size_t i=0;i<m_mewgtinfo.m_wass.size();++i)
    m_mewgtinfo.m_wass[i]=m_wass[i]/m_sfactor;
  p_partner->FillMEwgts(m_mewgtinfo);
  m_mewgtinfo*=m_Norm*m_sfactor;
  m_mewgtinfo.m_K = p_partner->LastK();

  const double kpterm(KPTerms(0,p_int->ISR()->PDF(0),p_int->ISR()->PDF(1)));
  m_lastkp = kpterm;
  m_mewgtinfo.m_KP = kpterm;

  p_partner->m_bvimode=bvimode;

  m_lastbxs*=m_Norm;
  size_t precision(msg->Out().precision());
  msg->SetPrecision(16);
  DEBUG_VAR(m_lastb);
  DEBUG_VAR(m_lastv);
  DEBUG_VAR(m_lasti);
  DEBUG_VAR(m_lastkp);
  msg->SetPrecision(precision);
  return m_lastxs = wgt * ( m_lastdxs + m_lastkp );
}

double Single_Virtual_Correction::Calc_B()
{
  return *p_reqborn;
}

double Single_Virtual_Correction::Calc_V(const ATOOLS::Vec4D_Vector &mom,
                                         Variations_Mode varmode)
{
  if (m_calcv==0) return 0.0;
  DEBUG_FUNC(p_loopme);
  double res(0.);
  if (!p_loopme) THROW(fatal_error,"No loop ME set.");
  p_loopme->SetRenScale(p_scale->Scale(stp::ren,1));
  p_loopme->SetPList(&p_LO_process->PartonListQCD());
  p_loopme->SetDSij(&m_dsijqcd);
  p_loopme->SetCalcAssContribs(varmode != Variations_Mode::nominal_only);
  p_loopme->Calc(mom);
  p_loopme->SetCalcAssContribs(true);
  double cplfac(1.), bornorderqcd(0), beta0qcd(0.);
  // assume alpha_qed fixed for now
  // -> otherwise add term m*beta0(QED)
  if (m_stype&sbt::qcd) {
    cplfac=p_partner->KPTermsQCD()->Coupling();
    beta0qcd=p_partner->KernelQCD()->Beta0QCD();
    bornorderqcd=MaxOrder(0)-1.;
  }
  else {
    // if no QCD subtraction exist, no QCD counterterm
    cplfac=p_partner->KPTermsEW()->Coupling();
  }
  msg_Debugging()<<"cpl="<<cplfac<<", Born O(as)="<<bornorderqcd
                 <<", \\beta_0(QCD)="<<beta0qcd<<std::endl;
  // p_loopme->Mode()&1 -> is returning full Re(M_B M_V^*)
  // virtual me2 is returning local nlo kfactor to born -> needs coupling
  if (p_loopme->Mode()==0) {
    res = m_lastb*cplfac*p_loopme->ME_Finite();
    if (m_murcoeffvirt && p_loopme->ProvidesPoles()) {
      if (m_sccmur) {
        m_cmur[0]+=(p_loopme->ME_E1()+bornorderqcd*beta0qcd)*m_lastb*cplfac;
        m_cmur[1]+=p_loopme->ME_E2()*m_lastb*cplfac;
      }
      else {
        m_cmur[0]+=m_lastb*cplfac*p_loopme->ScaleDependenceCoefficient(1);
        m_cmur[1]+=m_lastb*cplfac*p_loopme->ScaleDependenceCoefficient(2);
      }
    }
  }
  // virtual me2 is returning full Re(M_B M_V^*)
  else if (p_loopme->Mode()==1) {
    res = cplfac*p_loopme->ME_Finite();
    if (m_murcoeffvirt && p_loopme->ProvidesPoles()) {
      if (m_sccmur) {
        m_cmur[0]+=(p_loopme->ME_E1()+bornorderqcd*beta0qcd)*cplfac;
        m_cmur[1]+=p_loopme->ME_E2()*cplfac;
      }
      else {
        m_cmur[0]+=cplfac*p_loopme->ScaleDependenceCoefficient(1);
        m_cmur[1]+=cplfac*p_loopme->ScaleDependenceCoefficient(2);
      }
    }
  }
  else if (p_loopme->Mode()==2) {
    // loop ME already contains I
    res = m_lastb*cplfac*p_loopme->ME_Finite();
  }
  else if (p_loopme->Mode()==3) {
    // loop ME already contains I and is returning full Re(M_B M_V^*)+I
    res = cplfac*p_loopme->ME_Finite();
  }
  else THROW(not_implemented,"Unknown mode");
  return res;
}

double Single_Virtual_Correction::Calc_V_WhenMapped
(const Vec4D_Vector &mom, Variations_Mode varmode)
{
  if (m_loopmapped) return p_partner->m_lastv*m_sfactor;
  if ((m_stype!=sbt::none) && (m_pinfo.m_fi.m_nlotype&nlo_type::loop) &&
      (m_bvimode&4)) {
    Vec4D_Vector _mom(mom);
    Poincare cms;
    if (m_nin==2 && ((p_int->ISR() && p_int->ISR()->On()) || p_int->Beam()->On())) {
      cms=Poincare(_mom[0]+_mom[1]);
      for (size_t i(0);i<_mom.size();++i) cms.Boost(_mom[i]);
    }
    m_dsijqcd=p_partner->m_dsijqcd;
    return Calc_V(_mom, varmode);
  }
  return 0.;
}

double Single_Virtual_Correction::Calc_I(const ATOOLS::Vec4D_Vector &mom)
{
  DEBUG_FUNC("mode="<<(p_loopme?p_loopme->Mode():0));
  if (p_loopme && p_loopme->Mode()&2) return 0.;
  if (!(m_user_imode&cs_itype::I)) return 0.;
  m_finite=m_singlepole=m_doublepole=0.;
  if (m_stype&sbt::qcd) Calc_I(sbt::qcd,p_LO_process->PartonListQCD(),
                               p_kernel_qcd,p_kpterms_qcd,mom,m_dsijqcd);
  if (m_stype&sbt::qed) Calc_I(sbt::qed,p_LO_process->PartonListQED(),
                               p_kernel_ew,p_kpterms_ew,mom,m_dsijew);
  m_cmur[0]=m_singlepole;
  m_cmur[1]=m_doublepole;
  msg_Debugging()<<"I_fin = "<<m_Norm*m_finite<<std::endl;
  msg_Debugging()<<"I_e1  = "<<m_Norm*m_singlepole<<std::endl;
  msg_Debugging()<<"I_e2  = "<<m_Norm*m_doublepole<<std::endl;
  return m_finite;
}

double Single_Virtual_Correction::Calc_I(const ATOOLS::sbt::subtype st,
                                         const std::vector<size_t>& partonlist,
                                         PHASIC::Massive_Kernels * kernel,
                                         const PHASIC::KP_Terms * kpterms,
                                         const ATOOLS::Vec4D_Vector &mom,
                                         std::vector<std::vector<double> >& dsij)
{
  DEBUG_FUNC(st);
  double mur2(p_scale->Scale(stp::ren,1));
  double finite(0.), singlepole(0.), doublepole(0.);
  for (size_t i=0;i<partonlist.size();i++)
    msg_Debugging()<<m_flavs[partonlist[i]]<<": "
                   <<mom[partonlist[i]]<<std::endl;
  for (size_t i=0;i<partonlist.size();i++) {
    for (size_t k=i+1;k<partonlist.size();k++) {
      if (dsij[i][k]==0. && dsij[k][i]==0.) continue;
      ist::itype typei = AssignType(partonlist[i],st);
      ist::itype typek = AssignType(partonlist[k],st);
      double sik=2.*mom[partonlist[i]]
                   *mom[partonlist[k]];
      double mi=m_flavs[partonlist[i]].Mass();
      double mk=m_flavs[partonlist[k]].Mass();
      bool inii = partonlist[i]<m_nin;
      bool inik = partonlist[k]<m_nin;

      // I_ik
      double splf(0.),splf1(0.),splf2(0.);
      if (dsij[i][k]!=0.) {
        kernel->Calculate(typei,mur2,sik,mi,mk,inii,inik,m_drmode);
        splf  += kernel->I_Fin() * dsij[i][k];
        splf1 += kernel->I_E1() * dsij[i][k];
        splf2 += kernel->I_E2() * dsij[i][k];
        msg_Debugging()<<"I_"<<partonlist[i]<<partonlist[k]
                       <<"("<<m_flavs[partonlist[i]]<<","
                       <<m_flavs[partonlist[k]]<<")\n";
        msg_Debugging()<<"    splf_e2  = "<<kernel->I_E2()<<std::endl;
        msg_Debugging()<<"    splf_e1  = "<<kernel->I_E1()<<std::endl;
        msg_Debugging()<<"    splf_fin = "<<kernel->I_Fin()<<std::endl;
      }
      // I_ki
      if (dsij[k][i]!=0.) {
        kernel->Calculate(typek,mur2,sik,mk,mi,inik,inii,m_drmode);
        splf  += kernel->I_Fin() * dsij[k][i];
        splf1 += kernel->I_E1() * dsij[k][i];
        splf2 += kernel->I_E2() * dsij[k][i];
        msg_Debugging()<<"I_"<<partonlist[k]<<partonlist[i]
                       <<"("<<m_flavs[partonlist[k]]<<","
                       <<m_flavs[partonlist[i]]<<")\n";
        msg_Debugging()<<"    splf_e2  = "<<kernel->I_E2()<<std::endl;
        msg_Debugging()<<"    splf_e1  = "<<kernel->I_E1()<<std::endl;
        msg_Debugging()<<"    splf_fin = "<<kernel->I_Fin()<<std::endl;
      }


      double lsc(0.);
      if (!p_loopme || !(p_loopme->fixedIRscale())) 
        lsc = log(4.*M_PI*mur2/dabs(sik)/Eps_Scheme_Factor(mom));
      else{
        double irscale=p_loopme->IRscale();
        lsc = log(4.*M_PI*sqr(irscale)/dabs(sik)/Eps_Scheme_Factor(mom));
      }

      //double lsc(log(4.*M_PI*mur2/dabs(sik)/Eps_Scheme_Factor(mom)));
      msg_Debugging()<<"lsc="<<lsc<<std::endl;
      splf  += splf1*lsc + splf2*0.5*sqr(lsc);
      splf1 += splf2*lsc;

      finite     += splf;
      singlepole += splf1;
      doublepole += splf2;
    }
  }
  double cpl(kpterms->Coupling());
  finite*=-cpl;
  singlepole*=-cpl;
  doublepole*=-cpl;
  msg_Debugging()<<"I_"<<st<<"_fin = "<<m_Norm*finite<<std::endl;
  msg_Debugging()<<"I_"<<st<<"_e1  = "<<m_Norm*singlepole<<std::endl;
  msg_Debugging()<<"I_"<<st<<"_e2  = "<<m_Norm*doublepole<<std::endl;
  m_finite+=finite;
  m_singlepole+=singlepole;
  m_doublepole+=doublepole;
  return finite;
}

void Single_Virtual_Correction::Calc_KP(const ATOOLS::Vec4D_Vector &mom)
{
  if (!((m_user_imode&cs_itype::K) || (m_user_imode&cs_itype::P))) return;
  if (!p_LO_process->HasInitialStateEmitter()) return;
  DEBUG_FUNC("");
  m_x0=1.,m_x1=1.,m_eta0=1.,m_eta1=1.;
  double weight(1.);
  // dice eta0 and eta1, incorporate phase space weight in weight
  if (p_int->ISR()->PDF(0) && p_int->ISR()->PDF(0)->Contains(m_flavs[0])) {
    m_eta0=mom[0][3]>0.0?mom[0].PPlus()/rpa->gen.PBunch(0).PPlus():
      mom[0].PMinus()/rpa->gen.PBunch(1).PMinus();
    if (m_z0>0.) m_x0 = m_z0;
    else {
         m_x0 = m_eta0+p_fsmc->ERan("z_1")*(1.-m_eta0);
	 if (p_read) {
		 m_x0=m_eta0+p_read->SubEvt()->m_x1*(1.-m_eta0);
		 msg_Debugging()<<"read in x0 = "<<m_x0
			 <<" ("<<p_read->SubEvt()->m_x1<<") "<<m_eta0<<"\n";
	 }
    }
    weight*=(1.-m_eta0);
    msg_Debugging()<<"x0="<<m_x0<<std::endl;
  }
  if (p_int->ISR()->PDF(1) && p_int->ISR()->PDF(1)->Contains(m_flavs[1])) {
    m_eta1=mom[1][3]<0.0?mom[1].PMinus()/rpa->gen.PBunch(1).PMinus():
      mom[1].PPlus()/rpa->gen.PBunch(0).PPlus();
    if (m_z1>0.) m_x1 = m_z1;
    else {
        m_x1 = m_eta1+p_fsmc->ERan("z_2")*(1.-m_eta1);
	 if (p_read) {
		 m_x1=m_eta1+p_read->SubEvt()->m_x2*(1.-m_eta1);
		 msg_Debugging()<<"read in x1 = "<<m_x1
			 <<" ("<<p_read->SubEvt()->m_x2<<") "<<m_eta1<<"\n";
	 }
    }
    weight*=(1.-m_eta1);
    msg_Debugging()<<"x1="<<m_x1<<std::endl;
  }
  if (p_kpterms_qcd && p_LO_process->HasInitialStateQCDEmitter())
    p_kpterms_qcd->Calculate(p_int->Momenta(),m_dsijqcd,
                             m_x0,m_x1,m_eta0,m_eta1,weight);
  if (p_kpterms_ew  && p_LO_process->HasInitialStateQEDEmitter())
    p_kpterms_ew->Calculate(p_int->Momenta(),m_dsijew,
                            m_x0,m_x1,m_eta0,m_eta1,weight);
}

double Single_Virtual_Correction::KPTerms
(int mode, PDF::PDF_Base *pdfa, PDF::PDF_Base *pdfb, double scalefac2)
{
  // determine momentum fractions
  double eta0(0.), eta1(0.);
  if (mode == 0) {
    if (p_int->Momenta()[0][3]>0.0)
      eta0=p_int->Momenta()[0].PPlus()/rpa->gen.PBunch(0).PPlus();
    else eta0=p_int->Momenta()[0].PMinus()/rpa->gen.PBunch(1).PMinus();
    if (p_int->Momenta()[1][3]<0.0)
      eta1=p_int->Momenta()[1].PMinus()/rpa->gen.PBunch(1).PMinus();
    else eta1=p_int->Momenta()[1].PPlus()/rpa->gen.PBunch(0).PPlus();
  }
  else THROW(fatal_error,"Invalid call");
  // determine KP terms
  double kpterm(0.);
  if (p_partner->m_bvimode & 2) {
    kpterm = p_partner->Get_KPTerms(pdfa,
                                    pdfb,
                                    eta0, eta1,
                                    m_flavs[0],m_flavs[1],
                                    scalefac2);
  }
  // return normalised result
  if (p_partner != this) {
    kpterm *= m_sfactor;
  }
  return m_Norm * kpterm;
}

double Single_Virtual_Correction::Get_KPTerms(PDF_Base *pdfa, PDF_Base *pdfb,
                                              const double &eta0,
                                              const double &eta1,
                                              const ATOOLS::Flavour& fl0,
                                              const ATOOLS::Flavour& fl1,
                                              const double& sf)
{
  if (!((m_user_imode&cs_itype::K) || (m_user_imode&cs_itype::P))) return 0.;
  if (!p_LO_process->HasInitialStateEmitter()) return 0.;
  DEBUG_FUNC("");
  if (m_stype==sbt::none || !(m_pinfo.m_fi.m_nlotype&nlo_type::vsub))
    return 0.;
  double res(0.);
  double muf2(ScaleSetter()->Scale(stp::fac,1));
  if (p_kpterms_qcd && p_LO_process->HasInitialStateQCDEmitter())
    res+=p_kpterms_qcd->Get(pdfa,pdfb,m_x0,m_x1,eta0,eta1,muf2,muf2,sf,sf,fl0,fl1);
  if (p_kpterms_ew  && p_LO_process->HasInitialStateQEDEmitter())
    res+=p_kpterms_ew->Get(pdfa,pdfb,m_x0,m_x1,eta0,eta1,muf2,muf2,sf,sf,fl0,fl1);
  return res * p_partner->m_lastki;
}

void Single_Virtual_Correction::CheckBorn()
{
  if (!p_loopme) {
    msg_Info()<<"Didn't initialise virtual ME. Ignoring Born check."<<std::endl;
    return;
  }
  if (!p_reqborn) {
    msg_Info()<<"Didn't set Born ME. Ignoring Born check."<<std::endl;
    return;
  }
  double born(*p_reqborn);
  double sb(m_Norm*born), olpb(p_loopme->ME_Born());
  if (!m_checkthreshold || !ATOOLS::IsEqual(sb,olpb,m_checkthreshold)) {
    size_t precision(msg->Out().precision());
    msg->SetPrecision(16);
    msg_Out()<<"Born:         "
             <<"Sherpa = "<<sb<<" vs. OLP = "<<olpb
             <<"\n              rel. diff.: "<<(sb-olpb)/(sb+olpb)
             <<", ratio: "<<sb/olpb<<std::endl;
    msg->SetPrecision(precision);
  }
}

void Single_Virtual_Correction::CheckFinite(const double & I, const double & L)
{
  size_t precision(msg->Out().precision());
  msg->SetPrecision(16);
  // assume OLP returns I
  msg_Out()<<"Finite:       "
           <<"Sherpa = "<<m_Norm*I<<" vs. OLP = "<<-m_Norm*L
           <<", rel. diff. "<<(I-L)/(I+L)
           <<", ratio: "<<-I/L<<std::endl;
  msg->SetPrecision(precision);
}

void Single_Virtual_Correction::CheckPoleCancelation(const ATOOLS::Vec4D_Vector mom)
{
  if (m_calcv==0) return;
  DEBUG_FUNC(Name());
  if (!p_loopme) {
    msg_Info()<<"Didn't initialise virtual ME. Ignoring pole check."<<std::endl;
    return;
  }
  // multiply with norm and take out "-"sign
  double singlepole(-m_Norm*m_singlepole);
  double doublepole(-m_Norm*m_doublepole);
  double p1(m_Norm*p_loopme->ME_E1()),
         p2(m_Norm*p_loopme->ME_E2());
  double cplfac(1.);
  if (m_stype&sbt::qcd) cplfac=p_partner->KPTermsQCD()->Coupling();
  else                  cplfac=p_partner->KPTermsEW()->Coupling();
  if (p_loopme->Mode()==0) {
    p1*=cplfac*m_lastb;
    p2*=cplfac*m_lastb;
  }
  else {
    p1*=cplfac;
    p2*=cplfac;
  }
  size_t precision(msg->Out().precision());
  msg->SetPrecision(16);
  if (!m_checkthreshold || !ATOOLS::IsEqual(doublepole,p2,m_checkthreshold) ||
                           !ATOOLS::IsEqual(singlepole,p1,m_checkthreshold)) {
    msg_Out()<<"------------------------------------------------------------\n";
    msg_Out()<<"Process: "<<Name();
    if (p_partner!=this) {
      msg_Out()<<" -> "<<p_partner->Name();
      if (!m_loopmapped) msg_Out()<<" (loop not mapped)";
    }
    msg_Out()<<std::endl;
    for (size_t i=0;i<mom.size();i++) msg_Out()<<i<<": "<<mom[i]<<std::endl;
  }
  if (!m_checkthreshold || !ATOOLS::IsEqual(doublepole,p2,m_checkthreshold)) {
    msg_Out()<<"Double poles: "
             <<"Sherpa = "<<doublepole<<" vs. OLP = "<<p2
             <<"\n              rel. diff.: "<<(doublepole-p2)/(doublepole+p2)
             <<", ratio: "<<doublepole/p2<<std::endl;
  }
  if (!m_checkthreshold || !ATOOLS::IsEqual(singlepole,p1,m_checkthreshold)) {
    msg_Out()<<"Single poles: "
             <<"Sherpa = "<<singlepole<<" vs. OLP = "<<p1
	     <<"\n              rel. diff.: "<<(singlepole-p1)/(singlepole+p1)
             <<", ratio: "<<singlepole/p1<<std::endl;
  }
  msg->SetPrecision(precision);
}

double Single_Virtual_Correction::operator()(const ATOOLS::Vec4D_Vector &mom,
                                             Variations_Mode varmode,
                                             const int mode)
{
  DEBUG_FUNC("bvimode="<<m_bvimode);
  m_lastxs = m_lastdxs = m_lastbxs = 0.;
  if (p_partner!=this) {
    p_partner->Integrator()->SetMomenta(p_int->Momenta());
    return p_partner->operator()(mom,varmode,mode)*m_sfactor;
  }

  double B(0.),V(0.),I(0.);
  m_cmur[0]=m_cmur[1]=0.;

  Vec4D_Vector _mom(mom);
  Poincare cms;
  if (m_nin==2 && ((p_int->ISR() && p_int->ISR()->On()) || p_int->Beam()->On())) {
    cms=Poincare(_mom[0]+_mom[1]);
    for (size_t i(0);i<_mom.size();++i) cms.Boost(_mom[i]);
  }
  size_t precision(msg->Out().precision());
  msg->SetPrecision(16);
  msg_Debugging()<<"CMS momenta"<<std::endl;
  for (size_t i(0);i<_mom.size();++i) msg_Debugging()<<"i: "<<_mom[i]<<std::endl;
  msg->SetPrecision(precision);
  size_t num(0);
  p_LO_process->Calc_AllXS(p_int->Momenta(),&_mom.front(),
                           m_dsijqcd,m_dsijew,mode);
  AttachChargeFactors();
  PrintDSij();
  double kfactorb((m_bvimode&1)?KFactor(1|2):1.0);
  double kfactorvi(m_lastki=m_lastk=(m_bvimode&6)?KFactor((m_bvimode&1)?0:2):1.0);
  m_lastb=Calc_B()*kfactorb;

  if ((m_stype!=sbt::none) && (m_pinfo.m_fi.m_nlotype&nlo_type::born) &&
      (m_bvimode&1)) {
    m_lastk=kfactorb;
    B=m_lastb*m_lastk;
  }

  if ((m_stype!=sbt::none) && (m_pinfo.m_fi.m_nlotype&nlo_type::vsub) &&
      (m_bvimode&2)) {
    m_lastki=kfactorb;
    I=Calc_I(mom)*m_lastki;
    Calc_KP(mom);
  }

  if ((m_stype!=sbt::none) && (m_pinfo.m_fi.m_nlotype&nlo_type::loop) &&
      (m_bvimode&4)) {
    m_lastki=kfactorb;
    V=Calc_V(mom, varmode)*m_lastki;
  }

  if (p_loopme)
    for (size_t i(0);i<p_loopme->ME_AssContribs_Size();++i)
      m_wass[i]=m_dsijqcd[0][0]*p_kpterms_qcd->Coupling()*p_loopme->ME_AssContribs(i);

  if (m_checkpoles)  CheckPoleCancelation(mom);
  if (m_checkfinite) CheckFinite(I,V);
  if (m_checkborn)   CheckBorn();

  m_lastbxs=B;
  m_lastv=V;
  m_lasti=I;
  if (p_loopme)
    for (size_t i(0);i<p_loopme->ME_AssContribs_Size();++i)
      m_wass[i] *= m_lastk;
  double M2(B+V+I);


  if (IsBad(M2) || msg_LevelIsDebugging()) {
    size_t precision(msg->Out().precision());
    msg->SetPrecision(16);
    msg_Error()<<METHOD<<"("<<Name()<<"){\n  M2 = "<<M2
                       <<"\n  m_eta0 = "<<m_eta0<<" ,  m_x0 = "<<m_x0
                       <<" ,  m_eta1 = "<<m_eta1<<" ,  m_x1 = "<<m_x1<<"\n"
                       <<" ,  \\alpha_S = "<<m_cpls.Get("Alpha_QCD")->Default()
                                             *m_cpls.Get("Alpha_QCD")->Factor()
                       <<" ,  \\alpha = "<<m_cpls.Get("Alpha_QED")->Default()
                                           *m_cpls.Get("Alpha_QED")->Factor()
                       <<"\n"
                       <<" ,  m_dsijqcd[0][0] = "
                       <<(m_dsijqcd.size()?m_dsijqcd[0][0]:0.)
                       <<" ,  m_dsijew[0][0] = "
                       <<(m_dsijew.size()?m_dsijew[0][0]:0.)
                       <<"\n  B = "<<B<<" ,  V = "<<V
	               <<" ,  I = "<<I<<"\n  V+I = "<<V+I
                       <<" ,  \\delta = "<<(V+I)/B
                       <<"\n  K(B) = "<<kfactorb
                       <<" ,  K(VI) = "<<kfactorvi
                       <<"\n  norm = "<<m_Norm<<std::endl;
    for (size_t i=0;i<m_nin+m_nout;i++)
      msg_Error()<<"  "<<i<<": "<<mom[i]<<std::endl;
    msg_Error()<<"}\n";
    msg->SetPrecision(precision);
  }

  return M2;
}

void Single_Virtual_Correction::FillAmplitudes
(vector<METOOLS::Spin_Amplitudes>& amps,vector<vector<Complex> >& cols)
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

void Single_Virtual_Correction::SetNLOMC(PDF::NLOMC_Base *const mc)
{
  if (p_LO_process) p_LO_process->SetNLOMC(mc);
  if (p_kpterms_qcd) p_kpterms_qcd->SetNLOMC(mc);
  if (p_kpterms_ew) p_kpterms_ew->SetNLOMC(mc);
  p_nlomc=mc;
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
  wgtinfo.m_swap=p_int->Momenta()[0][3]<p_int->Momenta()[1][3];
  wgtinfo.m_y1=wgtinfo.m_swap?m_x1:m_x0;
  wgtinfo.m_y2=wgtinfo.m_swap?m_x0:m_x1;
  if (wgtinfo.m_type&mewgttype::VI)
    for (size_t i=0;i<2;i++) wgtinfo.m_wren[i]=m_cmur[i]*=m_lastki;
  if (p_kpterms_qcd) p_kpterms_qcd->FillMEwgts(wgtinfo);
  if (p_kpterms_ew)  p_kpterms_ew->FillMEwgts(wgtinfo);
  for (size_t i=2;i<wgtinfo.m_wren.size();++i) wgtinfo.m_wren[i]*=m_lastk;
  for (size_t i=0;i<wgtinfo.m_wfac.size();++i) wgtinfo.m_wfac[i]*=m_lastk;
}

void Single_Virtual_Correction::MPICollect(std::vector<double> &sv,size_t &i)
{
  sv.resize(i+4);
  sv[i+0]=m_mn;
  sv[i+1]=m_mbsum;
  sv[i+2]=m_mvsum;
  sv[i+3]=m_misum;
  i+=4;
}

void Single_Virtual_Correction::MPIReturn(std::vector<double> &sv,size_t &i)
{
  m_mn=sv[i+0];
  m_mbsum=sv[i+1];
  m_mvsum=sv[i+2];
  m_misum=sv[i+3];
  i+=4;
}

void Single_Virtual_Correction::MPISync(const int mode)
{
  Process_Base::MPISync(mode);
#ifdef USING__MPI
  m_n+=m_mn;
  m_bsum+=m_mbsum;
  m_vsum+=m_mvsum;
  m_isum+=m_misum;
  m_mn=m_mbsum=m_mvsum=m_misum=0.0;
#endif
}

Flavour Single_Virtual_Correction::ReMap(const Flavour &fl,
                                         const size_t &id) const
{
  return p_LO_process->ReMap(fl,id);
}

void Single_Virtual_Correction::AttachChargeFactors()
{
  if (m_stype&sbt::qed) {
    for (size_t i(0);i<m_dsijew.size();i++) {
      for (size_t k(0);k<m_dsijew[i].size();k++) {
        m_dsijew[i][k]*=m_Q2ij[i][k];
      }
    }
  }
}

void Single_Virtual_Correction::ComputeChargeFactors()
{
  size_t nijqed(p_LO_process->PartonListQED().size());
  std::vector<size_t> nPFFsplittings(nijqed,0);
  for (size_t i=0;i<nijqed;i++) {
    for (size_t j=0;j<nijqed;j++) {
      if (i==j) m_Q2ij[i][j]=1.;
      else if (m_flavs[p_LO_process->PartonListQED()[i]].IsPhoton()) {
        if (AllowAsSpecInPFF(i,j)) {
          m_Q2ij[i][j]=-1.;
          nPFFsplittings[i]++;
        }
      }
      else if (m_flavs[p_LO_process->PartonListQED()[i]].Charge() &&
               m_flavs[p_LO_process->PartonListQED()[i]].Charge()) {
        double Qi(m_flavs[p_LO_process->PartonListQED()[i]].Charge());
        double Qj(m_flavs[p_LO_process->PartonListQED()[j]].Charge());
        bool inii(p_LO_process->PartonListQED()[i]<m_nin);
        bool inij(p_LO_process->PartonListQED()[j]<m_nin);
        m_Q2ij[i][j]=(inii==inij?1.:-1)*Qi*Qj;
      }
    }
  }
  for (size_t i=0;i<nijqed;i++) {
    for (size_t j=0;j<nijqed;j++) {
      if (i==j) continue;
      if (m_flavs[p_LO_process->PartonListQED()[i]].IsPhoton()) {
        if (nPFFsplittings[i]==0) {
          if (p_LO_process->PartonListQED()[i]<m_nin) {
            THROW(fatal_error,"No spectator for Photon splitting assigned. "
                              +std::string("Try different DIPOLES:PFF_IS_RECOIL_SCHEME."));
          }
          else {
            THROW(fatal_error,"No spectator for Photon splitting assigned. "
                              +std::string("Try different DIPOLES:PFF_FS_RECOIL_SCHEME."));
          }
        }
        m_Q2ij[i][j]/=(double)nPFFsplittings[i];
      }
    }
  }
  if (msg_LevelIsDebugging()) {
    if (m_Q2ij.size()) {
      msg_Out()<<std::setw(4)<<"Charge factors";
      for (size_t j(0);j<m_Q2ij[0].size();++j) {
        if (j==0) msg_Out()<<std::setw(11)<<"j="<<std::setw(2)<<j;
        else      msg_Out()<<std::setw(21)<<"j="<<std::setw(2)<<j;
      }
      msg_Out()<<std::endl;
      for (size_t i(0);i<m_Q2ij.size();++i) {
        msg_Out()<<std::setw(2)<<"i="<<std::setw(2)<<i;
        for (size_t j(0);j<m_Q2ij[i].size();++j) {
          msg_Out()<<std::setw(23)<<m_Q2ij[i][j];
        }
        msg_Out()<<std::endl;
      }
    }
  }
}

ist::itype Single_Virtual_Correction::AssignType(const size_t& id,
                                                 const ATOOLS::sbt::subtype st)
{
  if (st==sbt::qcd) {
    // QCD subtraction
    if      (m_flavs[id].IsQuark() && !m_flavs[id].IsMassive())  return ist::q;
    else if (m_flavs[id].IsQuark() &&  m_flavs[id].IsMassive())  return ist::Q;
    else if (m_flavs[id].IsGluon())                              return ist::g;
    else if (IsSusy(m_flavs[id]) && m_flavs[id].IsScalar())     return ist::sQ;
    else if (IsSusy(m_flavs[id]) && m_flavs[id].IsFermion())    return ist::sG;
  }
  else if (st==sbt::qed) {
    // QED subtraction
    if      (m_flavs[id].IsQuark() && !m_flavs[id].IsMassive())  return ist::q;
    else if (m_flavs[id].IsQuark() &&  m_flavs[id].IsMassive())  return ist::Q;
    else if (m_flavs[id].IsLepton() && !m_flavs[id].IsMassive()) return ist::q;
    else if (m_flavs[id].IsLepton() &&  m_flavs[id].IsMassive()) return ist::Q;
    else if (m_flavs[id].IsPhoton())                             return ist::g;
    else if (m_flavs[id].IsVector() && m_flavs[id].IsMassive())  return ist::V;
    else if (IsSusy(m_flavs[id]) && m_flavs[id].IsScalar())     return ist::sQ;
    else if (IsSusy(m_flavs[id]) && m_flavs[id].IsFermion())    return ist::sG;
  }
  else THROW(not_implemented,"Cannot assign type for "+ToString(st)
                             +" subtraction.");
  return ist::none;
}

void Single_Virtual_Correction::FillProcessMap(NLOTypeStringProcessMap_Map *apmap)
{
  Process_Base::FillProcessMap(apmap);
  p_LO_process->SetProcMap(apmap);
}

void Single_Virtual_Correction::PrintDSij()
{
  if (msg_LevelIsDebugging()) {
    DEBUG_FUNC("");
    double precision(msg_Out().precision());
    msg->SetPrecision(16);
    if (m_dsijqcd.size()) {
      msg_Out()<<std::setw(4)<<"QCD";
      for (size_t j(0);j<m_dsijqcd[0].size();++j) {
        msg_Out()<<std::setw(21)<<"j="<<std::setw(2)<<j;
      }
      msg_Out()<<std::endl;
      for (size_t i(0);i<m_dsijqcd.size();++i) {
        msg_Out()<<std::setw(2)<<"i="<<std::setw(2)<<i;
        for (size_t j(0);j<m_dsijqcd[i].size();++j) {
          msg_Out()<<std::setw(23)<<m_dsijqcd[i][j];
        }
        msg_Out()<<std::endl;
      }
    }
    if (m_dsijew.size()) {
      msg_Out()<<std::setw(4)<<"EW";
      for (size_t j(0);j<m_dsijew[0].size();++j) {
        msg_Out()<<std::setw(21)<<"j="<<std::setw(2)<<j;
      }
      msg_Out()<<std::endl;
      for (size_t i(0);i<m_dsijew.size();++i) {
        msg_Out()<<std::setw(2)<<"i="<<std::setw(2)<<i;
        for (size_t j(0);j<m_dsijew[i].size();++j) {
          msg_Out()<<std::setw(23)<<m_dsijew[i][j];
        }
        msg_Out()<<std::endl;
      }
    }
    msg->SetPrecision(precision);
  }
}

bool Single_Virtual_Correction::AllowAsSpecInPFF(const size_t& i,
                                                 const size_t &k)
{
  if (p_LO_process->PartonListQED()[i]<m_nin) {
    switch (m_pspisrecscheme) {
    case 0:
      if (p_LO_process->PartonListQED()[k]<m_nin) return true;
      break;
    case 1:
      if (p_LO_process->PartonListQED()[k]>=m_nin) return true;
      break;
    case 2:
      if (m_flavs[p_LO_process->PartonListQED()[k]].Charge()) return true;
      break;
    case 3:
      if (!m_flavs[p_LO_process->PartonListQED()[k]].Charge()) return true;
      break;
    case 4:
      return true;
      break;
    default:
      THROW(fatal_error,"Unknown IS P->ff recoil scheme.")
      break;
    }
    return false;
  }
  else {
    switch (m_pspfsrecscheme) {
    case 0:
      if (p_LO_process->PartonListQED()[k]<m_nin) return true;
      break;
    case 1:
      if (p_LO_process->PartonListQED()[k]>=m_nin) return true;
      break;
    case 2:
      if (m_flavs[p_LO_process->PartonListQED()[k]].Charge()) return true;
      break;
    case 3:
      if (!m_flavs[p_LO_process->PartonListQED()[k]].Charge()) return true;
      break;
    case 4:
      return true;
      break;
    default:
      THROW(fatal_error,"Unknown FS P->ff recoil scheme.")
      break;
    }
    return false;
  }
  return false;
}

void Single_Virtual_Correction::SetCaller(PHASIC::Process_Base *const proc)
{
  p_caller=proc;
  p_LO_process->SetCaller(static_cast<Single_Virtual_Correction*>(proc)->p_LO_process);
}
