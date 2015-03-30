#include "AMEGIC++/DipoleSubtraction/Single_OSTerm.H"
#include "AMEGIC++/DipoleSubtraction/Single_LOProcess.H"
#include "AMEGIC++/DipoleSubtraction/Single_LOProcess_MHV.H"

#include "AMEGIC++/Main/Single_Process.H"
#include "AMEGIC++/Main/Single_Process_MHV.H"
#include "AMEGIC++/Main/Single_Process_External.H"

#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PDF/Main/ISR_Handler.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"

#include "AMEGIC++/DipoleSubtraction/FF_DipoleSplitting.H"
#include "AMEGIC++/DipoleSubtraction/FI_DipoleSplitting.H"
#include "AMEGIC++/DipoleSubtraction/IF_DipoleSplitting.H"
#include "AMEGIC++/DipoleSubtraction/II_DipoleSplitting.H"

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

Single_OSTerm::Single_OSTerm(const Process_Info &pinfo,size_t pi,size_t pj,size_t swit, Process_Integrator* pint) :
  p_partner(this), p_OS_mom(0), p_os_process(0), m_wwindow(5.), p_realint(pint)
{
  DEBUG_FUNC("");
  PHASIC::Process_Base::Init(pinfo, pint->Beam(), pint->ISR());
  AMEGIC::Process_Base::Init();

  m_pi = pi;
  m_pj = pj;
  m_switch = swit;
  m_pk=2;
  if (m_pi==m_pk || m_pj==m_pk) m_pk ++;
  if (m_pi==m_pk || m_pj==m_pk) m_pk ++;

  m_name+= "_OS"+ToString(m_pi)+"_"+ToString(m_pj)+"_"+ToString(swit);

  bool val=DetermineType();
  if (!val) return;

 //Construct OS Process
  Process_Info lopi(m_pinfo);

  lopi.m_fi.m_nloqcdtype=nlo_type::lo;
  lopi.m_fi.m_nloewtype=nlo_type::lo;
  lopi.m_fi.m_ps.erase(FindInInfo(lopi.m_fi, m_pi-m_nin));
  if (m_pj>m_pi) lopi.m_fi.m_ps.erase(FindInInfo(lopi.m_fi, m_pj-m_nin-1));
  else lopi.m_fi.m_ps.erase(FindInInfo(lopi.m_fi, m_pj-m_nin));
  vector<Subprocess_Info>::iterator part=FindInInfo(m_pinfo.m_fi, m_pi-m_nin);
  lopi.m_fi.m_ps.push_back(Subprocess_Info(m_flij,part->m_id, part->m_pol));
  lopi.m_fi.m_ps.back().m_id="osdecay";

  Subprocess_Info ACFS;
  ACFS.m_ps.push_back(lopi.m_fi.m_ps.back());
  DEBUG_VAR(ACFS);
  BuildDecay(ACFS);
  DEBUG_VAR(ACFS);
  for (size_t afsi(0);afsi<ACFS.m_ps.size();++afsi) {
    msg_Debugging()<<METHOD<<"(): Check N_max ("
		     <<"): {\n"<<ACFS.m_ps[afsi]<<"}\n";
    lopi.m_fi.m_ps.back().GetNMax(ACFS.m_ps[afsi]);
  }
  m_osinfo= lopi;

  double helpd(5.);
  Data_Reader reader(" ",";","!","=");
  reader.AddComment("#");
  reader.SetInputPath(rpa->GetPath());
  reader.SetInputFile(rpa->gen.Variable("ME_DATA_FILE"));
  if (reader.ReadFromFile(helpd,"OS_SUB_WINDOW")) {
    m_wwindow = helpd;
    msg_Tracking()<<"Set width window for os subtraction="<<m_wwindow<<"."<<std::endl;
  }
}


Single_OSTerm::~Single_OSTerm()
{
  p_scale=NULL;
  if (p_os_process) {delete p_os_process; p_os_process=0;}
  if (p_OS_mom)     {delete[] p_OS_mom; p_OS_mom=0;}
}

/*------------------------------------------------------------------------------
  
  Generic stuff for initialization of Single_OSTerms
      
  ------------------------------------------------------------------------------*/

void Single_OSTerm::BuildDecay
(Subprocess_Info &ACFS)
{
  int osf=1; 
  Subprocess_Info ACDIS, ACDFS;
  std::string decid("osdecay");
  ACDIS.m_ps.resize(1);
  ACDIS.m_ps.front().m_ps.clear();
  for (size_t i(0);i<ACDIS.m_ps.size();++i)
    ACDIS.m_ps[i].m_ps.push_back(Subprocess_Info(m_flij,decid,""));
  ACDFS.m_ps.resize(1);
  ACDFS.m_ps.front().m_ps.clear();
  for (size_t i(0);i<ACDFS.m_ps.size();++i)
	ACDFS.m_ps[i].m_ps.push_back(Subprocess_Info(m_fli,"",""));
  for (size_t i(0);i<ACDFS.m_ps.size();++i)
	ACDFS.m_ps[i].m_ps.push_back(Subprocess_Info(m_flj,"",""));
  DEBUG_VAR(ACDFS);
  if (ACDIS.m_ps.empty() || ACDFS.m_ps.empty())
    THROW(fatal_error,"Wrong decay specification");
  Subprocess_Info &CIS(ACDIS.m_ps.front());
  size_t oldsize(ACFS.m_ps.size()), cdsize(ACDFS.m_ps.size());
  ACFS.m_ps.resize(oldsize*cdsize);
  for (size_t cfss(1);cfss<cdsize;++cfss)
    for (size_t acfsi(0);acfsi<oldsize;++acfsi)
      ACFS.m_ps[cfss*oldsize+acfsi]=
	ACFS.m_ps[(cfss-1)*oldsize+acfsi];
  for (size_t acfsi(0);acfsi<oldsize;++acfsi) {
    for (size_t cfss(0);cfss<cdsize;++cfss) {
      Subprocess_Info &CFS(ACDFS.m_ps[cfss]);
      DEBUG_VAR(CFS);
      msg_Debugging()<<METHOD<<"(): Init decay {\n"<<CIS<<CFS<<"}\n";
      if (CIS.NExternal()!=1 || CFS.NExternal()<2)
        THROW(fatal_error,"Wrong number of particles in decay");
      if (!ACFS.m_ps[cfss*oldsize+acfsi].AddDecay(CIS,CFS,osf))
        cout<<"decay not matching..."<<endl;
      }
    }
}


bool Single_OSTerm::DetermineType() {
  m_valid=true;
  if (m_pj<m_nin) m_valid=false;
  if (!m_valid) return false;
  m_fli=m_flavs[m_pi];
  m_flj=m_flavs[m_pj];
  m_flk=m_flavs[m_pk];
  if (!m_flj.IsQuark()) return m_valid=false;
  int kfj = m_flj.Kfcode();
  bool anti = m_flj.IsAnti();
  if (IsGluino(m_fli) || IsNeutralino(m_fli)) {
    if (m_switch==0) {
      m_flij = Flavour((kf_code)(1000000+kfj));
    }
    else if (m_switch==1) {
      m_flij = Flavour((kf_code)(2000000+kfj));
    }
    else return m_valid = false;
    if (anti) m_flij = m_flij.Bar();
  }
  else if (IsSquark(m_fli)){
    if (((m_fli.Kfcode() - 1000000)==m_flj.Kfcode()) ||
         ((m_fli.Kfcode() - 2000000)==m_flj.Kfcode())) {
        if ((m_fli.IsAnti() && m_flj.IsAnti()) || 
             (!m_fli.IsAnti() && !m_flj.IsAnti())) return m_valid = false;
        if (m_switch==0) m_flij = Flavour(1000021);
        if (m_switch==1) m_flij = Flavour(1000022);
        if (m_switch==2) m_flij = Flavour(1000023);
        if (m_switch==3) m_flij = Flavour(1000025);
        if (m_switch==4) m_flij = Flavour(1000035);
      } 
    else return m_valid = false;
  }
  else return m_valid = false;
  if (m_flij.Mass()<(m_fli.Mass()+m_flj.Mass())) return m_valid=false;
  return m_valid;
}

vector<Subprocess_Info>::iterator
Single_OSTerm::FindInInfo(Subprocess_Info& fi, int idx) const
{
  // Find particle #idx in the given final state tree
  // and make sure it is not part of a decay for the time being
  int cnt=0;
  for (size_t i=0; i<fi.m_ps.size(); ++i) {
    cnt+=fi.m_ps[i].NExternal();
    if (idx<cnt) {
      if (fi.m_ps[i].NExternal()==1) {
        return fi.m_ps.begin()+i;
      }
      else {
        THROW(not_implemented,
              "Dipole subtraction for coloured particles in decays not implemented yet.");
      }
    }
  }
  THROW(fatal_error, "Internal Error");
  return fi.m_ps.end();
}

double Single_OSTerm::lambda(double x, double y, double z)
{
  return x*x +y*y +z*z - 2.*x*y - 2.*y*z - 2.*z*x;
}

void Single_OSTerm::SetLOMomenta(const Vec4D* moms,const ATOOLS::Poincare &cms)
{
  PHASIC::Process_Info osinfo(p_os_process->Info());
  int ip(-1), jp(-1), kp(-1),decayfound(0);
  for (size_t i=0; i<osinfo.m_fi.m_ps.size(); i++) if (osinfo.m_fi.m_ps[i].m_id == "osdecay") {
    if (osinfo.m_fi.m_ps[i].m_ps[0].m_fl==m_fli){
      ip = m_nin+i;
      jp = ip+1;
    }
    else {
      jp = m_nin+i;
      ip = jp+1;
    }
    kp=m_pk;
    if (i+m_nin<=m_pk && m_pk<m_pi) kp++;
    if (i+m_nin<=m_pk && m_pk<m_pi) kp++;
    decayfound=1;
  }
  if (decayfound==0) THROW(fatal_error,"os decay not found")

  int cnt=0;
  for (size_t i=0;i<m_nin+m_nout;i++) {
    for (;cnt==ip||cnt==jp||cnt==kp;) cnt++;
    if (i!=m_pi&&i!=m_pj&&i!=m_pk) {
      Vec4D momtemp(moms[i]);
      p_OS_mom[cnt] =momtemp;
      cnt++;
    }
  }
  Vec4D Q(moms[m_pi]+moms[m_pj]+moms[m_pk]);
  double mij2 = sqr(m_flij.Mass());
  double mk2 = sqr(m_flk.Mass());
  double q2 = Q*Q;
  if (sqrt(q2)<(m_flij.Mass()+m_flk.Mass())) THROW(fatal_error,"q2 too small!!!!");
  double pipj2 = (moms[m_pi]+moms[m_pj])*(moms[m_pi]+moms[m_pj]);
  Vec4D ptk = sqrt(lambda(q2,mij2,mk2)/lambda(q2,pipj2,mk2))*(moms[m_pk] - 
                     (Q*moms[m_pk])/q2*Q) +(q2+mk2-mij2)/(2.*q2)*Q;

  Vec4D ptij = Q-ptk;

  Vec4D pi = moms[m_pi];
  Vec4D pj = moms[m_pj];
  Poincare bstcms(ptij);
  bstcms.Boost(pi);
  bstcms.Boost(pj);
  bstcms.Boost(ptij);
  double mi2(sqr(m_fli.Mass())),mj2(sqr(m_flj.Mass()));
  double p2 = pi.PSpat2();
  double k2 = (mij2*mij2 + sqr(mi2-mj2) -2.*mij2*(mi2+mj2))/(4.*mij2*p2);
  double k= sqrt(k2);
  for (size_t i=1; i<4;i++) pi[i] *=k;
  pi[0]= sqrt(mi2+pi.PSpat2());
  
  pj = ptij -pi;
  bstcms.BoostBack(pi);
  bstcms.BoostBack(pj);
  Vec4D pti = pi;
  Vec4D ptj = pj;
  
  p_OS_mom[ip] = pi;
  p_OS_mom[jp] = pj;
  p_OS_mom[kp] = ptk;

  for (size_t i=0;i<m_nin+m_nout;++i) {
    p_OS_labmom[i] = p_OS_mom[i];
    cms.BoostBack(p_OS_labmom[i]);
  }
}

void Single_OSTerm::PrintLOmom()
{
  if (this!=p_partner) { p_partner->PrintLOmom();return;}
  for (size_t i=0;i<m_nin+m_nout-1;i++) cout<<i<<": "<<p_OS_mom[i]<<endl;
}

/*------------------------------------------------------------------------------

  Initializing libraries, amplitudes, etc.

  ------------------------------------------------------------------------------*/



int Single_OSTerm::InitAmplitude(Amegic_Model *model,Topology* top,
				    vector<Process_Base *> & links,
				    vector<Process_Base *> & errs)
{
  p_OS_mom = new Vec4D[m_nin+m_nout];
  p_OS_labmom.resize(m_nin+m_nout);

  if (m_osinfo.m_amegicmhv>0) {
    if (CF.MHVCalculable(m_osinfo)) p_os_process = new Single_Process_MHV();
    if (m_osinfo.m_amegicmhv==2) {m_valid=0; return 0;}
  }
  if (!p_os_process) p_os_process = new AMEGIC::Single_Process();

  int status;
  p_os_process->PHASIC::Process_Base::Init(m_osinfo,p_realint->Beam(),p_realint->ISR());
  Poincare cms;
  SetLOMomenta(p_testmoms,cms);
  p_os_process->SetTestMoms(p_OS_mom);
  status=p_os_process->InitAmplitude(model,top,links,errs);
  if (status<=0) { 
    m_valid=0;
    return status;
  }

  m_subevt.m_n    = m_nin+m_nout;
  m_subevt.p_fl   = &(p_os_process->Flavours().front());
  m_subevt.p_dec  = &m_decins;
  m_subevt.p_mom  = &p_OS_labmom.front();
  m_subevt.m_i    = m_pi;
  m_subevt.m_j    = m_pj;
  m_subevt.m_k    = m_pk;
  m_subevt.p_proc = p_os_process->Integrator();
  m_sids.resize(m_nin+m_nout);
  for (size_t i(0);i<m_nin+m_nout;++i) m_sids[i]=1<<i;
  m_subevt.p_id=&m_sids.front();
  m_subevt.m_pname = GenerateName(p_os_process->Info().m_ii,p_os_process->Info().m_fi);

  SetMaxOrders(p_os_process->MaxOrders());
  SetMinOrders(p_os_process->MinOrders());
  SetMaxOrder(0,p_os_process->MaxOrder(0)+1);
  SetMinOrder(0,p_os_process->MinOrder(0)+1);
  return 1;
}




/*------------------------------------------------------------------------------
  
  Phase space initialization
  
  ------------------------------------------------------------------------------*/


bool Single_OSTerm::SetUpIntegrator() 
{  
  bool res=p_os_process->SetUpIntegrator();
  if (res) return res;
  return true;
}


/*------------------------------------------------------------------------------
  
  Process management
  
  ------------------------------------------------------------------------------*/
void Single_OSTerm::SetLookUp(const bool lookup)
{
  m_lookup=lookup; 
  if (p_os_process) p_os_process->SetLookUp(lookup);
}

void Single_OSTerm::Minimize()
{
  if (p_partner==this) return;
  if (p_os_process) 
    {delete p_os_process; p_os_process=0;}
  if (p_OS_mom)     {delete[] p_OS_mom; p_OS_mom=0;}
  m_subevt.p_mom = p_partner->GetSubevt()->p_mom;
}


double Single_OSTerm::Partonic(const Vec4D_Vector &_moms,const int mode) { return 0.; }

double Single_OSTerm::operator()(const ATOOLS::Vec4D * mom,const ATOOLS::Poincare &cms,const int mode)
{
  ATOOLS::Vec4D dec(mom[m_pi]+mom[m_pj]);
  double invm2 = dec*dec;
  if (abs(sqr(m_flij.Mass()) - invm2) > m_wwindow*m_flij.Mass()*m_flij.Width() ) return 0.;
  if (sqrt((dec+mom[m_pk]).Abs2())<(m_flij.Mass()+m_flk.Mass())) return 0.;
  ResetLastXS();
  p_os_process->ResetLastXS();
  SetLOMomenta(mom,cms);

  bool trg(false);
  trg= p_os_process->Trigger(p_OS_labmom) || !p_os_process->Selector()->On();
  p_int->SetMomenta(p_OS_labmom);
  p_os_process->Integrator()->SetMomenta(p_OS_labmom);

  m_subevt.m_me = m_subevt.m_mewgt = m_subevt.m_result = 0.;

  if (!trg) return m_lastxs=m_subevt.m_me=m_subevt.m_mewgt=0.;

  ATOOLS::Vec4D_Vector lomoms;
  for (size_t i=0;i<m_nin+m_nout;i++) lomoms.push_back(p_OS_mom[i]);
  p_os_process->ScaleSetter()->CalculateScale(lomoms);
  double norm = p_os_process->Norm();
  double M2 =p_os_process->operator()(&lomoms.front())*norm;

  double mij2 = sqr(m_flij.Mass());
  double wij2 = sqr(m_flij.Width());
  if (wij2==0.) THROW (fatal_error,"width is zero for on shell decay");
  double factor= mij2*wij2/(sqr(invm2-mij2) + mij2*wij2);
  m_lastxs = factor*M2;
  m_subevt.m_me = m_subevt.m_mewgt = -m_lastxs;
  m_subevt.m_mu2[stp::fac] = p_scale->Scale(stp::fac);
  m_subevt.m_mu2[stp::ren] = p_scale->Scale(stp::ren);

  return m_lastxs;
}

void Single_OSTerm::SetSelector(const PHASIC::Selector_Key &key)
{
  p_os_process->SetSelector(key);
}

int Single_OSTerm::NumberOfDiagrams() { 
  if (p_partner==this) return p_os_process->NumberOfDiagrams(); 
  return p_partner->NumberOfDiagrams();
}

Point * Single_OSTerm::Diagram(int i) { 
  if (p_partner==this) return p_os_process->Diagram(i); 
  return p_partner->Diagram(i);
} 

void Single_OSTerm::AddChannels(std::list<std::string>*psln)
{
  p_os_process->AddChannels(psln);
}

void Single_OSTerm::PrintProcessSummary(int it)
{
  for(int i=0;i<it;i++) cout<<"  ";
  cout<<m_pi<<"-"<<m_pj<<"-"<<m_pk<<" ("<<p_os_process->Name()<<")";
  if (p_partner!=this) {
    cout<<"; partner (*"<<m_sfactor<<"): ";
    p_partner->PrintProcessSummary(0);
    return;
  }
  cout<<endl;
} 

void Single_OSTerm::SetScale(const Scale_Setter_Arguments &args)
{
  if (!p_os_process->IsMapped()) p_os_process->SetScale(args);
  p_scale=p_os_process->Partner()->ScaleSetter();
}
