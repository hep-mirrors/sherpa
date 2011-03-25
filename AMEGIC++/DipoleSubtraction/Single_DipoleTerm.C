#include "AMEGIC++/DipoleSubtraction/Single_DipoleTerm.H"
#include "AMEGIC++/DipoleSubtraction/Single_LOProcess.H"
#include "AMEGIC++/DipoleSubtraction/Single_LOProcess_MHV.H"

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

Single_DipoleTerm::Single_DipoleTerm(const Process_Info &pinfo,size_t pi,size_t pj,size_t pk, Process_Integrator* pint) :
  m_maxgsmass(0.), p_partner(this), p_LO_process(0), p_LO_mom(0), m_ftype(0), p_dipole(0), p_realint(pint)
{
  DEBUG_FUNC("");
  PHASIC::Process_Base::Init(pinfo, pint->Beam(), pint->ISR());
  AMEGIC::Process_Base::Init();

  m_pi = pi;
  m_pj = pj;
  m_pk = pk;

  m_name+= "_RS"+ToString(m_pi)+"_"+ToString(m_pj)+"_"+ToString(m_pk);

  bool val=DetermineType();
  if (!val) return;

  m_LOpij = m_nin+m_nout-2;
  if (m_pi<m_nin) m_LOpij = m_pi;
  m_LOpk  = m_pk;
  if (m_pi<m_pk&&m_pi>=m_nin) m_LOpk--;
  if (m_pj<m_pk) m_LOpk--;

 //Construct LO Process
  Process_Info lopi=m_pinfo;
  if (m_pi<m_nin) {
    lopi.m_ii.m_ps[m_pi]=Subprocess_Info(m_flij,m_pinfo.m_ii.m_ps[m_pi].m_id,
					 m_pinfo.m_ii.m_ps[m_pi].m_pol);
    lopi.m_ii.m_ps[m_pi].m_tag=1;
    lopi.m_fi.m_ps.erase(FindInInfo(lopi.m_fi, m_pj-m_nin));
  }
  else {
    lopi.m_fi.m_ps.erase(FindInInfo(lopi.m_fi, m_pi-m_nin));
    lopi.m_fi.m_ps.erase(FindInInfo(lopi.m_fi, m_pj-m_nin-1));
    vector<Subprocess_Info>::iterator part=FindInInfo(m_pinfo.m_fi, m_pi-m_nin);
    lopi.m_fi.m_ps.push_back(Subprocess_Info(m_flij,part->m_id, part->m_pol));
    lopi.m_fi.m_ps.back().m_tag=1;
  }
  if (m_LOpk<m_nin) lopi.m_ii.m_ps[m_LOpk].m_tag=2;
  else FindInInfo(lopi.m_fi, m_LOpk-m_nin)->m_tag=2;
  DEBUG_VAR(lopi);

  if (lopi.m_amegicmhv>0) {
    if (CF.MHVCalculable(lopi))
      p_LO_process = new Single_LOProcess_MHV(lopi, p_int->Beam(), p_int->ISR());
    if (lopi.m_amegicmhv==2) { m_valid=0; return; }
  }
  if (!p_LO_process)
    p_LO_process = new Single_LOProcess(lopi, p_int->Beam(), p_int->ISR());
  if (!p_LO_process) THROW(fatal_error,"LO process unknown");

  p_LO_mom = new Vec4D[m_nin+m_nout-1];
  p_LO_labmom.resize(m_nin+m_nout-1); 
  p_LO_process->SetTestMoms(p_LO_mom);

  m_subevt.m_n    = m_nin+m_nout-1;
  m_subevt.p_fl   = &(p_LO_process->Flavours().front());
  m_subevt.p_dec  = &m_decins;
  m_subevt.p_mom  = &p_LO_labmom.front();
  m_subevt.m_i    = m_pi;
  m_subevt.m_j    = m_pj;
  m_subevt.m_k    = m_pk;
  m_subevt.p_proc = p_LO_process->Integrator();

  m_sids.resize(m_nin+m_nout-1);
  size_t em=p_LO_process->GetEmit();
  size_t sp=p_LO_process->GetSpect();
  for (size_t cnt=0, i=0;i<m_nin+m_nout;i++) {
    for (;cnt==em||cnt==sp;) cnt++;
    if (i!=m_pi&&i!=m_pj&&i!=m_pk) {
      m_sids[cnt] = 1<<i;
      cnt++;
    }
  }
  m_sids[em]=(1<<m_pi)|(1<<m_pj);
  m_sids[sp]=1<<m_pk;
  m_subevt.m_ijt=em;
  m_subevt.m_kt=sp;
  m_subevt.p_id=&m_sids.front();
  m_subevt.m_pname=GenerateName(p_LO_process->Info().m_ii,p_LO_process->Info().m_fi);
  m_subevt.m_pname=m_subevt.m_pname.substr(0,m_subevt.m_pname.rfind("__"));

  m_dalpha = 1.;
  double helpd;
  Data_Reader reader(" ",";","!","=");
  reader.AddComment("#");
  reader.SetInputPath(rpa.GetPath());
  reader.SetInputFile(rpa.gen.Variable("ME_DATA_FILE"));
  if (reader.ReadFromFile(helpd,"DIPOLE_ALPHA")) {
    m_dalpha = helpd;
    msg_Tracking()<<"Set dipole cut alpha="<<m_dalpha<<"."<<std::endl;
  }

  int helpi;
  if (reader.ReadFromFile(helpi,"DIPOLE_NF_GSPLIT")) {
    msg_Tracking()<<"Set number of flavours from gluon splitting="<<helpi<<"."<<std::endl;
    Flavour flav((kf_code)(helpi));
    m_maxgsmass=flav.Mass();
  }
}


Single_DipoleTerm::~Single_DipoleTerm()
{
  p_scale=NULL;
  if (p_LO_process) {delete p_LO_process; p_LO_process=0;}
  if (p_LO_mom)     {delete[] p_LO_mom; p_LO_mom=0;}
  if (p_dipole)     {delete p_dipole; p_dipole=0;}
}

/*------------------------------------------------------------------------------
  
  Generic stuff for initialization of Single_DipoleTermes
      
  ------------------------------------------------------------------------------*/
bool Single_DipoleTerm::DetermineType() {
  m_valid=true;
  if (m_pi>=m_pj) m_valid=false;
  if (m_pj<m_nin) m_valid=false;
  if (!m_valid) return false;
  bool massive=0;
  bool massiveini=0;
  m_fli=m_flavs[m_pi];
  m_flj=m_flavs[m_pj];
  m_flk=m_flavs[m_pk];
  if (m_flk.IsMassive()) massive=1;
  if (massive&&m_pk<m_nin) massiveini=1; 
  if (m_pi>=m_nin) {
    if (m_fli.IsMassive()||m_flj.IsMassive()) massive=1;
    if (!massive) {
      if (m_pk>=m_nin) m_dipoletype = dpt::f_f;
      else m_dipoletype = dpt::f_i;
    }
    else {
      if (m_pk>=m_nin) m_dipoletype = dpt::f_fm;
      else m_dipoletype = dpt::f_im;
    }
  }
  else {
    if (m_fli.IsMassive()) massiveini=1;
    if (massiveini==0 && m_flj.IsMassive()) {
      m_valid=false;
      return m_valid;
    }
    if (!massive) {
      if (m_pk>=m_nin) m_dipoletype = dpt::i_f;
      else m_dipoletype = dpt::i_i;
    }
    else {
      if (m_pk>=m_nin) m_dipoletype = dpt::i_fm;
    }
  }

  if (massiveini) {
    msg_Error()<<METHOD<<" Cannot handle massive initial state! Abort."<<endl;
    abort();
  }

  switch (m_dipoletype) {
  case dpt::f_f:
  case dpt::f_fm:
  case dpt::f_i:
  case dpt::f_im:
    if (m_fli==Flavour(kf_gluon)) {
      m_flij = m_flj;
      if (m_flj==m_fli) m_ftype = 4;
      else if (!m_fli.IsSusy()) m_ftype = 2;
      else if (m_fli.IsGluino()) m_ftype = 6;
      else if (m_fli.IsSquark()) m_ftype = 8;
      else THROW(fatal_error,"SUSY particle in dipole term, but not squark or gluino");
    }
    else if (m_flj==Flavour(kf_gluon)) {
      m_flij = m_fli;
      if (m_flj==m_fli) m_ftype = 4;
      else if (!m_fli.IsSusy()) m_ftype = 1;
      else if (m_fli.IsGluino()) m_ftype = 5;
      else if (m_fli.IsSquark()) m_ftype = 7;
      else THROW(fatal_error,"SUSY particle in dipole term, but not squark or gluino");
    }
    else if (m_flj==m_fli.Bar()) {
      if (m_flj.Mass()>m_maxgsmass) {
	m_ftype = 0;
	break;
      }
      if (!m_fli.IsSusy()) m_ftype = 3;
      else {
        m_ftype = 0;
        break;
      }
      m_flij = Flavour(kf_gluon);
    }
    break;
  case dpt::i_f:
  case dpt::i_fm:
  case dpt::i_i:
    if (m_fli==Flavour(kf_gluon)) {
      m_flij = m_flj.Bar();
      if (m_flj==m_fli) m_ftype = 4;
      else m_ftype = 2;
    }
    else if (m_flj==Flavour(kf_gluon)) {
      m_flij = m_fli;
      if (m_flj==m_fli) m_ftype = 4;
      else m_ftype = 1;
    }
    else if (m_flj==m_fli) {
      m_ftype = 3;
      m_flij = Flavour(kf_gluon);
    }
    break;
  default:
    m_ftype = 0;
  }

  if (m_ftype==0) m_valid=false;
  return m_valid;
}

vector<Subprocess_Info>::iterator
Single_DipoleTerm::FindInInfo(Subprocess_Info& fi, int idx) const
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

void Single_DipoleTerm::SetLOMomenta(const Vec4D* moms,const ATOOLS::Poincare &cms)
{
  size_t cnt=0;
  size_t em=p_LO_process->GetEmit();
  size_t sp=p_LO_process->GetSpect();
  if (em==sp) {
    THROW(fatal_error,"Incorrect emitter and spectator assignments.");
  }
  for (size_t i=0;i<m_nin+m_nout;i++) {
    for (;cnt==em||cnt==sp;) cnt++;
    if (i!=m_pi&&i!=m_pj&&i!=m_pk) {
      p_LO_labmom[cnt] = p_LO_mom[cnt] = (*(p_dipole->GetMomenta()))[i];
      cnt++;
    }
  }
  
  p_LO_labmom[em] = p_LO_mom[em] = p_dipole->Getptij();
  p_LO_labmom[sp] = p_LO_mom[sp] = p_dipole->Getptk();

  if (p_LO_mom[0][3]<0.) {
    for (size_t i=0;i<m_nin+m_nout-1;i++) 
      p_LO_mom[i]=Vec4D(p_LO_mom[i][0],-p_LO_mom[i][1],-p_LO_mom[i][2],-p_LO_mom[i][3]);
  }
  Poincare bst(p_LO_mom[0]+p_LO_mom[1]);
  for (size_t i=0;i<m_nin+m_nout-1;i++) bst.Boost(p_LO_mom[i]);
  size_t ndip=(p_dipole->GetDiPolarizations())->size();
  for (size_t i=0;i<ndip;i++) bst.Boost((*(p_dipole->GetDiPolarizations()))[i]);

  p_LO_mom[0]=Vec4D(p_LO_mom[0][0],0.,0.,p_LO_mom[0][0]);
  p_LO_mom[1]=Vec4D(p_LO_mom[0][0],0.,0.,-p_LO_mom[0][0]);

  for (size_t i=0;i<m_nin+m_nout-1;++i) cms.BoostBack(p_LO_labmom[i]);
}

bool Single_DipoleTerm::CompareLOmom(const ATOOLS::Vec4D* p)
{
  for (size_t i=0;i<m_nin+m_nout-1;i++) if (!(p[i]==p_LO_mom[i])) return 0;
  return 1;
}

void Single_DipoleTerm::PrintLOmom()
{
  if (this!=p_partner) { p_partner->PrintLOmom();return;}
  for (size_t i=0;i<m_nin+m_nout-1;i++) cout<<i<<": "<<p_LO_mom[i]<<endl;
}

/*------------------------------------------------------------------------------

  Initializing libraries, amplitudes, etc.

  ------------------------------------------------------------------------------*/



int Single_DipoleTerm::InitAmplitude(Model_Base *model,Topology* top,
				    vector<Process_Base *> & links,
				    vector<Process_Base *> & errs)
{
  switch (m_dipoletype) {
  case dpt::f_f: 
    p_dipole = new FF_DipoleSplitting(m_ftype,m_nin+m_nout-1,m_pi,m_pj,m_pk);
    break;
  case dpt::f_i: 
    p_dipole = new FI_DipoleSplitting(m_ftype,m_nin+m_nout-1,m_pi,m_pj,m_pk);
    break;
  case dpt::i_f: 
    p_dipole = new IF_DipoleSplitting(m_ftype,m_nin+m_nout-1,m_pi,m_pj,m_pk);
    break;
  case dpt::i_i: 
    p_dipole = new II_DipoleSplitting(m_ftype,m_nin+m_nout-1,m_pi,m_pj,m_pk);
    break;
  case dpt::f_fm: 
    p_dipole = new FF_MassiveDipoleSplitting(m_ftype,m_nin+m_nout-1,m_pi,m_pj,m_pk,
					     m_fli.Mass(),m_flj.Mass(),m_flk.Mass(),m_flij.Mass());
    break;
  case dpt::f_im: 
    p_dipole = new FI_MassiveDipoleSplitting(m_ftype,m_nin+m_nout-1,m_pi,m_pj,m_pk,
					     m_fli.Mass(),m_flj.Mass(),m_flij.Mass());
    break;
  case dpt::i_fm: 
    p_dipole = new IF_MassiveDipoleSplitting(m_ftype,m_nin+m_nout-1,m_pi,m_pj,m_pk);
    break;
  default:
    p_dipole=NULL;
  }
  if (!p_dipole) {
    msg_Error()<<"ERROR in Single_DipoleTerm::InitAmplitude : Dipol type not implemented: "<<m_dipoletype
	       <<" ("<<m_pi<<","<<m_pj<<","<<m_pk<<")"<<std::endl;   
    abort();
  }
  p_dipole->SetMomenta(p_testmoms);
  p_dipole->CalcDiPolarizations();
  Poincare cms;
  SetLOMomenta(p_testmoms,cms);

  int status=p_LO_process->InitAmplitude(model,top,links,errs,
					 p_dipole->GetDiPolarizations(),p_dipole->GetFactors());
  if (status<=0) { 
    m_valid=0;
    return status;
  }
  SetOrderQCD(p_LO_process->OrderQCD()+1);
  SetOrderEW(p_LO_process->OrderEW());

  p_dipole->SetCoupling(((Single_LOProcess*)p_LO_process->Partner())->CouplingMap());
  p_dipole->SetAlpha(m_dalpha);

  map<string,Complex> dummy;
  for (size_t j=0;j<links.size();j++) if (Type()==links[j]->Type()) {
    Single_DipoleTerm* check = (Single_DipoleTerm*)links[j];
    bool canmap=0;
    m_sfactor=1.;
    if (GetSplitConfID()==check->GetSplitConfID() && CompareLOmom(check->GetLOmom())) {
      if (p_LO_process->Name()==check->GetLOProcess()->Name()) {
	canmap=1;
	m_sfactor=check->GetSFactor();
      }
      else {
	if (p_LO_process->LibName()==check->GetLOProcess()->LibName()) {
	  if (GetAmplitudeHandler()->CompareAmplitudes(check->GetAmplitudeHandler(),m_sfactor,dummy)) canmap=1;
	  m_sfactor=sqr(m_sfactor);
	}
      }
    }
    if (canmap) {
      msg_Tracking()<<"Can map Dipole Term: "<<Name()<<" -> "<<check->Name()<<" Factor: "<<m_sfactor<<endl;
      p_partner = check;
      Minimize();  //can be switched on if actually mapped!
      m_subevt.m_nco = p_partner->m_subevt.m_nco;
      return 1;
    }
  }
  links.push_back(this);
  Complex C(p_LO_process->GetAmplitudeHandler()->CommonColorFactor());
  if (C!=Complex(0.,0.)) {
    if (C.real()<0) m_subevt.m_nco=1;
    else m_subevt.m_nco=2;
  }

  return 1;
}




/*------------------------------------------------------------------------------
  
  Phase space initialization
  
  ------------------------------------------------------------------------------*/


bool Single_DipoleTerm::SetUpIntegrator() 
{  
  bool res=p_LO_process->SetUpIntegrator();
  if (res) return res;
  return true;
}


/*------------------------------------------------------------------------------
  
  Process management
  
  ------------------------------------------------------------------------------*/
void Single_DipoleTerm::SetLookUp(const bool lookup)
{
  m_lookup=lookup; 
  if (p_LO_process) p_LO_process->SetLookUp(lookup);
}

void Single_DipoleTerm::Minimize()
{
  if (p_partner==this) return;
  if (m_pinfo.m_amegicmhv==0&&p_LO_process) 
    {delete p_LO_process; p_LO_process=0;}
  if (p_LO_mom)     {delete[] p_LO_mom; p_LO_mom=0;}
  if (p_dipole)     {delete p_dipole; p_dipole=0;}
  m_subevt.p_mom = p_partner->GetSubevt()->p_mom;
}


double Single_DipoleTerm::Partonic(const Vec4D_Vector &_moms,const int mode) { return 0.; }

double Single_DipoleTerm::operator()(const ATOOLS::Vec4D * mom,const ATOOLS::Poincare &cms,const int mode)
{
  if (p_partner!=this) {
    if (m_lookup) m_lastxs = p_partner->LastXS()*m_sfactor;
    else m_lastxs = p_partner->operator()(mom,cms,mode)*m_sfactor;
    m_subevt.m_result = m_subevt.m_last[0] = m_subevt.m_last[1] = 0.;
    m_subevt.m_me = m_subevt.m_mewgt = -m_lastxs;
    m_subevt.m_muf2 = p_partner->GetSubevt()->m_muf2;
    m_subevt.m_mur2 = p_partner->GetSubevt()->m_mur2;
    return m_lastxs;
  }

  ResetLastXS();
  p_LO_process->ResetLastXS();
  p_dipole->SetMomenta(mom);
  p_dipole->CalcDiPolarizations();
  SetLOMomenta(mom,cms);

  bool trg(false);
  trg= p_LO_process->Trigger(p_LO_labmom) || !p_LO_process->Selector()->On();
  p_int->SetMomenta(p_LO_labmom);
  p_LO_process->Integrator()->SetMomenta(p_LO_labmom);

  double M2 =trg ? p_LO_process->operator()
    (p_LO_labmom,p_LO_mom,p_dipole->GetFactors(),
     p_dipole->GetDiPolarizations(),mode) : 0.0;
  double df = p_dipole->GetF();
  m_subevt.m_me = m_subevt.m_mewgt = m_subevt.m_result =
    m_subevt.m_last[0] = m_subevt.m_last[1] = 0.;

  if (!(df>0.)&& !(df<0.)) return m_lastxs=df;

  if (!trg) return m_lastxs=m_subevt.m_me=0.;

  m_lastxs = M2 * df * KFactor();
  m_subevt.m_me = m_subevt.m_mewgt = -m_lastxs;
  m_subevt.m_muf2 = p_scale->Scale(stp::fac);
  m_subevt.m_mur2 = p_scale->Scale(stp::ren);
  return m_lastxs;
}

void Single_DipoleTerm::SetSelector(const PHASIC::Selector_Key &key)
{
  p_LO_process->SetSelector(key);
}

int Single_DipoleTerm::NumberOfDiagrams() { 
  if (p_partner==this) return p_LO_process->NumberOfDiagrams(); 
  return p_partner->NumberOfDiagrams();
}

Point * Single_DipoleTerm::Diagram(int i) { 
  if (p_partner==this) return p_LO_process->Diagram(i); 
  return p_partner->Diagram(i);
} 

void Single_DipoleTerm::AddChannels(std::list<std::string>*psln)
{
  p_LO_process->AddChannels(psln);
}

void Single_DipoleTerm::FillAmplitudes(METOOLS::Amplitude_Tensor* atensor,double sfactor)
{
  if (p_partner==this) p_LO_process->FillAmplitudes(atensor,sfactor);
  else p_partner->FillAmplitudes(atensor,sfactor*sqrt(m_sfactor));
}

std::string Single_DipoleTerm::GetSplitConfID()
{
  return ToString(m_pi)+"_"+ToString(m_pj)+"_"+ToString(m_pk)+"_"+ToString(m_dipoletype)+"_"+ToString(m_ftype);
}

void Single_DipoleTerm::PrintProcessSummary(int it)
{
  for(int i=0;i<it;i++) cout<<"  ";
  cout<<m_pi<<"-"<<m_pj<<"-"<<m_pk<<" ("<<p_LO_process->Name()<<")";
  if (p_partner!=this) {
    cout<<"; partner (*"<<m_sfactor<<"): ";
    p_partner->PrintProcessSummary(0);
    return;
  }
  cout<<endl;
} 

void Single_DipoleTerm::SetScale(const Scale_Setter_Arguments &args)
{
  if (!p_LO_process->IsMapped()) p_LO_process->SetScale(args);
  p_scale=p_LO_process->Partner()->ScaleSetter();
}
