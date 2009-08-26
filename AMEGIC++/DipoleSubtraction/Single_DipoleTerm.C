#include "AMEGIC++/DipoleSubtraction/Single_DipoleTerm.H"
#include "AMEGIC++/DipoleSubtraction/Single_LOProcess.H"
#include "AMEGIC++/DipoleSubtraction/Single_LOProcess_MHV.H"

#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
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
  p_partner(this), p_LO_process(0), p_LO_mom(0), m_ftype(0), p_dipole(0), p_realint(pint)
{
  m_pinfo=pinfo;
  m_nin=m_pinfo.m_ii.NExternal();
  m_nout=m_pinfo.m_fi.NExternal();
  m_flavs.resize(m_nin+m_nout);
  if (m_pinfo.m_ii.m_ps.size()>0 && m_pinfo.m_fi.m_ps.size()>0) {
    SortFlavours(m_pinfo);
    m_new=0;
    m_nqcd=0;
    std::vector<Flavour> fl;
    m_pinfo.m_ii.GetExternal(fl);
    m_pinfo.m_fi.GetExternal(fl);
    if (fl.size()!=m_nin+m_nout) THROW(fatal_error,"Internal error");
    for (size_t i(0);i<fl.size();++i) {
      m_flavs[i]=fl[i];
      if (m_flavs[i].Strong()) ++m_nqcd;
      else ++m_new;
    }
    m_name=GenerateName(m_pinfo.m_ii,m_pinfo.m_fi);
  }

  Init();

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
    vector<Subprocess_Info>::iterator sit=lopi.m_fi.m_ps.begin()+m_pj-2;
    lopi.m_fi.m_ps.erase(sit);
  }
  else {
    vector<Subprocess_Info>::iterator sit=lopi.m_fi.m_ps.begin()+m_pi-2;
    lopi.m_fi.m_ps.erase(sit);
    sit=lopi.m_fi.m_ps.begin()+m_pj-m_nin-1;
    lopi.m_fi.m_ps.erase(sit);
    lopi.m_fi.m_ps.push_back(Subprocess_Info(m_flij,m_pinfo.m_fi.m_ps[m_pi-m_nin].m_id,
					     m_pinfo.m_fi.m_ps[m_pi-m_nin].m_pol));
    lopi.m_fi.m_ps[m_LOpij-m_nin].m_tag=1;
  }
  if (m_LOpk<m_nin) lopi.m_ii.m_ps[m_LOpk].m_tag=2;
  else lopi.m_fi.m_ps[m_LOpk-m_nin].m_tag=2;

  if (lopi.m_amegicmhv>0) {
    vector<ATOOLS::Flavour> flin;
    vector<ATOOLS::Flavour> flout;
    lopi.m_ii.GetExternal(flin);
    lopi.m_fi.GetExternal(flout);
    if (CF.MHVCalculable(flin,flout)) p_LO_process = new Single_LOProcess_MHV(lopi);
    if (lopi.m_amegicmhv==2) { m_valid=0; return; }
  }
  if (!p_LO_process) p_LO_process = new Single_LOProcess(lopi);
  if (!p_LO_process) THROW(fatal_error,"LO process unknown");
  if (m_pinfo.m_nlomode==0)
    p_LO_process->Get<PHASIC::Process_Base>()
      ->Init(p_LO_process->Info(),p_realint->Beam(),p_realint->ISR());

  p_LO_mom = new Vec4D[m_nin+m_nout-1];
  p_LO_labmom.resize(m_nin+m_nout-1); 
  p_LO_process->SetTestMoms(p_LO_mom);

  m_subevt.n      = m_nin+m_nout-1;
  m_subevt.p_fl   = &(p_LO_process->Flavours().front());
  m_subevt.p_mom  = &p_LO_labmom.front();
  m_subevt.m_ID   = ToString(m_pi)+"_"+ToString(m_pj)+"_"+ToString(m_pk);

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

  m_maxgsmass=0.;
  int helpi;
  if (reader.ReadFromFile(helpi,"DIPOLE_NF_GSPLIT")) {
    msg_Tracking()<<"Set number of flavours from gluon splitting="<<helpi<<"."<<std::endl;
    Flavour flav((kf_code)(helpi));
    m_maxgsmass=flav.Mass();
  }
}


Single_DipoleTerm::~Single_DipoleTerm()
{
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
      else m_ftype = 2;
    }
    else if (m_flj==Flavour(kf_gluon)) {
      m_flij = m_fli;
      if (m_flj==m_fli) m_ftype = 4;
      else m_ftype = 1;
    }
    else if (m_flj==m_fli.Bar()) {
      if (m_flj.Mass()>m_maxgsmass) {
	m_ftype = 0;
	break;
      }
      m_ftype = 3;
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

void Single_DipoleTerm::SetLOMomenta(const Vec4D* moms) 
{
  size_t cnt=0;
  size_t em=p_LO_process->GetEmit();
  size_t sp=p_LO_process->GetSpect();
  if (em==sp) {
    em=2; sp=3;
    PRINT_INFO("?????????????");
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

   p_realint->ISR()->BoostInLab(&p_LO_labmom.front(),m_nin+m_nout-1);
//    Vec4D msum(0.,0.,0.,0.);
//   cout<<p_dipole->GetType()<<endl;
//    cout<<"Moms: for "<<m_name<<endl;
   
//        for (size_t i=0;i<m_nin+m_nout;i++) {
// 	 cout<<i<<": "<<moms[i]<<" "<<moms[i].Abs2()<<endl;
//        }
//         cout<<"LOMoms: ("<<p_LO_process->Name()<<")"<<endl;
//          for (size_t i=0;i<m_nin+m_nout-1;i++) {
//  	  if (i<2) msum-=p_LO_mom[i];
//  	  else msum+=p_LO_mom[i];
//            cout<<i<<": "<<p_LO_mom[i]<<" "<<p_LO_mom[i].Abs2()<<endl;
//         }
//  	cout<<"sum: "<<msum<<endl;
 
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



int Single_DipoleTerm::InitAmplitude(Model_Base * model,Topology* top,
				    vector<Process_Base *> & links,
				    vector<Process_Base *> & errs)
{
  switch (m_dipoletype) {
  case dpt::f_f: 
    p_dipole = new FF_DipoleSplitting(model,m_ftype,m_nin+m_nout-1,m_pi,m_pj,m_pk);
    break;
  case dpt::f_i: 
    p_dipole = new FI_DipoleSplitting(model,m_ftype,m_nin+m_nout-1,m_pi,m_pj,m_pk);
    break;
  case dpt::i_f: 
    p_dipole = new IF_DipoleSplitting(model,m_ftype,m_nin+m_nout-1,m_pi,m_pj,m_pk);
    break;
  case dpt::i_i: 
    p_dipole = new II_DipoleSplitting(model,m_ftype,m_nin+m_nout-1,m_pi,m_pj,m_pk);
    break;
  case dpt::f_fm: 
    p_dipole = new FF_MassiveDipoleSplitting(model,m_ftype,m_nin+m_nout-1,m_pi,m_pj,m_pk,
					     m_fli.Mass(),m_flj.Mass(),m_flk.Mass(),m_flij.Mass());
    break;
  case dpt::f_im: 
    p_dipole = new FI_MassiveDipoleSplitting(model,m_ftype,m_nin+m_nout-1,m_pi,m_pj,m_pk,
					     m_fli.Mass(),m_flj.Mass(),m_flij.Mass());
    break;
  case dpt::i_fm: 
    p_dipole = new IF_MassiveDipoleSplitting(model,m_ftype,m_nin+m_nout-1,m_pi,m_pj,m_pk);
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
  SetLOMomenta(p_testmoms);

  int status=p_LO_process->InitAmplitude(model,top,links,errs,
					 p_dipole->GetDiPolarizations(),p_dipole->GetFactors());
  if (status<=0) { 
    m_valid=0;
    return status;
  }
  SetOrderQCD(p_LO_process->OrderQCD()+1);
  SetOrderEW(p_LO_process->OrderEW());

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
      return 1;
    }
  }
  links.push_back(this);

  return 1;
}




/*------------------------------------------------------------------------------
  
  Phase space initialization
  
  ------------------------------------------------------------------------------*/


bool Single_DipoleTerm::SetUpIntegrator() 
{  
  return 1;
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


double Single_DipoleTerm::Differential(const Vec4D_Vector &_moms) { return 0.; }

double Single_DipoleTerm::Differential2() { return 0.; }


double Single_DipoleTerm::GetDPSF(const ATOOLS::Vec4D * mom)
{
  p_dipole->SetMomenta(mom);
  SetLOMomenta(mom);
  bool trg = p_realint->Process()->JetTrigger(p_LO_labmom,p_LO_process->Flavours(),m_nout-1);
  return trg?p_dipole->GetDPSF():-p_dipole->GetDPSF();
}

double Single_DipoleTerm::operator()(const ATOOLS::Vec4D * mom)
{
  if (p_partner!=this) {
    if (m_lookup) m_lastxs = p_partner->LastXS()*m_sfactor;
    else m_lastxs = p_partner->operator()(mom)*m_sfactor;
    m_subevt.m_me = m_subevt.m_result = -m_lastxs;
    m_subevt.m_scale = p_partner->GetSubevt()->m_scale;
    m_subevt.m_alpha = p_partner->GetSubevt()->m_alpha;
    return m_lastxs;
  }

  ResetLastXS();
  p_LO_process->ResetLastXS();
  p_dipole->SetMomenta(mom);
  p_dipole->CalcDiPolarizations();
  SetLOMomenta(mom);

  double df = p_dipole->GetF();
  if (!(df>0.)&& !(df<0.))
    return m_lastxs=df;

  bool trg(false);
  if (m_pinfo.m_nlomode==0) trg=p_LO_process->JetTrigger(p_LO_labmom);
  else trg=p_realint->Process()->JetTrigger(p_LO_labmom,p_LO_process->Flavours(),m_nout-1);
  if (!trg) return m_lastxs=0.;

  m_subevt.m_scale = p_scale->CalculateScale(p_LO_labmom);
  double M2 = p_LO_process->operator()(p_LO_mom,p_dipole->GetFactors(),p_dipole->GetDiPolarizations());

  m_lastxs = M2 * df * KFactor();
  m_subevt.m_me = m_subevt.m_result = -m_lastxs;
  m_subevt.m_alpha = p_dipole->GetDPSF();
  return m_lastxs;
}



int Single_DipoleTerm::NumberOfDiagrams() { 
  if (p_partner==this) return p_LO_process->NumberOfDiagrams(); 
  return p_partner->NumberOfDiagrams();
}

Point * Single_DipoleTerm::Diagram(int i) { 
  if (p_partner==this) return p_LO_process->Diagram(i); 
  return p_partner->Diagram(i);
} 

void Single_DipoleTerm::AddChannels(std::list<std::string>*) { }

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
