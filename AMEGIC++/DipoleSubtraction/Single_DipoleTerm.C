#include "AMEGIC++/DipoleSubtraction/Single_DipoleTerm.H"
#include "AMEGIC++/DipoleSubtraction/Single_LOProcess.H"
#include "AMEGIC++/DipoleSubtraction/Single_LOProcess_MHV.H"
#include "AMEGIC++/DipoleSubtraction/Single_LOProcess_External.H"

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
  DEBUG_FUNC("("<<pi<<","<<pj<<") <-> "<<pk);
  PHASIC::Process_Base::Init(pinfo, pint->Beam(), pint->ISR());
  AMEGIC::Process_Base::Init();

  m_pi = pi;
  m_pj = pj;
  m_pk = pk;

  m_name+= "_RS"+ToString(m_pi)+"_"+ToString(m_pj)+"_"+ToString(m_pk);

  Data_Reader reader(" ",";","!","=");
  reader.AddComment("#");
  reader.SetInputPath(rpa->GetPath());
  reader.SetInputFile(rpa->gen.Variable("ME_DATA_FILE"));
  int helpi=reader.GetValue<int>("DIPOLE_NF_GSPLIT",Flavour(kf_quark).Size()/2);
  Flavour flav((kf_code)(helpi));
  m_maxgsmass=flav.Mass();

  bool val=DetermineType();
  if (!val) return;

  m_LOpij = m_nin+m_nout-2;
  if (m_pi<m_nin) m_LOpij = m_pi;
  m_LOpk  = m_pk;
  if (m_pi<m_pk&&m_pi>=m_nin) m_LOpk--;
  if (m_pj<m_pk) m_LOpk--;

  //Construct LO Process
  Process_Info lopi=m_pinfo;
  for (size_t i(0);i<m_nin;++i) lopi.m_ii.m_ps[i].m_tag=i;
  int tag=m_nin;
  lopi.m_fi.SetTags(tag);
  if (tag!=m_nin+m_nout) {
    THROW(fatal_error, "Internal error");
  }

  if (m_pi<m_nin) {
    lopi.m_ii.m_ps[m_pi]=Subprocess_Info(m_flij,m_pinfo.m_ii.m_ps[m_pi].m_id,
					 m_pinfo.m_ii.m_ps[m_pi].m_pol);
    lopi.m_ii.m_ps[m_pi].m_tag=-1;
    lopi.m_fi.m_ps.erase(FindInInfo(lopi.m_fi, m_pj-m_nin));
  }
  else {
    lopi.m_fi.m_ps.erase(FindInInfo(lopi.m_fi, m_pi-m_nin));
    lopi.m_fi.m_ps.erase(FindInInfo(lopi.m_fi, m_pj-m_nin-1));
    vector<Subprocess_Info>::iterator part=FindInInfo(m_pinfo.m_fi, m_pi-m_nin);
    lopi.m_fi.m_ps.push_back(Subprocess_Info(m_flij,part->m_id, part->m_pol));
    lopi.m_fi.m_ps.back().m_tag=-1;
  }
  if (m_LOpk<m_nin) lopi.m_ii.m_ps[m_LOpk].m_tag=-2;
  else FindInInfo(lopi.m_fi, m_LOpk-m_nin)->m_tag=-2;
  DEBUG_VAR(lopi);
  Flavour f0(lopi.m_ii.m_ps[0].m_fl);
  SortFlavours(lopi);

  if (lopi.m_amegicmhv>0) {
    if (lopi.m_amegicmhv==10)
      p_LO_process = new Single_LOProcess_External(lopi, p_int->Beam(), p_int->ISR());
    else if (CF.MHVCalculable(lopi))
      p_LO_process = new Single_LOProcess_MHV(lopi, p_int->Beam(), p_int->ISR());
    if (lopi.m_amegicmhv==2) { m_valid=0; return; }
  }
  if (!p_LO_process)
    p_LO_process = new Single_LOProcess(lopi, p_int->Beam(), p_int->ISR());
  if (!p_LO_process) THROW(fatal_error,"LO process unknown");

  p_LO_mom = new Vec4D[m_nin+m_nout-1];
  p_LO_labmom.resize(m_nin+m_nout-1); 
  p_LO_process->SetTestMoms(p_LO_mom);

  m_lofl=p_LO_process->Flavours();

  m_subevt.m_n    = m_nin+m_nout-1;
  m_subevt.p_fl   = &m_lofl.front();
  m_subevt.p_dec  = &m_decins;
  m_subevt.p_mom  = &p_LO_labmom.front();
  m_subevt.m_i    = m_pi;
  m_subevt.m_j    = m_pj;
  m_subevt.m_k    = m_pk;
  m_subevt.p_proc = this;

  m_sids.resize(m_nin+m_nout-1);
  size_t em=p_LO_process->GetEmit();
  size_t sp=p_LO_process->GetSpect();
  for (size_t i=0;i<m_nin+m_nout;++i) {
    int cnt=p_LO_process->SRMap()[i];
    if (cnt>=0) m_sids[cnt] = 1<<i;
  }
  m_sids[em]=(1<<m_pi)|(1<<m_pj);
  m_sids[sp]=1<<m_pk;
  m_subevt.m_ijt=em;
  m_subevt.m_kt=sp;
  m_subevt.p_id=&m_sids.front();
  Process_Info cpi(p_LO_process->Info());
  m_subevt.m_pname=GenerateName(cpi.m_ii,cpi.m_fi);
  m_subevt.m_pname=m_subevt.m_pname.substr(0,m_subevt.m_pname.rfind("__"));

  m_dalpha = 1.;
  double helpd;
  m_dkt2max = std::numeric_limits<double>::max();
  if (reader.ReadFromFile(helpd,"DIPOLE_ALPHA")) {
    m_dalpha = helpd;
    msg_Tracking()<<"Set dipole cut alpha="<<m_dalpha<<"."<<std::endl;
  }
  if (reader.ReadFromFile(helpd,"DIPOLE_KT2MAX")) {
    m_dkt2max = helpd;
    msg_Tracking()<<"Set dipole cut kt2max="<<m_dkt2max<<"."<<std::endl;
  }
  switch (m_dipoletype) {
  case dpt::f_f:
  case dpt::f_fm:
    if (reader.ReadFromFile(helpd,"DIPOLE_ALPHA_FF")) {
      m_dalpha = helpd;
      msg_Tracking()<<"Set ff dipole cut alpha="<<m_dalpha<<" . "<<std::endl;
    }
    break;
  case dpt::f_i:
  case dpt::f_im:
    if (reader.ReadFromFile(helpd,"DIPOLE_ALPHA_FI")) {
      m_dalpha = helpd;
      msg_Tracking()<<"Set fi dipole cut alpha="<<m_dalpha<<" . "<<std::endl;
    }
    break;
  case dpt::i_f:
  case dpt::i_fm:
    if (reader.ReadFromFile(helpd,"DIPOLE_ALPHA_IF")) {
      m_dalpha = helpd;
      msg_Tracking()<<"Set if dipole cut alpha="<<m_dalpha<<" . "<<std::endl;
    }
    break;
  case dpt::i_i:
    if (reader.ReadFromFile(helpd,"DIPOLE_ALPHA_II")) {
      m_dalpha = helpd;
      msg_Tracking()<<"Set ii dipole cut alpha="<<m_dalpha<<" . "<<std::endl;
    }
    break;
  default:
    break;
  }

}


Single_DipoleTerm::~Single_DipoleTerm()
{
  p_selector=NULL;
  p_kfactor=NULL;
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

  if ((!Flavour(kf_jet).Includes(m_fli) ||
       !Flavour(kf_jet).Includes(m_flj)) &&
      Flavour(kf_jet).Includes(m_flij)) {
    m_valid=false;
  }
  if (!Flavour(kf_jet).Includes(m_flj) && m_pi<m_nin) {
    m_valid=false;
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
    int cnt=p_LO_process->SRMap()[i];
    if (cnt<0) continue;
    p_LO_labmom[cnt] = p_LO_mom[cnt] = (*(p_dipole->GetMomenta()))[i];
  }
  
  p_LO_labmom[em] = p_LO_mom[em] = p_dipole->Getptij();
  p_LO_labmom[sp] = p_LO_mom[sp] = p_dipole->Getptk();

  Poincare bst(p_LO_mom[0]+p_LO_mom[1]);
  for (size_t i=0;i<m_nin+m_nout-1;i++) bst.Boost(p_LO_mom[i]);
  size_t ndip=(p_dipole->GetDiPolarizations())->size();
  for (size_t i=0;i<ndip;i++) bst.Boost((*(p_dipole->GetDiPolarizations()))[i]);

  if (m_subevt.m_i<m_nin &&
      m_subevt.m_ijt!=m_subevt.m_i) {
    for (size_t i=0;i<m_nin+m_nout-1;++i) cms.Boost(p_LO_labmom[i]);
  }
  else {
    for (size_t i=0;i<m_nin+m_nout-1;++i) cms.BoostBack(p_LO_labmom[i]);
  }

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
  p_dipole->SetKt2Max(m_dkt2max);

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
  if (p_partner!=this) p_partner->SetLookUp(lookup);
}

void Single_DipoleTerm::Minimize()
{
  if (p_partner==this) return;
  if (p_LO_mom)     {delete[] p_LO_mom; p_LO_mom=0;}
  if (p_dipole)     {delete p_dipole; p_dipole=0;}
  m_subevt.p_mom = p_partner->GetSubevt()->p_mom;
}

bool Single_DipoleTerm::Trigger(const ATOOLS::Vec4D_Vector &p)
{
  return true;
}

double Single_DipoleTerm::Partonic(const Vec4D_Vector &_moms,const int mode)
{
  p_int->SetMomenta(_moms);
  Poincare cms;
  Vec4D_Vector pp(_moms);
  if (m_nin==2 && p_int->ISR() && p_int->ISR()->On()) {
    cms=Poincare(pp[0]+pp[1]);
    for (size_t i(0);i<pp.size();++i) cms.Boost(pp[i]);
  }
  return operator()(&pp.front(),cms,mode);
}

double Single_DipoleTerm::operator()(const ATOOLS::Vec4D * mom,const ATOOLS::Poincare &cms,const int _mode)
{
  int mode(_mode&~2);
  if (mode==1) return m_lastxs;
  if (p_partner!=this) THROW(not_implemented,"No!!!");

  ResetLastXS();
  p_LO_process->ResetLastXS();
  p_dipole->SetMomenta(mom);
  p_dipole->CalcDiPolarizations();
  SetLOMomenta(mom,cms);

  p_scale->SetCaller((_mode&2)?p_LO_process->Partner():p_LO_process);

  if (p_LO_process->Selector()->On())
    m_subevt.m_trig=p_dipole->KinCheck()?p_LO_process->Trigger(p_LO_labmom):0;
  else m_subevt.m_trig=true;
  p_LO_process->Integrator()->SetMomenta(p_LO_labmom);

  int calc=m_subevt.m_trig;
  if (m_smth) {
    double a=m_smth>0.0?p_dipole->KT2():p_dipole->LastAlpha();
    if (a<dabs(m_smth)) calc=1;
  }
  double M2 = calc ? p_LO_process->operator()
    (p_LO_labmom,p_LO_mom,p_dipole->GetFactors(),
     p_dipole->GetDiPolarizations(),mode) : 0.0;

  if (m_subevt.p_ampl) m_subevt.p_ampl->Delete();
  m_subevt.p_ampl=NULL;

  if (m_subevt.m_trig && p_scale->FixedScales().empty() && (_mode&2)) {
    ClusterAmplitude_Vector &ampls
      (ScaleSetter(1)->Amplitudes());
    if (ampls.size()) {
      m_subevt.p_ampl = Cluster_Amplitude::New();
      m_subevt.p_ampl->SetNext(ampls.front()->CopyAll());
      m_subevt.p_ampl->SetMS(m_subevt.p_ampl->Next()->MS());
      m_subevt.p_ampl->Next()->SetKin(1);
      m_subevt.p_ampl->SetNIn(m_nin);
      for (size_t i(0);i<m_nin;++i)
	m_subevt.p_ampl->CreateLeg(-p_int->Momenta()[i],Flavours()[i].Bar(),ColorID(),1<<i);
      for (size_t i(m_nin);i<m_nin+m_nout;++i)
	m_subevt.p_ampl->CreateLeg(p_int->Momenta()[i],Flavours()[i],ColorID(),1<<i);
      m_subevt.p_ampl->SetMuR2(m_subevt.p_ampl->Next()->MuR2());
      m_subevt.p_ampl->SetMuF2(m_subevt.p_ampl->Next()->MuF2());
      m_subevt.p_ampl->SetQ2(m_subevt.p_ampl->Next()->Q2());
      m_subevt.p_ampl->SetKT2(p_dipole->KT2());
      m_subevt.p_ampl->SetMu2(p_dipole->KT2());
      m_subevt.p_ampl->SetOrderEW(m_subevt.p_ampl->Next()->OrderEW());
      m_subevt.p_ampl->SetOrderQCD(m_subevt.p_ampl->Next()->OrderQCD()+1);
      std::vector<int> rsm(m_nin+m_nout-1);
      for (size_t i(0);i<rsm.size();++i) {
        int cnt=p_LO_process->RSMap()[i];
	if (cnt==-1) cnt=(1<<m_pi)|(1<<m_pj);
	else if (cnt==-2) cnt=1<<m_pk;
	else {
	  m_subevt.p_ampl->Leg(cnt)->SetCol
	    (m_subevt.p_ampl->Next()->Leg(i)->Col());
	  cnt=1<<cnt;
	}
	rsm[i]=cnt;
      }
      m_subevt.p_ampl->Leg(m_subevt.m_k)->SetCol
	     (m_subevt.p_ampl->Next()->Leg(m_subevt.m_kt)->Col());
      m_subevt.p_ampl->Next()->Leg(m_subevt.m_ijt)->SetK(1<<m_subevt.m_kt);
      Cluster_Amplitude::SetColours
	(m_subevt.p_ampl->Next()->Leg(m_subevt.m_ijt),
	 m_subevt.p_ampl->Leg(m_pi),
	 m_subevt.p_ampl->Leg(m_pj));
      for (Cluster_Amplitude *campl(m_subevt.p_ampl->Next());campl;campl=campl->Next()) {
	for (size_t i(0);i<campl->Legs().size();++i) {
	  campl->Leg(i)->SetFlav(ReMap(campl->Leg(i)->Flav(),campl->Leg(i)->Id()));
	  std::vector<int> ids(ID(campl->Leg(i)->Id()));
	  size_t id(0);
	  for (size_t j(0);j<ids.size();++j) id|=rsm[ids[j]];
          campl->Leg(i)->SetId(id);
	       if (campl->Leg(i)->K()) {
		 std::vector<int> ids(ID(campl->Leg(i)->K()));
		 size_t id(0);
		 for (size_t j(0);j<ids.size();++j) id|=rsm[ids[j]];
		      campl->Leg(i)->SetK(id);
	       }
	}
      }     
    }
  }

  p_dipole->SetMCMode(m_mcmode);
  if (m_subevt.m_trig && m_mcmode) {
    p_dipole->SetKt2Max(p_scale->Scale(stp::res));
    if (p_scale->Scales().size()>(stp::size+stp::res)) {
      p_dipole->SetKt2Max(p_scale->Scale(stp::id(stp::size+stp::res)));
    }
  }

  double df = p_dipole->KinCheck()?p_dipole->GetF():nan;
  if (!(df>0.)&& !(df<0.)) {
    m_subevt.m_me = m_subevt.m_mewgt = 0.;
    m_subevt.m_trig = false;
    return m_lastxs=(m_mcmode&1)?0.0:df;
  }

  if (m_mcmode && p_dipole->MCSign()<0) df=-df;

  m_lastxs = M2 * df * p_dipole->SPFac() * KFactor() * Norm();
  m_subevt.m_me = m_subevt.m_mewgt = -m_lastxs;
  m_subevt.m_mu2[stp::fac] = p_scale->Scale(stp::fac);
  m_subevt.m_mu2[stp::ren] = p_scale->Scale(stp::ren);
  m_subevt.m_mu2[stp::res] = p_scale->Scale(stp::res);
  m_subevt.m_kt2=p_dipole->KT2();
  if (!m_subevt.m_trig) m_lastxs=0.0;
  DEBUG_VAR(m_lastxs);
  return m_lastxs;
}

void Single_DipoleTerm::SetSelector(const PHASIC::Selector_Key &key)
{
  if (p_LO_process==NULL) return;
  p_LO_process->SetSelector(key);
  p_selector=p_LO_process->Selector();
}

void Single_DipoleTerm::SetShower(PDF::Shower_Base *const ps)
{
  p_shower=ps;
  if (p_LO_process) p_LO_process->SetShower(ps);
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
  if (p_LO_process==NULL) return;
  p_LO_process->AddChannels(psln);
}

std::string Single_DipoleTerm::GetSplitConfID()
{
  return ToString(m_pi)+"_"+ToString(m_pj)+"_"+ToString(m_pk)+"_"+ToString(m_dipoletype)+"_"+ToString(m_ftype);
}

void Single_DipoleTerm::PrintProcessSummary(int it)
{
  for(int i=0;i<it;i++) cout<<"  ";
  cout<<m_pi<<"-"<<m_pj<<"-"<<m_pk<<" ("<<p_partner->p_LO_process->Name()<<")";
  if (p_partner!=this) {
    cout<<"; partner (*"<<m_sfactor<<"): ";
    p_partner->PrintProcessSummary(0);
    return;
  }
  cout<<endl;
} 

void Single_DipoleTerm::SetScale(const Scale_Setter_Arguments &args)
{
  if (p_LO_process==NULL) return;
  if (!p_LO_process->IsMapped()) p_LO_process->SetScale(args);
  p_scale=p_LO_process->Partner()->ScaleSetter();
}

void Single_DipoleTerm::SetKFactor(const KFactor_Setter_Arguments &args)
{
  if (!p_LO_process->IsMapped()) p_LO_process->SetKFactor(args);
  p_kfactor=p_LO_process->Partner()->KFactorSetter();
}

size_t Single_DipoleTerm::SetMCMode(const size_t mcmode)
{
  size_t cmcmode(p_LO_process->SetMCMode(mcmode));
  m_mcmode=mcmode;
  return cmcmode;
}

size_t Single_DipoleTerm::SetClusterMode(const size_t cmode)
{
  size_t ccmode(p_LO_process->SetClusterMode(cmode));
  m_cmode=cmode;
  return ccmode;
}


ATOOLS::Flavour Single_DipoleTerm::ReMap(const ATOOLS::Flavour &fl,const size_t &id) const
{
  return p_LO_process->ReMap(fl,id);
}
