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
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_MPI.H"

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

/*----------------------------------------------------------------------------

  Constructors

  ----------------------------------------------------------------------------*/

Single_DipoleTerm::Single_DipoleTerm(const Process_Info &pinfo,
                                     size_t pi,size_t pj,size_t pk,
                                     ATOOLS::sbt::subtype stype,
                                     Process_Integrator* pint) :
  m_iresult(0.), m_valid(true), m_noISclustertolepton(true), m_dalpha(1.),
  m_dkt2max(std::numeric_limits<double>::max()),
  m_maxgsmass(0.), p_partner(this), p_LO_process(NULL), p_LO_mom(NULL),
  m_stype(stype), m_dtype(dpt::none), m_ftype(spt::none),
  m_pi(pi), m_pj(pj), m_pk(pk), m_LOpij(0), m_LOpk(0),
  p_dipole(NULL), m_smth(0.), m_nphotonsplits(0),
  m_pspissplitscheme(1), m_pspfssplitscheme(1),
  p_realint(pint)
{
  DEBUG_FUNC(m_stype<<"[("<<pi<<","<<pj<<");"<<pk<<"]");
  PHASIC::Process_Base::Init(pinfo, pint->Beam(), pint->ISR());
  AMEGIC::Process_Base::Init();

  m_name+= "_RS"+ToString(m_pi)+"_"+ToString(m_pj)+"_"+ToString(m_pk);

  // read in g->QQ option
  Flavour flav((kf_code)(ToType<int>(rpa->gen.Variable("DIPOLE_NF_GSPLIT"))));
  m_maxgsmass=flav.Mass();

  // read in P->FF split scheme
  // 0 - do not cluster P->FF
  // 1 - cluster P->QQ, but not P->LL
  // 2 - cluster all P->FF
  m_pspissplitscheme
    =ToType<size_t>(rpa->gen.Variable("DIPOLE_PFF_IS_SPLIT_SCHEME"));
  m_pspfssplitscheme
    =ToType<size_t>(rpa->gen.Variable("DIPOLE_PFF_FS_SPLIT_SCHEME"));
  m_noISclustertolepton
    =ToType<size_t>(rpa->gen.Variable("DIPOLE_IS_CLUSTER_TO_LEPTONS"));

  if (!DetermineType()) return;

  // determine LO indices
  m_LOpij = m_pi;
  size_t LOpj=m_pj;
  m_LOpk  = m_pk;
  if (m_pi>m_nin && m_flavs[m_pi].IsVector() && m_flavs[m_pj].IsFermion())
    std::swap<size_t>(m_LOpij,LOpj);

  // Construct LO Process
  msg_Debugging()<<m_pinfo<<std::endl;
  Process_Info lopi(m_pinfo);
  // reset couplings order
  if (m_stype&sbt::qcd) {
    lopi.m_maxcpl[0]--;
    lopi.m_mincpl[0]--;
  }
  if (m_stype&sbt::qed) {
    lopi.m_maxcpl[1]--;
    lopi.m_mincpl[1]--;
  }

  // set tags
  for (size_t i(0);i<m_nin;++i) lopi.m_ii.m_ps[i].m_tag=i;
  int tag=m_nin;
  lopi.m_fi.SetTags(tag);
  std::vector<int> tags;
  lopi.m_fi.GetTags(tags);
  if (tag!=m_nin+m_nout) {
    THROW(fatal_error, "Internal error");
  }
  if (m_LOpij<m_nin) lopi.m_ii.m_ps[m_LOpij].m_tag=-1;
  else tags[m_LOpij-m_nin]=-1;
  if (m_LOpk<m_nin) lopi.m_ii.m_ps[m_LOpk].m_tag=-2;
  else tags[m_LOpk-m_nin]=-2;
  lopi.m_fi.SetTags(tags);

  // set LO flavours
  if (m_pi<m_nin) {
    lopi.m_ii.m_ps[m_pi].m_fl=m_flij;
  }
  else {
    Flavour_Vector flavs(lopi.m_fi.GetExternal());
    flavs[m_LOpij-m_nin]=m_flij;
    lopi.m_fi.SetExternal(flavs);
  }
  // remove extra parton
  bool found(false);
  for (std::vector<Subprocess_Info>::iterator
         it(lopi.m_fi.m_ps.begin());it!=lopi.m_fi.m_ps.end();++it)
    if (it->m_tag==LOpj) {
      lopi.m_fi.m_ps.erase(it);
      found=true;
      break;
    }
  if (!found) THROW(fatal_error,"Could not identify which particle to remove.");

  lopi.m_fi.m_nlotype=nlo_type::rsub;
  msg_Debugging()<<lopi<<std::endl;

  if (lopi.m_amegicmhv>0) {
    if (lopi.m_amegicmhv==12)
      p_LO_process = new Single_LOProcess_External(lopi, p_int->Beam(),
                                                         p_int->ISR(), m_stype);
    else if (CF.MHVCalculable(lopi))
      p_LO_process = new Single_LOProcess_MHV(lopi, p_int->Beam(),
                                                    p_int->ISR(), m_stype);
    if (lopi.m_amegicmhv==2) { m_valid=false; return; }
  }
  if (!p_LO_process)
    p_LO_process = new Single_LOProcess(lopi, p_int->Beam(),
                                              p_int->ISR(), m_stype);
  if (!p_LO_process) THROW(fatal_error,"LO process unknown");

  if (!(m_valid=p_LO_process->IsValid())) return;

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
  m_subevt.m_oqcd = m_maxcpl[0]/2;
  m_subevt.m_oew  = m_maxcpl[1]/2;
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

  p_LO_process->SetSubEvt(&m_subevt);

  m_dalpha = ToType<double>(rpa->gen.Variable("DIPOLE_ALPHA"));
  m_dkt2max = ToType<double>(rpa->gen.Variable("DIPOLE_KT2MAX"));
  switch (m_dtype) {
  case dpt::f_f:
  case dpt::f_fm:
    m_dalpha = ToType<double>(rpa->gen.Variable("DIPOLE_ALPHA_FF"));
    break;
  case dpt::f_i:
  case dpt::f_im:
    m_dalpha = ToType<double>(rpa->gen.Variable("DIPOLE_ALPHA_FI"));
    break;
  case dpt::i_f:
  case dpt::i_fm:
    m_dalpha = ToType<double>(rpa->gen.Variable("DIPOLE_ALPHA_IF"));
    break;
  case dpt::i_i:
    m_dalpha = ToType<double>(rpa->gen.Variable("DIPOLE_ALPHA_II"));
    break;
  default:
    break;
  }
  msg_Debugging()<<METHOD<<"(): "<<m_dtype<<": dalpha="<<m_dalpha
                         <<" ,  dkt2max="<<m_dkt2max<<std::endl;
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

/*----------------------------------------------------------------------------

  Generic stuff for initialization of Single_DipoleTermes

  ----------------------------------------------------------------------------*/
bool Single_DipoleTerm::DetermineType() {
  DEBUG_FUNC("");
  if (m_pi>=m_pj) { msg_Debugging()<<"i>j\n"; m_valid=false; }
  if (m_pj<m_nin) { msg_Debugging()<<"j in IS\n"; m_valid=false; }
  if (!m_valid) return false;
  bool massive(false);
  bool massiveini(false);
  m_fli=m_flavs[m_pi];
  m_flj=m_flavs[m_pj];
  m_flk=m_flavs[m_pk];
  // only allow massless emittees for now
  if (m_flj.IsMassive()) return m_valid=false;
  if (m_flk.IsMassive()) massive=true;
  if (massive && m_pk<m_nin) massiveini=true;
  if (m_pi>=m_nin) {
    if (m_fli.IsMassive()||m_flj.IsMassive()) massive=true;
    msg_Debugging()<<"massive: "<<massive<<std::endl;
    if (!massive) {
      if (m_pk>=m_nin) m_dtype = dpt::f_f;
      else m_dtype = dpt::f_i;
    }
    else {
      if (m_pk>=m_nin) m_dtype = dpt::f_fm;
      else m_dtype = dpt::f_im;
    }
  }
  else {
    if (m_fli.IsMassive()) massiveini=true;
    if (m_flj.IsMassive()) massiveini=true;
    msg_Debugging()<<"massive: "<<massive<<std::endl;
    msg_Debugging()<<"massiveini: "<<massiveini<<std::endl;
    if (!massive) {
      if (m_pk>=m_nin) m_dtype = dpt::i_f;
      else m_dtype = dpt::i_i;
    }
    else {
      if (m_pk>=m_nin) m_dtype = dpt::i_fm;
    }
  }

  if (massiveini) {
    msg_Error()<<METHOD<<" Massive intial state subtraction not implemented. Abort."<<endl;
    Abort();
  }

  msg_Debugging()<<"dtype: "<<m_dtype<<std::endl;
  if      (m_stype==sbt::qcd) return DetermineQCDType();
  else if (m_stype==sbt::qed) return DetermineEWType();
  return m_valid=false;
}

bool Single_DipoleTerm::DetermineQCDType()
{
  DEBUG_FUNC(m_dtype<<"[("<<m_fli<<","<<m_flj<<");"<<m_flk<<"]");
  switch (m_dtype) {
  case dpt::f_f:
  case dpt::f_fm:
  case dpt::f_i:
  case dpt::f_im:
    if (m_fli==Flavour(kf_gluon)) {
      m_flij = m_flj;
      if (m_flj==m_fli) m_ftype = spt::g2gg;
      else if (!IsSusy(m_fli)) m_ftype = spt::q2gq;
      else if (IsGluino(m_fli)) m_ftype = spt::s2gs;
      else if (IsSquark(m_fli)) m_ftype = spt::G2gG;
      else THROW(fatal_error,"Unknown dipole for "+m_flij.IDName()
                             +" -> "+m_fli.IDName()+" "+m_flj.IDName());
    }
    else if (m_flj==Flavour(kf_gluon)) {
      m_flij = m_fli;
      if (m_flj==m_fli) m_ftype = spt::g2gg;
      else if (!IsSusy(m_fli)) m_ftype = spt::q2qg;
      else if (IsGluino(m_fli)) m_ftype = spt::s2sg;
      else if (IsSquark(m_fli)) m_ftype = spt::G2Gg;
      else THROW(fatal_error,"Unknown dipole for "+m_flij.IDName()
                             +" -> "+m_fli.IDName()+" "+m_flj.IDName());
    }
    else if (m_flj==m_fli.Bar()) {
      if (m_flj.Mass()>m_maxgsmass) {
        m_ftype = spt::none;
	break;
      }
      if (!IsSusy(m_fli)) m_ftype = spt::g2qq;
      else {
        m_ftype = spt::none;
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
      if (m_flj==m_fli) m_ftype = spt::g2gg;
      else m_ftype = spt::q2gq;
    }
    else if (m_flj==Flavour(kf_gluon)) {
      m_flij = m_fli;
      if (m_flj==m_fli) m_ftype = spt::g2gg;
      else m_ftype = spt::q2qg;
    }
    else if (m_flj==m_fli) {
      m_ftype = spt::g2qq;
      m_flij = Flavour(kf_gluon);
    }
    break;
  default:
    m_ftype = spt::none;
  }
  msg_Debugging()<<"ftype: "<<m_ftype<<std::endl;

  // consistency check
  if (m_ftype==spt::none) {
    return m_valid=false;
  }
  if ((!Flavour(kf_jet).Includes(m_fli) ||
       !Flavour(kf_jet).Includes(m_flj)) &&
      Flavour(kf_jet).Includes(m_flij)) {
    return m_valid=false;
  }
  if (!Flavour(kf_jet).Includes(m_flj) && m_pi<m_nin) {
    return m_valid=false;
  }

  return m_valid;
}

bool Single_DipoleTerm::DetermineEWType()
{
  DEBUG_FUNC(m_dtype<<"[("<<m_fli<<","<<m_flj<<");"<<m_flk<<"]");
  switch (m_dtype) {
  case dpt::f_f:
  case dpt::f_fm:
  case dpt::f_i:
  case dpt::f_im:
    if (m_fli==Flavour(kf_photon)) {
      m_flij = m_flj;
      if      (m_flj.IsFermion())          m_ftype = spt::q2gq;
      else if (m_flj.Kfcode()==kf_Wplus)   m_ftype = spt::V2gV;
      else THROW(not_implemented,
                 "QED Subtraction not implemented for "+m_flj.IDName());
    }
    else if (m_flj==Flavour(kf_photon)) {
      m_flij = m_fli;
      if      (m_fli.IsFermion())          m_ftype = spt::q2qg;
      else if (m_fli.Kfcode()==kf_Wplus)   m_ftype = spt::V2Vg;
      else THROW(not_implemented,
                 "QED Subtraction not implemented for "+m_fli.IDName());
    }
    else if (m_flj==m_fli.Bar()) {
      if (m_flj.Mass()>m_maxgsmass) {
        m_ftype = spt::none;
        break;
      }
      if (m_pspfssplitscheme==0) {
        m_ftype = spt::none;
        break;
      }
      if (m_pspfssplitscheme==1 && m_fli.IsLepton()) {
        m_ftype = spt::none;
        break;
      }
      m_flij = Flavour(kf_photon);
      if (m_fli.IntCharge()) m_ftype = spt::g2qq;
    }
    break;
  case dpt::i_f:
  case dpt::i_fm:
  case dpt::i_i:
    if (m_fli==Flavour(kf_photon)) {
      m_flij = m_flj.Bar();
      if (m_flj.IntCharge()) m_ftype = spt::q2gq;
    }
    else if (m_flj==Flavour(kf_photon)) {
      m_flij = m_fli;
      m_ftype = spt::q2qg;
    }
    else if (m_flj==m_fli) {
      if (m_pspissplitscheme==0) {
        m_ftype = spt::none;
        break;
      }
      m_ftype = spt::g2qq;
      m_flij = Flavour(kf_photon);
    }
    break;
  default:
    m_ftype = spt::none;
  }
  msg_Debugging()<<"ftype: "<<m_ftype<<std::endl;

  // disallow initial state leptons where there was none
  if (m_noISclustertolepton && m_flij.IsLepton() && !m_fli.IsLepton() &&
      (m_dtype==dpt::i_i || m_dtype==dpt::i_f || m_dtype==dpt::i_fm))
    m_valid=false;

  if (m_ftype==spt::none) m_valid=false;
  return m_valid;
}

void Single_DipoleTerm::SetLOMomenta(const Vec4D* moms,
                                     const ATOOLS::Poincare &cms)
{
  p_dipole->SetMomenta(moms);
  p_dipole->CalcDiPolarizations();
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

/*----------------------------------------------------------------------------

  Initializing libraries, amplitudes, etc.

  ----------------------------------------------------------------------------*/



int Single_DipoleTerm::InitAmplitude(Amegic_Model *model,Topology* top,
				    vector<Process_Base *> & links,
				    vector<Process_Base *> & errs)
{
  p_LO_process->SetSubevtList(p_subevtlist);
  DEBUG_FUNC("");
  switch (m_dtype) {
  case dpt::f_f: 
    p_dipole = new FF_DipoleSplitting(m_stype,m_ftype,m_nin+m_nout-1,
                                      m_pi,m_pj,m_pk);
    break;
  case dpt::f_i: 
    p_dipole = new FI_DipoleSplitting(m_stype,m_ftype,m_nin+m_nout-1,
                                      m_pi,m_pj,m_pk);
    break;
  case dpt::i_f: 
    p_dipole = new IF_DipoleSplitting(m_stype,m_ftype,m_nin+m_nout-1,
                                      m_pi,m_pj,m_pk);
    break;
  case dpt::i_i: 
    p_dipole = new II_DipoleSplitting(m_stype,m_ftype,m_nin+m_nout-1,
                                      m_pi,m_pj,m_pk);
    break;
  case dpt::f_fm: 
    p_dipole = new FF_MassiveDipoleSplitting(m_stype,m_ftype,m_nin+m_nout-1,
                                             m_pi,m_pj,m_pk,
                                             m_fli.Mass(),m_flj.Mass(),
                                             m_flk.Mass(),m_flij.Mass());
    break;
  case dpt::f_im: 
    p_dipole = new FI_MassiveDipoleSplitting(m_stype,m_ftype,m_nin+m_nout-1,
                                             m_pi,m_pj,m_pk,
                                             m_fli.Mass(),m_flj.Mass(),
                                             m_flij.Mass());
    break;
  case dpt::i_fm: 
    p_dipole = new IF_MassiveDipoleSplitting(m_stype,m_ftype,m_nin+m_nout-1,
                                             m_pi,m_pj,m_pk);
    break;
  default:
    p_dipole=NULL;
  }
  if (!p_dipole) {
    MyStrStream stream;
    stream << "Dipole type not implemented: " << m_dtype;
    stream << " (" << m_pi << "," << m_pj << "," << m_pk << ")" << std::endl;
    THROW(not_implemented, stream.str());
  }
  p_dipole->SetSubevt(&m_subevt);
  msg_Debugging()<<"Initialised dipole "<<*p_dipole<<std::endl;
  Poincare cms;
  SetLOMomenta(p_testmoms,cms);
  p_dipole->CalcDiPolarizations();

  int status=p_LO_process->InitAmplitude(model,top,links,errs,
                                         p_dipole->GetDiPolarizations(),
                                         p_dipole->GetFactors());
  if (status<=0) { 
    m_valid=false;
    return status;
  }
  SetMaxOrders(p_LO_process->MaxOrders());
  SetMinOrders(p_LO_process->MinOrders());
  if (m_stype==sbt::qcd) {
    SetMaxOrder(0,p_LO_process->MaxOrder(0)+1);
    SetMinOrder(0,p_LO_process->MinOrder(0)+1);
  }
  if (m_stype==sbt::qed) {
    SetMaxOrder(1,p_LO_process->MaxOrder(1)+1);
    SetMinOrder(1,p_LO_process->MinOrder(1)+1);
  }
  msg_Debugging()<<"maxorders: "<<MaxOrders()<<std::endl;
  msg_Debugging()<<"minorders: "<<MinOrders()<<std::endl;

  p_dipole->SetCoupling(((Single_LOProcess*)p_LO_process->Partner())
                        ->CouplingMap());
  p_dipole->SetAlpha(m_dalpha);
  p_dipole->SetKt2Max(m_dkt2max);

  return 1;
}




/*----------------------------------------------------------------------------
  
  Phase space initialization
  
  ----------------------------------------------------------------------------*/


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
  SetLOMomenta(&pp.front(),cms);
  return m_mewgtinfo.m_B=operator()(&pp.front(),cms,mode);
}

double Single_DipoleTerm::operator()(const ATOOLS::Vec4D * mom,
                                     const ATOOLS::Poincare &cms,
                                     const int _mode)
{
  DEBUG_FUNC("mode="<<_mode);
  int mode(_mode&~2);
  if (mode==1) return m_lastxs;
  if (p_partner!=this) THROW(not_implemented,"No!!!");

  ResetLastXS();
  p_LO_process->ResetLastXS();

  ((_mode&2)?p_LO_process->Partner():p_LO_process)->SetSubevtList(p_subevtlist);
  p_scale->SetCaller((_mode&2)?p_LO_process->Partner():p_LO_process);

  if (p_LO_process->Selector()->On())
    m_subevt.m_trig=p_dipole->KinCheck()?p_LO_process->Trigger(p_LO_labmom):0;
  else m_subevt.m_trig=true;
  msg_Debugging()<<"Trigger: "<<m_subevt.m_trig
                 <<" ("<<p_dipole->KinCheck()<<")"<<std::endl;
  p_LO_process->Integrator()->SetMomenta(p_LO_labmom);

  int calc=m_subevt.m_trig;
  if (m_smth) {
    double a=m_smth>0.0?p_dipole->KT2():p_dipole->LastAlpha();
    if (a<dabs(m_smth)) calc=m_subevt.p_real->m_trig;
  }
  double M2 = calc ? p_LO_process->operator()
    (p_LO_labmom,p_LO_mom,p_dipole->GetFactors(),
     p_dipole->GetDiPolarizations(),mode) : 0.0;

  if (m_subevt.p_ampl) m_subevt.p_ampl->Delete();
  m_subevt.p_ampl=NULL;

  if (p_scale->FixedScales().empty()) {
    ClusterAmplitude_Vector &ampls
      (ScaleSetter(1)->Amplitudes());
    if (ampls.size()) {
      m_subevt.p_ampl = ampls.front()->CopyAll();
      if (m_subevt.p_ampl->Splitter()==NULL) {
      std::vector<int> rsm(m_nin+m_nout-1);
      for (size_t i(0);i<rsm.size();++i) {
        int cnt=p_LO_process->RSMap()[i];
	if (cnt==-1) cnt=(1<<m_pi)|(1<<m_pj);
	else if (cnt==-2) cnt=1<<m_pk;
	else cnt=1<<cnt;
	rsm[i]=cnt;
      }
      m_subevt.p_ampl->Leg(m_subevt.m_ijt)->SetK(1<<m_subevt.m_kt);
      for (Cluster_Amplitude *campl(m_subevt.p_ampl);campl;campl=campl->Next()) {
	for (size_t i(0);i<campl->Legs().size();++i) {
	  if (p_LO_process->Partner()!=p_LO_process) {
	    Flavour fl(campl->Leg(i)->Flav());
	    fl=ReMap(i<m_nin?fl.Bar():fl,campl->Leg(i)->Id());
	    campl->Leg(i)->SetFlav(i<m_nin?fl.Bar():fl);
	  }
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
	if (campl->IdNew()) {
	  std::vector<int> ids(ID(campl->IdNew()));
	  size_t id(0);
	  for (size_t j(0);j<ids.size();++j) id|=rsm[ids[j]];
	  campl->SetIdNew(id);
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
    m_subevt.m_mu2[stp::fac] = p_scale->Scale(stp::fac);
    m_subevt.m_mu2[stp::ren] = p_scale->Scale(stp::ren);
    m_subevt.m_mu2[stp::res] = p_scale->Scale(stp::res);
    m_subevt.m_kt2=p_dipole->KT2();
    m_subevt.m_alpha = 0.;
    m_subevt.m_trig = false;
    m_subevt.m_K = 1.0;
    return m_lastxs=(m_mcmode&1)?0.0:df;
  }

  if (m_mcmode && p_dipole->MCSign()<0) df=-df;

  if (msg_LevelIsDebugging()) {
    size_t prec(msg->Precision());
    msg->SetPrecision(16);
    msg_Debugging()<<"  |M|^2 = "<<M2*Norm()<<"\n"
                   <<"+  splf = "<<df<<" * "<<p_dipole->SPFac()
                   <<" = "<<df*p_dipole->SPFac()<<"\n"
                   <<"--------------------------------------------------\n"
                   <<"          "<<M2 * df * p_dipole->SPFac() * Norm()<<"\n"
                   <<"kfac ="<<KFactor()<<" ,  norm="<<Norm()<<std::endl;
    msg->SetPrecision(prec);
  }
  m_lastxs = M2 * df * p_dipole->SPFac() * Norm();
  if (m_lastxs) m_lastxs*=m_lastk=KFactor((m_mcmode&2)?4:0);
  m_subevt.m_K = m_lastk;
  m_subevt.m_me = m_subevt.m_mewgt = -m_lastxs;
  m_subevt.m_mu2[stp::fac] = p_scale->Scale(stp::fac);
  m_subevt.m_mu2[stp::ren] = p_scale->Scale(stp::ren);
  m_subevt.m_mu2[stp::res] = p_scale->Scale(stp::res);
  m_subevt.m_alpha = p_dipole->LastAlpha();
  m_subevt.m_kt2=p_dipole->KT2();
  if (!m_subevt.m_trig) m_lastxs=0.0;
  return m_lastxs;
}

void Single_DipoleTerm::SetChargeFactors()
{
  if (m_stype==sbt::qed) {
    if (msg_LevelIsDebugging()) {
      msg_Out()<<"Set charge factors for "
               <<(m_LOpij<m_nin?m_flij.Bar():m_flij)
               <<" -> "<<(m_pi<m_nin?m_fli.Bar():m_fli)<<" "<<m_flj<<" :  ";
      if (m_flij.IsPhoton())
        msg_Out()<<(m_LOpij<m_nin?-1.:1.)
                   *(m_fli.Strong()?abs(m_fli.StrongCharge()):1.)
                   *sqr(m_fli.Charge())/m_nphotonsplits;
      else
        msg_Out()<<(m_LOpij<m_nin?-1.:1.)*(m_LOpk<m_nin?-1.:1.)
                                  *m_flij.Charge()*m_flk.Charge();
      msg_Out()<<std::endl;
    }
    if (m_flij.IsPhoton())
      p_dipole->SetChargeFactor((m_LOpij<m_nin?-1.:1.)
                                *(m_fli.Strong()?abs(m_fli.StrongCharge()):1.)
                                *sqr(m_fli.Charge())/m_nphotonsplits);
    else
      p_dipole->SetChargeFactor((m_LOpij<m_nin?-1.:1.)*(m_LOpk<m_nin?-1.:1.)
                                *m_flij.Charge()*m_flk.Charge());
  }
  else p_dipole->SetChargeFactor(1.);
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

void Single_DipoleTerm::SetNLOMC(PDF::NLOMC_Base *const mc)
{
  p_nlomc=mc;
  if (p_dipole) p_dipole->SetNLOMC(mc);
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
  size_t cmcmode(p_LO_process->Partner()->SetMCMode(mcmode));
  m_mcmode=mcmode;
  return cmcmode;
}

size_t Single_DipoleTerm::SetClusterMode(const size_t cmode)
{
  size_t ccmode(p_LO_process->SetClusterMode(cmode));
  m_cmode=cmode;
  return ccmode;
}

void Single_DipoleTerm::SetVariationWeights(Variation_Weights *const vw)
{
  Process_Base::SetVariationWeights(vw);
  if (!p_LO_process->IsMapped()) p_LO_process->SetVariationWeights(vw);
  p_LO_process->Partner()->SetVariationWeights(vw);
}

ATOOLS::Flavour Single_DipoleTerm::ReMap(const ATOOLS::Flavour &fl,const size_t &id) const
{
  return p_LO_process->ReMap(fl,id);
}

void Single_DipoleTerm::FillProcessMap(NLOTypeStringProcessMap_Map *apmap)
{
  p_apmap=apmap;
  p_LO_process->SetProcMap(p_apmap);
  if (p_apmap->find(nlo_type::rsub)==p_apmap->end())
    (*p_apmap)[nlo_type::rsub] = new StringProcess_Map();
  std::string fname(m_name);
  size_t len(m_pinfo.m_addname.length());
  if (len) fname=fname.erase(fname.rfind(m_pinfo.m_addname),len);
  (*(*p_apmap)[nlo_type::rsub])[fname]=this;
}
