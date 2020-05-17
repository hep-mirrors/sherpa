#include "AMEGIC++/DipoleSubtraction/Single_Real_Correction.H"
#include "AMEGIC++/Main/Single_Process.H"
#include "AMEGIC++/Main/Single_Process_MHV.H"
#include "AMEGIC++/Main/Single_Process_External.H"
#include "AMEGIC++/DipoleSubtraction/Single_DipoleTerm.H"
#include "AMEGIC++/DipoleSubtraction/Single_OSTerm.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "PDF/Main/ISR_Handler.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Combined_Selector.H"

#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace AMEGIC;
using namespace MODEL;
using namespace PHASIC;
using namespace PDF;
using namespace BEAM;
using namespace ATOOLS;
using namespace std;

/*-------------------------------------------------------------------------------

  Constructors

  ------------------------------------------------------------------------------- */

Single_Real_Correction::Single_Real_Correction() :
  m_newlib(false), m_no_tree(false), m_listdips(false),
  m_iresult(0.), m_smear_threshold(0.), m_smear_power(0.),
  m_libnumb(0), m_ossubon(0),
  m_pspisrecscheme(0), m_pspfsrecscheme(0),
  p_partner(this), p_tree_process(NULL)
{
  DEBUG_FUNC("");
  m_Norm = 1.;
  static bool addcite(false);
  if (!addcite) {
    addcite=true;
    rpa->gen.AddCitation(1,"The automated generation of Catani-Seymour dipole\
 terms in Amegic is published under \\cite{Gleisberg:2007md}.");
  }
  m_smear_threshold=ToType<double>(rpa->gen.Variable("NLO_SMEAR_THRESHOLD"));
  m_smear_power=ToType<double>(rpa->gen.Variable("NLO_SMEAR_POWER"));

  auto dipolesettings = Settings::GetMainSettings()["DIPOLES"];
  m_listdips = dipolesettings["LIST"].Get<size_t>();
  m_pspisrecscheme = dipolesettings["PFF_IS_RECOIL_SCHEME"].Get<size_t>();
  m_pspfsrecscheme = dipolesettings["PFF_FS_RECOIL_SCHEME"].Get<size_t>();
}


Single_Real_Correction::~Single_Real_Correction()
{
  p_scale=NULL;
  p_selector=NULL;
  if (p_tree_process) delete p_tree_process;
  for (size_t i=0;i<m_subtermlist.size();i++) delete m_subtermlist[i];
  for (size_t i=0;i<m_subostermlist.size();i++) delete m_subostermlist[i];
  for (std::map<void*,DM_Info>::const_iterator it(m_dfmap.begin());
       it!=m_dfmap.end();++it) {
    delete it->second.p_fl;
    delete it->second.p_id;
  }
}



/*------------------------------------------------------------------------------

  Initializing libraries, amplitudes, etc.

  ------------------------------------------------------------------------------*/

int Single_Real_Correction::InitAmplitude(Amegic_Model * model,Topology* top,
					vector<Process_Base *> & links,
					vector<Process_Base *> & errs)
{
  DEBUG_FUNC(m_name);
  Init();
  if (!model->p_model->CheckFlavours(m_nin,m_nout,&m_flavs.front())) return 0;

  m_newlib   = false;
  if (m_pinfo.m_amegicmhv>0) {
    if (m_pinfo.m_amegicmhv==10 || m_pinfo.m_amegicmhv==12) {
      p_tree_process = new Single_Process_External();
      m_no_tree=0;
    } else if (CF.MHVCalculable(m_pinfo)) {
      p_tree_process = new Single_Process_MHV();
    }
    if (m_pinfo.m_amegicmhv==2)
      return 0;
  }
  if (!p_tree_process)
    p_tree_process = new AMEGIC::Single_Process();

  p_tree_process->SetSubevtList(&m_subevtlist);

  int status;

  Process_Info rinfo(m_pinfo);
  rinfo.m_fi.m_nlotype=nlo_type::real;
  p_tree_process->PHASIC::Process_Base::Init(rinfo,p_int->Beam(),p_int->ISR());
  p_tree_process->SetTestMoms(p_testmoms);
  p_tree_process->SetPrintGraphs(rinfo.m_gpath);

  // build tree process
  if (m_no_tree) {
    p_tree_process->Init();
    if (dynamic_cast<AMEGIC::Single_Process*>(p_tree_process))
      p_tree_process->Get<AMEGIC::Single_Process>()->PolarizationNorm();
    if (dynamic_cast<AMEGIC::Single_Process_MHV*>(p_tree_process))
      p_tree_process->Get<AMEGIC::Single_Process_MHV>()->PolarizationNorm();
    status=1;
  }
  else {
    status = p_tree_process->InitAmplitude(model,top,links,errs);

    m_maxcpl=p_tree_process->MaxOrders();
    m_mincpl=p_tree_process->MinOrders();
    if (p_tree_process->Partner()->NewLibs()) m_newlib = 1;

    m_iresult=p_tree_process->Result();
    if (status==0) {
      return status;
    }

  if (p_tree_process!=p_tree_process->Partner()) {
    for (size_t j=0;j<links.size();j++) if (Type()==links[j]->Type()) {
      if (p_tree_process->Partner()->Flavours()==links[j]->Flavours()) {
	msg_Tracking()<<"Can map full real process: "<<Name()<<" -> "
		      <<links[j]->Name()<<" Factor: "<<p_tree_process->GetSFactor()<<endl;
	p_mapproc = p_partner = (Single_Real_Correction*)links[j];
	m_sfactor = p_tree_process->GetSFactor();
	// return 1;
      }
    }
  }
  }

  m_real_momenta.resize(m_nin+m_nout);

  m_realevt.m_n    = m_nin+m_nout;
  m_realevt.p_fl   = &p_tree_process->Flavours().front();
  m_realevt.p_dec  = &m_decins;

  m_realevt.p_mom  = &m_real_momenta.front();
  m_realevt.m_i    = m_realevt.m_j = m_realevt.m_k = 0;
  m_realevt.m_oqcd = p_tree_process->MaxOrder(0);
  m_realevt.m_oew  = p_tree_process->MaxOrder(1);

  m_sids.resize(m_nin+m_nout);
  for (size_t i(0);i<m_nin+m_nout;++i) m_sids[i]=1<<i;
  m_realevt.p_id=&m_sids.front();
  m_realevt.m_pname = GenerateName(m_pinfo.m_ii,m_pinfo.m_fi);
  m_realevt.m_pname = m_realevt.m_pname.substr(0,m_realevt.m_pname.rfind("__"));
  m_realevt.p_proc = this;
  m_realevt.p_real = &m_realevt;

  if (p_tree_process==p_tree_process->Partner()) {
  Process_Info cinfo(m_pinfo);
  cinfo.m_fi.m_nlotype=nlo_type::rsub;
  size_t nPFFsplittings(0);
  for (size_t i=0;i<m_flavs.size();i++) {
    for (size_t j=i+1;j<m_flavs.size();j++) {
      for (size_t k=0;k<m_flavs.size();k++) if (k!=i&&k!=j&&i!=j) {
        if (j<m_nin) continue;
        std::vector<sbt::subtype> stypes;
        bool isPFFsplitting(false);
        if (cinfo.m_maxcpl[0]>=1.) {
          if (m_flavs[i].IsQCD() && m_flavs[j].IsQCD() && m_flavs[k].IsQCD()) {
            stypes.push_back(sbt::qcd);
          }
        }
        else msg_Debugging()<<"No QCD subtraction possible."<<std::endl;
        if (cinfo.m_maxcpl[1]>=1.) {
          if (m_flavs[i].IsPhoton() && m_flavs[j].Charge() &&
              m_flavs[k].Charge()) {
            stypes.push_back(sbt::qed);
          }
          else if (m_flavs[i].Charge() && m_flavs[j].IsPhoton() &&
              m_flavs[k].Charge()) {
            stypes.push_back(sbt::qed);
          }
          else if (i<m_nin &&
                   m_flavs[i].Charge() && m_flavs[j].Charge() &&
                   m_flavs[i]==m_flavs[j] && AllowAsSpecInISPFF(k)) {
            stypes.push_back(sbt::qed);
            isPFFsplitting=true;
          }
          else if (i>=m_nin &&
                   m_flavs[i].Charge() && m_flavs[j].Charge() &&
                   m_flavs[i]==m_flavs[j].Bar() && AllowAsSpecInFSPFF(k)) {
            stypes.push_back(sbt::qed);
            isPFFsplitting=true;
          }
        }
        else msg_Debugging()<<"No QED subtraction possible."<<std::endl;
        std::string ststr("");
        for (size_t s(0);s<stypes.size();++s) {
          ststr+=ToString(stypes[s])+" ";
        }
        msg_Debugging()<<"[("<<i<<","<<j<<");"<<k<<"] : "
                       <<(Combinable(1<<i,1<<j)?"":"not ")<<"combinable"
                       <<", types: "<<ststr<<std::endl;
        for (size_t s(0);s<stypes.size();++s) {
	  if (m_pinfo.m_ckkw && stypes[s]==sbt::qed) continue;
          Single_DipoleTerm *pdummy
            = new Single_DipoleTerm(cinfo,i,j,k,stypes[s],p_int);
          msg_Debugging()<<stypes[s]<<"[("<<i<<","<<j<<");"<<k<<"]("
                         <<stypes[s]<<") -> "
                         <<(pdummy->IsValid()?"":"in")<<"valid";
	  if (pdummy->IsValid()) {
	    pdummy->SetRealSubevt(&m_realevt);
            pdummy->SetTestMoms(p_testmoms);
            pdummy->SetNPhotonSplittings(1); // TODO: fix me
            int st=pdummy->InitAmplitude(model,top,links,errs);
            if (!pdummy->IsValid()) msg_Debugging()<<" -> invalid";
            if (pdummy->IsValid()) {
              status=Min(st,status);
              if (pdummy->Partner()->NewLibs()) m_newlib = 1;
              m_subtermlist.push_back(pdummy);
              m_subtermlist.back()->SetNorm(p_tree_process->Norm());
              m_subtermlist.back()->SetSmearThreshold(m_smear_threshold);
              if (isPFFsplitting) nPFFsplittings++;
              m_subevtlist.push_back(pdummy->GetSubevt());
            }
            else delete pdummy;
          }
          else delete pdummy;
          msg_Debugging()<<"\n";
        }
        msg_Debugging()<<"---------------------------------------------\n";
      }
    }
  }
  if (nPFFsplittings)
    for (size_t i(0);i<m_subtermlist.size();++i)
      m_subtermlist[i]->SetNPhotonSplittings(nPFFsplittings);
  for (size_t i(0);i<m_subtermlist.size();++i)
    m_subtermlist[i]->SetChargeFactors();
  if (m_listdips || msg_LevelIsTracking()) {
    msg_Out()<<m_name<<" @ O"<<m_maxcpl<<": "<<m_subtermlist.size()<<" dipoles:\n";
    if (m_listdips || msg_LevelIsDebugging()) {
      for (size_t i(0);i<m_subtermlist.size();++i) {
        Single_DipoleTerm * dt(m_subtermlist[i]);
        msg_Out()<<"  "<<dt->Name()<<":  "
                 <<dt->GetLOProcess()->Flavours()
                 <<" @ O"<<dt->GetLOProcess()->MaxOrders()<<" "
                 <<dt->GetSubtractionType()
                 <<"[("<<dt->Li()<<","<<dt->Lj()<<");"<<dt->Lk()<<"] "
                 <<dt->GetDipoleType()
                 <<"[("<<dt->Flavours()[dt->Li()]<<","<<dt->Flavours()[dt->Lj()]
                 <<");"<<dt->Flavours()[dt->Lk()]<<"]"<<std::endl;
      }
    }
  }
  if (m_ossubon) {
    Process_Info sinfo(m_pinfo);
    sinfo.m_fi.m_nlotype=nlo_type::lo;
    for (size_t i=0;i<m_flavs.size();i++) if (IsSusy(m_flavs[i])){
      for (size_t j=0;j<m_flavs.size();j++) if (i!=j) {
        for (size_t swit=0;swit<5;swit++) {
          Single_OSTerm *pdummy = new Single_OSTerm(
              sinfo,i,j,swit,p_int);
	  if (pdummy->IsValid()) {
            pdummy->SetTestMoms(p_testmoms);
            int st=pdummy->InitAmplitude(model,top,links,errs);
            if (pdummy->IsValid()) {
              status=Min(st,status);
              if (pdummy->NewLibs()) m_newlib = 1;
              m_subostermlist.push_back(pdummy);
              m_subostermlist.back()->SetNorm(p_tree_process->Norm());
	      m_subevtlist.push_back(pdummy->GetSubevt());
            }
            else delete pdummy;
	  }
	  else delete pdummy;
        }
      }
    }
  }
  if (m_no_tree)
    if (m_subtermlist.empty() && m_subostermlist.empty()) return 0;
  }
  m_subevtlist.push_back(&m_realevt);

  if (p_mapproc && !p_partner->NewLibs()) Minimize();
  if (p_partner==this) msg_Info()<<"."<<std::flush;

  if (status>=0) links.push_back(this);
  if (status<0) errs.push_back(this);
  return status;
}



/*------------------------------------------------------------------------------
  
  Phase space initialization
  
  ------------------------------------------------------------------------------*/

bool AMEGIC::Single_Real_Correction::FillIntegrator
(PHASIC::Phase_Space_Handler *const psh)
{
  if (p_partner!=this) return true;
  My_In_File::OpenDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/");
  if (!SetUpIntegrator()) THROW(fatal_error,"No integrator");
  My_In_File::CloseDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/");
  if (m_pinfo.m_nlomode==nlo_mode::powheg) return true;
  return p_tree_process->FillIntegrator(psh);
}


bool AMEGIC::Single_Real_Correction::Combinable
(const size_t &idi,const size_t &idj)
{
  return p_tree_process->Combinable(idi, idj);
}

  
const Flavour_Vector &AMEGIC::Single_Real_Correction::CombinedFlavour
(const size_t &idij)
{
  return p_tree_process->CombinedFlavour(idij);
}

  
bool Single_Real_Correction::SetUpIntegrator() 
{  
  if (m_nin==2) {
    if ( (m_flavs[0].Mass() != p_int->ISR()->Flav(0).Mass()) ||
	 (m_flavs[1].Mass() != p_int->ISR()->Flav(1).Mass()) ) p_int->ISR()->SetPartonMasses(m_flavs);
  }
  return p_tree_process->SetUpIntegrator();
}


/*------------------------------------------------------------------------------
  
  Process management
  
  ------------------------------------------------------------------------------*/
void Single_Real_Correction::SetLookUp(const bool lookup)
{
  m_lookup=lookup; 
  if (p_tree_process) p_tree_process->SetLookUp(false);
  for (size_t i=0;i<m_subtermlist.size();i++) 
    m_subtermlist[i]->SetLookUp(false);
}

void Single_Real_Correction::Minimize()
{
  if (p_partner==this) return;
  p_tree_process->Minimize();
  for (size_t i=0;i<m_subtermlist.size();i++) 
    m_subtermlist[i]->Minimize();
  for (size_t i=0;i<m_subostermlist.size();i++) 
    m_subostermlist[i]->Minimize();
  m_subevtlist.clear();
  for (size_t i=0;i<p_partner->m_subtermlist.size();i++)
    if (p_partner->m_subtermlist[i]->IsValid()) {
      m_subevtlist.push_back(new NLO_subevt(*p_partner->m_subtermlist[i]->GetSubevt()));
      ReMapFlavs(m_subevtlist.back(),1);
    }
  m_subevtlist.push_back(new NLO_subevt(p_partner->m_realevt));
  ReMapFlavs(m_subevtlist.back(),1);
  for (size_t i=0;i<m_subtermlist.size();++i)
    m_subevtlist[i]->p_proc=m_subtermlist[i];
  m_subevtlist.back()->p_proc=this;
  for (size_t i=0;i<m_subevtlist.size();++i)
    m_subevtlist[i]->p_real=m_subevtlist.back();
}

void Single_Real_Correction::ReMapFlavs(NLO_subevt *const sub,const int mode)
{
  if (mode==0) {
    std::map<void*,DM_Info>::const_iterator it(m_dfmap.find((void*)sub->p_fl));
    if (it==m_dfmap.end()) THROW(fatal_error,"Internal error");
    sub->p_fl=&it->second.p_fl->front();
    sub->p_id=&it->second.p_id->front();
    sub->m_pname=it->second.m_pname;
    return;
  }
  Cluster_Amplitude *ampl(Cluster_Amplitude::New());
  ampl->SetNIn(m_nin);
  Flavour_Vector *fls(new Flavour_Vector());
  for (size_t i(0);i<sub->m_n;++i) {
    fls->push_back(p_tree_process->ReMap(sub->p_fl[i],ToString(sub->p_id[i])));
    ampl->CreateLeg(Vec4D(),i<m_nin?fls->back().Bar():fls->back(),ColorID(),sub->p_id[i]);
  }
  ampl->Decays()=*sub->p_dec;
  SortFlavours(ampl);
  std::vector<size_t> *ids(new std::vector<size_t>(sub->m_n));
  for (size_t i(0);i<sub->m_n;++i) {
    (*fls)[i]=i<m_nin?ampl->Leg(i)->Flav().Bar():ampl->Leg(i)->Flav();
    (*ids)[i]=ampl->Leg(i)->Id();
  }
  m_dfmap[(void*)sub->p_fl]=DM_Info(fls,ids,GenerateName(ampl));
  ampl->Delete();
  ReMapFlavs(sub,0);
}

double Single_Real_Correction::Partonic(const ATOOLS::Vec4D_Vector &moms,
                                        int mode)
{
  DEBUG_FUNC("mode="<<mode);
  if (mode==1) return m_lastxs;
  m_lastxs = 0.;

  // fill m_subevtlist
  if (p_partner == this) m_lastdxs = operator()(moms,mode);
  else {
    if (m_lookup) m_lastdxs = p_partner->LastDXS()*m_sfactor;
    else m_lastdxs = p_partner->operator()(moms,mode)*m_sfactor;
    std::vector<NLO_subevt*>* partnerlist=p_partner->GetSubevtList();
    if (partnerlist->size()!=m_subevtlist.size()) THROW(fatal_error,"Internal error");
    for (size_t i=0;i<partnerlist->size();++i) {
      m_subevtlist[i]->CopyXSData((*partnerlist)[i]);
      if ((*partnerlist)[i]->p_ampl) {
	if (i+1<partnerlist->size() && m_subevtlist[i]->p_ampl->IdNew()) {
	  m_subevtlist[i]->p_ampl->Delete();
	  m_subevtlist[i]->p_ampl=(*partnerlist)[i]->p_ampl->First()->CopyAll();
	  if (m_subevtlist[i]->p_ampl->Next()) {
	    m_subevtlist[i]->p_ampl=m_subevtlist[i]->p_ampl->Next();
	    m_subevtlist[i]->p_ampl->DeletePrev();
	  }
	}
      }
    }
    m_subevtlist.Mult(m_sfactor);
  }
  msg_Debugging()<<"RS="<<m_lastdxs<<std::endl;
  return m_lastxs=m_lastdxs;
}

double Single_Real_Correction::operator()(const ATOOLS::Vec4D_Vector &_mom,const int mode)
{
  DEBUG_FUNC(m_name);
  p_tree_process->Integrator()->SetMomenta(_mom);
  for (size_t i=0; i<m_real_momenta.size(); ++i) m_real_momenta[i]=_mom[i];

  Vec4D_Vector mom(_mom);
  Poincare cms;
  if (m_nin==2 && p_int->ISR() && p_int->ISR()->On()) {
    cms=Poincare(mom[0]+mom[1]);
    for (size_t i(0);i<mom.size();++i) cms.Boost(mom[i]);
  }

  bool res(true);
  for (size_t i=0;i<m_subtermlist.size();i++) if (m_subtermlist[i]->IsValid()) {
    m_subtermlist[i]->Integrator()->SetMomenta(_mom);
    double test = (*m_subtermlist[i])(&mom.front(),cms,mode|2);
    if (IsBad(test)) res=false;
    m_subtermlist[i]->GetSubevt()->m_oqcd = m_subtermlist[i]->MaxOrder(0);
    m_subtermlist[i]->GetSubevt()->m_oew = m_subtermlist[i]->MaxOrder(1);
    m_subtermlist[i]->GetSubevt()->p_real=&m_realevt;
  }

  if (m_ossubon){
    for (size_t i=0;i<m_subostermlist.size();i++) if (m_subostermlist[i]->IsValid()) {
      double test = (*m_subostermlist[i])(&mom.front(),cms,mode|2);
      if (IsBad(test)) res=false;
      m_subtermlist[i]->GetSubevt()->m_oqcd = m_subostermlist[i]->MaxOrder(0);
      m_subtermlist[i]->GetSubevt()->m_oew = m_subostermlist[i]->MaxOrder(1);
      m_subostermlist[i]->GetSubevt()->p_real=&m_realevt;
    }
  }
  if (m_smear_threshold!=0.0) SmearCounterEvents(m_subevtlist);

  m_realevt.m_me = m_realevt.m_mewgt = 0.0;
  m_realevt.m_trig = false;
  m_realevt.m_K = 1.0;

  if (m_realevt.p_ampl) m_realevt.p_ampl->Delete();
  m_realevt.p_ampl=NULL;

  bool realtrg(true);

  if (!m_no_tree) {
    realtrg=p_tree_process->Trigger(_mom);
    if (res && realtrg) {
      p_tree_process->ScaleSetter()->CalculateScale(_mom,m_cmode);
      m_realevt.m_mu2[stp::fac]=p_tree_process->ScaleSetter()->Scale(stp::fac);
      m_realevt.m_mu2[stp::ren]=p_tree_process->ScaleSetter()->Scale(stp::ren);
      m_realevt.m_mu2[stp::res]=p_tree_process->ScaleSetter()->Scale(stp::res);
      if (p_tree_process->ScaleSetter(1)->Amplitudes().size() &&
          p_tree_process->ScaleSetter(1)->FixedScales().empty()) {
        m_realevt.p_ampl = p_tree_process->ScaleSetter(1)->Amplitudes().front()->CopyAll();
        m_realevt.p_ampl->SetProc(this);
      }
      double real = p_tree_process->operator()(&mom.front())*p_tree_process->Norm();
      if (IsBad(real) || real == 0. ) res=false;
      m_realevt.m_me = real;
      m_realevt.m_mewgt = real;
      m_realevt.m_K = p_tree_process->LastK();
    }
  }
  m_realevt.m_trig=realtrg;
  for (size_t i(0);i<m_subevtlist.size();++i) {
    if (!m_subevtlist[i]->m_trig || !res)
      m_subevtlist[i]->m_me=m_subevtlist[i]->m_mewgt=0.0;
  }

  m_lastdxs = m_realevt.m_me;

  if (msg_LevelIsDebugging()) {
    size_t prec(msg->Precision());
    msg->SetPrecision(16);
    for (size_t k(0);k<mom.size();++k) {
      msg_Out()<<std::setw(4)<<m_flavs[k]<<": p["<<k<<"] = "<<mom[k]<<std::endl;
    }
    msg_Out()<<"dipoles (i,j;k):"<<std::endl;
    for (size_t k=0;k<m_subtermlist.size();++k) {
      msg_Out()<<std::setw(4)<<k<<" ("
               <<std::setw(2)<<m_flavs[m_subtermlist[k]->Li()]<<", "
               <<std::setw(2)<<m_flavs[m_subtermlist[k]->Lj()]<<"; "
               <<std::setw(2)<<m_flavs[m_subtermlist[k]->Lk()]<<")"
               <<", alpha="<<m_subevtlist[k]->m_alpha
               <<", me="<<m_subevtlist[k]->m_me<<std::endl;
    }
    msg_Out()<<"real:                                        "
             <<"me="<<m_realevt.m_me<<std::endl;
    msg->SetPrecision(prec);
  }
  m_mewgtinfo.m_bkw = p_tree_process->GetMEwgtinfo()->m_bkw;

  return m_lastdxs;
}

void Single_Real_Correction::FillAmplitudes(vector<METOOLS::Spin_Amplitudes>& amps,
                                            vector<vector<Complex> >& cols)
{
  p_tree_process->FillAmplitudes(amps, cols);
}

bool Single_Real_Correction::Trigger(const ATOOLS::Vec4D_Vector &p)
{
  return true;
}

size_t Single_Real_Correction::SetMCMode(const size_t mcmode)
{
  size_t cmcmode(p_tree_process->SetMCMode(mcmode));
  for (size_t i(0);i<m_subtermlist.size();++i)
    m_subtermlist[i]->SetMCMode(mcmode);
  m_mcmode=mcmode;
  return cmcmode;
}

size_t Single_Real_Correction::SetClusterMode(const size_t cmode)
{
  size_t ccmode(p_tree_process->SetClusterMode(cmode));
  for (size_t i(0);i<m_subtermlist.size();++i)
    m_subtermlist[i]->SetClusterMode(cmode);
  m_cmode=cmode;
  return ccmode;
}

void Single_Real_Correction::SetScale(const Scale_Setter_Arguments &args)
{
  if (!m_no_tree) {
    p_tree_process->SetScale(args);
    p_scale=p_tree_process->ScaleSetter();
  }
  for (size_t i(0);i<m_subtermlist.size();++i) {
    m_subtermlist[i]->SetScale(args);
  }
  for (size_t i(0);i<m_subostermlist.size();++i) {
    m_subostermlist[i]->SetScale(args);
  }
}

void Single_Real_Correction::SetKFactor(const KFactor_Setter_Arguments &args)
{
  p_tree_process->SetKFactor(args);
  for (size_t i(0);i<m_subtermlist.size();++i) {
    m_subtermlist[i]->SetKFactor(args);
  }
}


int Single_Real_Correction::NumberOfDiagrams() { 
  return m_subtermlist.size()+1;
}

Point * Single_Real_Correction::Diagram(int i) { 
  if (p_partner==this) return p_tree_process->Diagram(i); 
  return p_partner->Diagram(i);
} 

void Single_Real_Correction::AddChannels(std::list<std::string>* list) 
{ 
  if (m_pinfo.m_nlomode==nlo_mode::powheg) {
    for (size_t i(0);i<m_subtermlist.size();++i)
      m_subtermlist[i]->AddChannels(list);
  }
  p_tree_process->AddChannels(list);
}


/*------------------------------------------------------------------------------
  
  Helpers
  
  ------------------------------------------------------------------------------*/


void Single_Real_Correction::PrintProcessSummary(int it)
{
  Process_Base::PrintProcessSummary(it);
  if (p_partner!=this) {
    for(int i=0;i<it;i++) cout<<"  ";
    cout<<"  (partner process: "<<p_partner->Name()<<" *"<<m_sfactor<<")"<<endl;
//     p_partner->PrintProcessSummary(it+1);
    return;
  }
  for(int i=0;i<it+1;i++) cout<<"  ";
  cout<<"++++real term+++++++++++++++++++++++++++++"<<endl;
  p_tree_process->PrintProcessSummary(it+1);
  for(int i=0;i<it+1;i++) cout<<"  ";
  cout<<"----dipole terms--------------------------"<<endl;
  for (size_t i=0;i<m_subtermlist.size();++i) 
    if (m_subtermlist[i]->IsValid()) m_subtermlist[i]->PrintProcessSummary(it+1);
  for(int i=0;i<it+1;i++) cout<<"  ";
  cout<<"++++++++++++++++++++++++++++++++++++++++++"<<endl;
} 

void Single_Real_Correction::PrintSubevtSummary()
{
  cout<<"Subevent summary: "<<Name()<<endl;
  for (size_t i=0;i<m_subevtlist.size();++i) {
    std::cout<<m_subevtlist[i];
    for (size_t j=0;j<m_subevtlist[i]->m_n;++j)
      cout<<"Mom "<<j<<": "<<m_subevtlist[i]->p_mom[j]<<" ("<<m_subevtlist[i]->p_fl[j]<<")"<<endl; 
  }
}

void Single_Real_Correction::SetSelector(const Selector_Key &key)
{
  p_tree_process->SetSelector(key);
  for (size_t i=0;i<m_subtermlist.size();++i) {
    m_subtermlist[i]->SetSelector(key);
  }
  for (size_t i=0;i<m_subostermlist.size();++i) {
    m_subostermlist[i]->SetSelector(key);
  }
  p_selector=p_tree_process->Selector();
}

void Single_Real_Correction::SetGenerator(ME_Generator_Base *const gen) 
{ 
  if (p_tree_process==NULL) {
    p_gen=gen;
    return;
  }
  p_tree_process->SetGenerator(gen);
  for (size_t i=0;i<m_subtermlist.size();++i) {
    if (m_subtermlist[i]->GetLOProcess())
    m_subtermlist[i]->GetLOProcess()->SetGenerator(gen);
  }
  for (size_t i=0;i<m_subostermlist.size();++i) {
    m_subostermlist[i]->GetOSProcess()->SetGenerator(gen);
  }
}

void Single_Real_Correction::SetShower(PDF::Shower_Base *const ps)
{
  p_shower=ps;
  p_tree_process->SetShower(ps);
  for (size_t i=0;i<m_subtermlist.size();++i) {
    if (m_subtermlist[i]->GetLOProcess())
    m_subtermlist[i]->SetShower(ps);
  }
}

void Single_Real_Correction::SetNLOMC(PDF::NLOMC_Base *const mc)
{
  p_nlomc=mc;
  p_tree_process->SetNLOMC(mc);
  for (size_t i=0;i<m_subtermlist.size();++i) {
    if (m_subtermlist[i]->GetLOProcess())
    m_subtermlist[i]->SetNLOMC(mc);
  }
}

void Single_Real_Correction::SetFixedScale(const std::vector<double> &s)
{
  p_tree_process->SetFixedScale(s);
  for (size_t i=0;i<m_subtermlist.size();++i)
    if (m_subtermlist[i]->GetLOProcess())
    m_subtermlist[i]->GetLOProcess()->SetFixedScale(s);
  for (size_t i=0;i<m_subostermlist.size();++i)
    m_subostermlist[i]->GetOSProcess()->SetFixedScale(s);
}

void Single_Real_Correction::SetSelectorOn(const bool on)
{
  p_tree_process->SetSelectorOn(on);
  for (size_t i=0;i<m_subtermlist.size();++i)
    if (m_subtermlist[i]->GetLOProcess())
    m_subtermlist[i]->GetLOProcess()->SetSelectorOn(on);
  for (size_t i=0;i<m_subostermlist.size();++i)
    m_subostermlist[i]->GetOSProcess()->SetSelectorOn(on);
}

void Single_Real_Correction::FillProcessMap(NLOTypeStringProcessMap_Map *apmap)
{
  if (!m_no_tree) {
    Process_Base::FillProcessMap(apmap);
    p_tree_process->FillProcessMap(apmap);
  }
  for (size_t i=0;i<m_subtermlist.size();++i)
    m_subtermlist[i]->FillProcessMap(apmap);
}

ATOOLS::Flavour Single_Real_Correction::ReMap(const ATOOLS::Flavour &fl,const size_t &id) const
{
  return p_tree_process->ReMap(fl,id);
}

AMEGIC::Process_Base *AMEGIC::Single_Real_Correction::GetReal()
{
  return p_tree_process;
}

void Single_Real_Correction::SmearCounterEvents(NLO_subevtlist& subevtlist)
{
  if (m_smear_threshold==0.0 || m_subtermlist.size()==0) return;
  DEBUG_FUNC(m_smear_threshold);

  DEBUG_VAR(m_realevt.m_me);
  for (size_t i=0;i<m_subtermlist.size();i++) {
    if (!m_subtermlist[i]->IsValid()) continue;
    double alpha = m_smear_threshold>0.0 ? m_subtermlist[i]->Dipole()->KT2() :
                                           m_subtermlist[i]->Dipole()->LastAlpha();
    if (alpha<dabs(m_smear_threshold)) {
      double x=pow(alpha/dabs(m_smear_threshold), m_smear_power);

      DEBUG_VAR(alpha);
      DEBUG_VAR(x);
      DEBUG_INFO("me = "<<m_subtermlist[i]->GetSubevt()->m_me<<" --> "<<m_subtermlist[i]->GetSubevt()->m_me*x);

      m_realevt.m_me += (1.0-x)*m_subtermlist[i]->GetSubevt()->m_me;
      m_subtermlist[i]->GetSubevt()->m_me *= x;

      m_realevt.m_mewgt += (1.0-x)*m_subtermlist[i]->GetSubevt()->m_mewgt;
      m_subtermlist[i]->GetSubevt()->m_mewgt *= x;
    }
  }
  DEBUG_VAR(m_realevt.m_me);
}

bool Single_Real_Correction::AllowAsSpecInISPFF(const size_t &k)
{
  switch (m_pspisrecscheme) {
  case 0:
    if (k<m_nin) return true;
    break;
  case 1:
    if (k>=m_nin) return true;
    break;
  case 2:
    if (m_flavs[k].Charge()) return true;
    break;
  case 3:
    if (!m_flavs[k].Charge()) return true;
    break;
  case 4:
    return true;
    break;
  default:
    THROW(fatal_error,"Unknown IS P->ff recoil scheme.")
    break;
  }
  msg_Debugging()<<k<<" not allowed as spectator for recombined IS photon."
                 <<std::endl;
  return false;
}

bool Single_Real_Correction::AllowAsSpecInFSPFF(const size_t &k)
{
  switch (m_pspfsrecscheme) {
  case 0:
    if (k<m_nin) return true;
    break;
  case 1:
    if (k>=m_nin) return true;
    break;
  case 2:
    if (m_flavs[k].Charge()) return true;
    break;
  case 3:
    if (!m_flavs[k].Charge()) return true;
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

void Single_Real_Correction::SetCaller(PHASIC::Process_Base *const proc)
{
  p_caller=proc;
  p_tree_process->SetCaller(static_cast<Single_Real_Correction*>(proc)->p_tree_process);
}
