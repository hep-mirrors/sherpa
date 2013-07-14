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
#include "ATOOLS/Org/Data_Reader.H"

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
  m_ossubon(0), p_partner(this), p_tree_process(NULL)
{
  m_Norm = 1.;  
  static bool addcite(false);
  if (!addcite) {
    addcite=true;
  rpa->gen.AddCitation(1,"The automated generation of Catani-Seymour dipole\
 terms in Amegic is published under \\cite{Gleisberg:2007md}.");
  }
  int helpi;
  Data_Reader reader(" ",";","!","=");
  reader.AddComment("#");
  reader.SetInputPath(rpa->GetPath());
  reader.SetInputFile(rpa->gen.Variable("ME_DATA_FILE"));
  if (reader.ReadFromFile(helpi,"OS_SUB")) {
    m_ossubon = helpi;
    if (m_ossubon==1) msg_Tracking()<<"Set on shell subtraction on. "<<std::endl;
  }
  m_smear_threshold=reader.GetValue<double>("NLO_SMEAR_THRESHOLD",0.0);
  m_smear_power=reader.GetValue<double>("NLO_SMEAR_POWER",4.0);
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

int Single_Real_Correction::InitAmplitude(Model_Base * model,Topology* top,
					vector<Process_Base *> & links,
					vector<Process_Base *> & errs)
{
  Init();
  if (!model->CheckFlavours(m_nin,m_nout,&m_flavs.front())) return 0;

  m_valid    = true;
  m_newlib   = false;
//   m_name+= "_REAL";
  if (m_pinfo.m_amegicmhv>0) {
    if (m_pinfo.m_amegicmhv==10) p_tree_process = new Single_Process_External();
    else if (CF.MHVCalculable(m_pinfo)) p_tree_process = new Single_Process_MHV();
    if (m_pinfo.m_amegicmhv==2) return 0;
  }
  if (!p_tree_process) p_tree_process = new AMEGIC::Single_Process();

  p_tree_process->SetSubevtList(&m_subevtlist);

  int status;

  Process_Info rinfo(m_pinfo);
  rinfo.m_fi.m_nloqcdtype&=(nlo_type::code)~nlo_type::rsub;
  rinfo.m_fi.m_nloewtype&=(nlo_type::code)~nlo_type::rsub;
  p_tree_process->PHASIC::Process_Base::Init(rinfo,p_int->Beam(),p_int->ISR());
  p_tree_process->SetTestMoms(p_testmoms);

  status = p_tree_process->InitAmplitude(model,top,links,errs);

  SetOrderQCD(p_tree_process->OrderQCD());
  SetOrderEW(p_tree_process->OrderEW());
  if (p_tree_process->Partner()->NewLibs()) m_newlib = 1;

  m_iresult=p_tree_process->Result();
  if (status==0) {
    return status;
  }

  if (p_tree_process!=p_tree_process->Partner()) {
    string partnerID=p_tree_process->Partner()->Name();
    partnerID.erase(partnerID.find("("),3);
    for (size_t j=0;j<links.size();j++) if (Type()==links[j]->Type()) {
      string lname=links[j]->Name();
      lname.erase(lname.find("("),4);
      if (partnerID==lname) {
	msg_Tracking()<<"Can map full real process: "<<Name()<<" -> "<<partnerID<<" Factor: "<<p_tree_process->GetSFactor()<<endl;
	p_mapproc = p_partner = (Single_Real_Correction*)links[j];
	m_sfactor = p_tree_process->GetSFactor();
	// return 1;
      }
    }
  }

  m_real_momenta.resize(m_nin+m_nout);

  m_realevt.m_n    = m_nin+m_nout;
  m_realevt.p_fl   = &p_tree_process->Flavours().front();
  m_realevt.p_dec  = &m_decins;

  m_realevt.p_mom  = &m_real_momenta.front();
  m_realevt.m_i = m_realevt.m_j = m_realevt.m_k = 0;

  m_sids.resize(m_nin+m_nout);
  for (size_t i(0);i<m_nin+m_nout;++i) m_sids[i]=1<<i;
  m_realevt.p_id=&m_sids.front();
  m_realevt.m_pname = GenerateName(m_pinfo.m_ii,m_pinfo.m_fi);
  m_realevt.m_pname = m_realevt.m_pname.substr(0,m_realevt.m_pname.rfind("__"));
  m_realevt.p_proc = this;
  m_realevt.p_real = &m_realevt;

  Process_Info cinfo(m_pinfo);
  cinfo.m_fi.m_nloqcdtype&=(nlo_type::code)~nlo_type::real;
  cinfo.m_fi.m_nloewtype&=(nlo_type::code)~nlo_type::real;
  vector<int> partlist;
  for (size_t i=0;i<m_nin+m_nout;i++) {
    if (m_flavs[i].Strong()) partlist.push_back(i);
  }
  for (size_t i=0;i<partlist.size();i++) {
    for (size_t j=0;j<partlist.size();j++) {
      for (size_t k=0;k<partlist.size();k++) if (k!=i&&k!=j&&i!=j) {
	Single_DipoleTerm *pdummy = new Single_DipoleTerm(cinfo,partlist[i],partlist[j],partlist[k],p_int);
	if (pdummy->IsValid()) {
          pdummy->SetTestMoms(p_testmoms);
          int st=pdummy->InitAmplitude(model,top,links,errs);
          if (pdummy->IsValid()) {
            status=Min(st,status);
            if (pdummy->Partner()->NewLibs()) m_newlib = 1;
            m_subtermlist.push_back(pdummy);
            m_subtermlist.back()->SetNorm(p_tree_process->Norm());
	    m_subtermlist.back()->SetSmearThreshold(m_smear_threshold);
          }
          else delete pdummy;
	}
	else delete pdummy;
      }
    }
  }
  if (m_ossubon){
    Process_Info sinfo(m_pinfo);
    sinfo.m_fi.m_nloqcdtype=nlo_type::lo;
    sinfo.m_fi.m_nloewtype=nlo_type::lo;
    for (size_t i=0;i<m_flavs.size();i++) if (m_flavs[i].IsSusy()){
      for (size_t j=0;j<partlist.size();j++) if (i!=partlist[j]) {
        for (size_t swit=0;swit<5;swit++) {
  	  Single_OSTerm *pdummy = new Single_OSTerm(sinfo,i,partlist[j],swit,p_int);
	  if (pdummy->IsValid()) {
            pdummy->SetTestMoms(p_testmoms);
            int st=pdummy->InitAmplitude(model,top,links,errs);
            if (pdummy->IsValid()) {
              status=Min(st,status);
              if (pdummy->NewLibs()) m_newlib = 1;
              m_subostermlist.push_back(pdummy);
              m_subostermlist.back()->SetNorm(p_tree_process->Norm());
            }
            else delete pdummy;
	  }
	  else delete pdummy;
        }
      }
    }
  }

  if (p_mapproc && !p_partner->NewLibs()) Minimize();

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
  if (!SetUpIntegrator()) THROW(fatal_error,"No integrator");
  if (m_pinfo.m_nlomode==2) return true;
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
  if (p_tree_process) p_tree_process->SetLookUp(lookup);
  for (size_t i=0;i<m_subtermlist.size();i++) 
    m_subtermlist[i]->SetLookUp(lookup);
}

void Single_Real_Correction::Minimize()
{
//   p_tree_process->Minimize();
//   for (size_t i=0;i<m_subtermlist.size();i++) 
//     m_subtermlist[i]->Minimize();
//   for (size_t i=0;i<m_subostermlist.size();i++) 
//     m_subostermlist[i]->Minimize();
  for (size_t i=0;i<m_subevtlist.size();++i) delete m_subevtlist[i];
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
  for (size_t i=0;i<m_subevtlist.size()-1;++i)
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

double Single_Real_Correction::Partonic(const ATOOLS::Vec4D_Vector &moms,const int mode)
{
  if (mode==1) return m_lastxs;
  m_lastxs = 0.;
    // So far only massless partons!!!!

  // fill m_subevtlist
  if (p_partner == this) m_lastdxs = operator()(moms,mode);
  else {
    if (m_lookup) m_lastdxs = p_partner->LastDXS()*m_sfactor;
    else m_lastdxs = p_partner->operator()(moms,mode)*m_sfactor;
    std::vector<NLO_subevt*>* partnerlist=p_partner->GetSubevtList();
    if (partnerlist->size()!=m_subevtlist.size()) THROW(fatal_error,"Internal error");
    for (size_t i=0;i<partnerlist->size();++i) {
      m_subevtlist[i]->CopyXSData((*partnerlist)[i]);
      if (m_subevtlist[i]->p_ampl) m_subevtlist[i]->p_ampl->Delete();
      m_subevtlist[i]->p_ampl=NULL;
      if ((*partnerlist)[i]->p_ampl) {
	m_subevtlist[i]->p_ampl = (*partnerlist)[i]->p_ampl->CopyAll();
	for (Cluster_Amplitude *campl(m_subevtlist[i]->p_ampl);campl;campl=campl->Next()) {
	  for (size_t i(0);i<campl->Legs().size();++i)
	    campl->Leg(i)->SetFlav(ReMap(campl->Leg(i)->Flav(),campl->Leg(i)->Id()));
	}
      }
    }
    m_subevtlist.Mult(m_sfactor);
  }
  DEBUG_VAR(m_lastdxs);
  return m_lastxs=m_lastdxs;
}

double Single_Real_Correction::operator()(const ATOOLS::Vec4D_Vector &_mom,const int mode)
{
  m_subevtlist.clear();
  p_tree_process->Integrator()->SetMomenta(_mom);
  for (size_t i=0; i<m_real_momenta.size(); ++i) m_real_momenta[i]=_mom[i];

  Vec4D_Vector mom(_mom);
  Poincare cms;
  if (m_nin==2 && p_int->ISR() && p_int->ISR()->On()) {
    cms=Poincare(mom[0]+mom[1]);
    for (size_t i(0);i<mom.size();++i) cms.Boost(mom[i]);
  }

  bool res=true;
  for (size_t i=0;i<m_subtermlist.size();i++) if (m_subtermlist[i]->IsValid()){
    m_subtermlist[i]->Integrator()->SetMomenta(_mom);
    double test = (*m_subtermlist[i])(&mom.front(),cms,mode|2);
    if (IsBad(test)) res=false;
    m_subevtlist.push_back(m_subtermlist[i]->GetSubevt());
    m_subevtlist.back()->p_real=&m_realevt;
  }

  if (m_ossubon){
    for (size_t i=0;i<m_subostermlist.size();i++) if (m_subostermlist[i]->IsValid()){
      double test = (*m_subostermlist[i])(&mom.front(),cms,mode|2);
      if (IsBad(test)) res=false;
      m_subevtlist.push_back(m_subostermlist[i]->GetSubevt());
      m_subevtlist.back()->p_real=&m_realevt;
    }
  }

  m_subevtlist.push_back(&m_realevt);
  m_realevt.m_me = m_realevt.m_mewgt = 0.0;
  m_realevt.m_trig = false;

  bool trg(false);

  if (res) {
  m_realevt.m_mu2[stp::fac]=p_tree_process->ScaleSetter()->CalculateScale(_mom,m_cmode);
  m_realevt.m_mu2[stp::ren]=p_tree_process->ScaleSetter()->Scale(stp::ren);
  m_realevt.m_mu2[stp::res]=p_tree_process->ScaleSetter()->Scale(stp::res);
  }
  if (m_realevt.p_ampl) m_realevt.p_ampl->Delete();
  m_realevt.p_ampl=NULL;
  if (p_tree_process->ScaleSetter(1)->Amplitudes().size() &&
      p_tree_process->ScaleSetter(1)->FixedScales().empty())
    m_realevt.p_ampl = p_tree_process->ScaleSetter(1)->Amplitudes().front()->CopyAll();
  trg=p_tree_process->Selector()->JetTrigger(_mom,&m_subevtlist);
  trg|=!p_tree_process->Selector()->On();
  m_realevt.m_trig=trg;
  if (m_smear_threshold!=0.0) SmearCounterEvents(m_subevtlist);
  if (res && trg) {
    double real = p_tree_process->operator()(&mom.front())*p_tree_process->Norm();
    m_realevt.m_me += real;
    m_realevt.m_mewgt += real;
    if (IsBad(m_realevt.m_me)) res=false;
  }
  for (size_t i(0);i<m_subevtlist.size();++i) {
    if (m_subevtlist[i]->m_trig==false || !res)
      m_subevtlist[i]->m_me=m_subevtlist[i]->m_mewgt=0.0;
  }

  m_lastdxs = m_realevt.m_me;

  if (msg_LevelIsDebugging()) {
    msg->SetPrecision(16);
    msg_Out() << "// Single_Real_Correction for " << Name() << endl;
    msg_Out() << "std::vector<Vec4D> p(" << mom.size() << ");" << endl;
    for (size_t k=0;k<mom.size();++k) {
      msg_Out() << "p[" << k << "] = Vec4D" << mom[k] << "; // " << m_flavs[k] << endl;
    }
    msg_Out() << "double realme = " << m_realevt.m_me << ";" << endl;
    msg_Out() << "std::vector<double> dipoles(" << m_subtermlist.size() << ");" << endl;
    for (size_t k=0;k<m_subtermlist.size();++k) {
      msg_Out() << "dipoles[" << k << "] = " << m_subtermlist[k]->LastXS() << "; // ("
                << m_flavs[m_subtermlist[k]->Li()] << ", "
                << m_flavs[m_subtermlist[k]->Lj()] << "; "
                << m_flavs[m_subtermlist[k]->Lk()] << ")" << endl;
    }
    msg->SetPrecision(6);
  }

  return m_lastdxs;
}

void Single_Real_Correction::FillAmplitudes(vector<METOOLS::Spin_Amplitudes>& amps,
                                            vector<vector<Complex> >& cols)
{
  p_tree_process->FillAmplitudes(amps, cols);
}

bool Single_Real_Correction::Trigger(const ATOOLS::Vec4D_Vector &p)
{
//   if (p_tree_process->IsMapped() && p_tree_process->LookUp())
//     return p_tree_process->Selector()->Result();
  return p_tree_process->Selector()->NoJetTrigger(p);
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
  p_tree_process->SetScale(args);
  p_scale=p_tree_process->ScaleSetter();
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
  if (m_pinfo.m_nlomode==2) {
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
