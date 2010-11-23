#include "AMEGIC++/DipoleSubtraction/Single_Real_Correction.H"
#include "AMEGIC++/Main/Single_Process.H"
#include "AMEGIC++/Main/Single_Process_MHV.H"
#include "AMEGIC++/Main/Single_Process_External.H"
#include "AMEGIC++/DipoleSubtraction/Single_DipoleTerm.H"
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
  p_partner(this), p_tree_process(NULL)
{
  m_Norm = 1.;  
  static bool addcite(false);
  if (!addcite) {
    addcite=true;
  rpa.gen.AddCitation(1,"The automated generation of Catani-Seymour dipole\
 terms in Amegic is published under \\cite{Gleisberg:2007md}.");
  }
}


Single_Real_Correction::~Single_Real_Correction()
{
  p_selector=NULL;
  if (p_tree_process) delete p_tree_process;
  for (size_t i=0;i<m_subtermlist.size();i++) delete m_subtermlist[i];
  for (std::map<void*,ATOOLS::Flavour_Vector*>::const_iterator it(m_dfmap.begin());
       it!=m_dfmap.end();++it) delete it->second;
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

  int status;

  Process_Info rinfo(m_pinfo);
  rinfo.m_fi.m_nloqcdtype&=(nlo_type::code)~nlo_type::rsub;
  rinfo.m_fi.m_nloewtype&=(nlo_type::code)~nlo_type::rsub;
  p_tree_process->PHASIC::Process_Base::Init(rinfo,p_int->Beam(),p_int->ISR());
  p_tree_process->SetTestMoms(p_testmoms);

  status = p_tree_process->InitAmplitude(model,top,links,errs);

  SetOrderQCD(p_tree_process->OrderQCD());
  SetOrderEW(p_tree_process->OrderEW());
  if (p_tree_process->NewLibs()) m_newlib = 1;

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
	return 1;
      }
    }
  }

  Subprocess_Info info(m_pinfo.m_ii);
  info.Add(m_pinfo.m_fi);
  m_decinfos=info.GetDecayInfos();

  m_real_momenta.resize(m_nin+m_nout);

  m_realevt.m_n    = m_nin+m_nout;
  m_realevt.p_fl   = &Flavours().front();
  m_realevt.p_dec  = &m_decinfos;

  m_realevt.p_mom  = &m_real_momenta.front();
  m_realevt.m_i = m_realevt.m_j = m_realevt.m_k = 0;

  m_sids.resize(m_nin+m_nout);
  for (size_t i(0);i<m_nin+m_nout;++i) m_sids[i]=1<<i;
  m_realevt.p_id=&m_sids.front();
  m_realevt.m_pname = GenerateName(m_pinfo.m_ii,m_pinfo.m_fi);
  m_realevt.m_pname = m_realevt.m_pname.substr(0,m_realevt.m_pname.rfind("__"));
  m_realevt.p_proc = p_tree_process->Integrator();

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
            if (pdummy->NewLibs()) m_newlib = 1;
            m_subtermlist.push_back(pdummy);
          }
          else delete pdummy;
	}
	else delete pdummy;
      }
    }
  }

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
}

void Single_Real_Correction::ReMapFlavs(NLO_subevt *const sub)
{
  std::map<void*,Flavour_Vector*>::const_iterator it(m_dfmap.find((void*)sub->p_fl));
  if (it!=m_dfmap.end()) {
    sub->p_fl=&it->second->front();
    sub->m_pname=m_dnmap[(void*)sub->p_fl];
    return;
  }
  Flavour_Vector *fls(new Flavour_Vector());
  for (size_t i(0);i<sub->m_n;++i)
    fls->push_back(p_tree_process->ReMap(sub->p_fl[i],ToString(sub->p_id[i])));
  m_dfmap[(void*)sub->p_fl]=fls;
  std::string name(sub->m_pname);
  for (size_t pos(0), i(0);i<fls->size();++i) {
    std::string ofln(ToString(sub->p_fl[i])), nfln(ToString((*fls)[i]));
    pos=name.find(ofln,pos);
    name.replace(pos,ofln.length(),nfln);
    pos+=nfln.length();
  }
  sub->p_fl=&fls->front();
  sub->m_pname=m_dnmap[(void*)sub->p_fl]=name;
}

double Single_Real_Correction::Partonic(const ATOOLS::Vec4D_Vector &moms,const int mode)
{
  if (mode==1 && !p_partner->p_tree_process->ScaleSetter()->Scale2()) return m_lastxs;
  m_lastxs = 0.;
    // So far only massless partons!!!!

  // fill m_subevtlist
  if (p_partner == this) m_lastdxs = operator()(moms,mode);
  else {
    if (m_lookup) m_lastdxs = p_partner->LastDXS()*m_sfactor;
    else m_lastdxs = p_partner->operator()(moms,mode)*m_sfactor;
    for (size_t i=0;i<m_subevtlist.size();++i) delete m_subevtlist[i];
    m_subevtlist.clear();
    std::vector<NLO_subevt*>* partnerlist=p_partner->GetSubevtList();
    if (partnerlist->size()==0) return 0.;
    for (size_t i=0;i<partnerlist->size();++i) {
      NLO_subevt* cpevt=new NLO_subevt((*partnerlist)[i]);
      ReMapFlavs(cpevt);
      m_subevtlist.push_back(cpevt);
    }
    for (size_t i=0;i<m_subevtlist.size()-1;++i)
      m_subevtlist[i]->p_real=m_subevtlist.back();
    m_subevtlist.Mult(m_sfactor);
  }
  return m_lastxs=m_lastdxs;
}

double Single_Real_Correction::operator()(const ATOOLS::Vec4D_Vector &_mom,const int mode)
{
  m_subevtlist.clear();
  p_tree_process->Integrator()->SetMomenta(_mom);
  for (size_t i=0; i<m_real_momenta.size(); ++i) m_real_momenta[i]=_mom[i];

  Vec4D_Vector mom(_mom);
  if (m_nin==2 && p_int->ISR() && p_int->ISR()->On()) {
    Poincare cms(mom[0]+mom[1]);
    for (size_t i(0);i<mom.size();++i) cms.Boost(mom[i]);
  }

  bool res=true;
  for (size_t i=0;i<m_subtermlist.size();i++) if (m_subtermlist[i]->IsValid()){
    if (IsBad((*m_subtermlist[i])(&mom.front(),mode))) res=false;
    m_subevtlist.push_back(m_subtermlist[i]->GetSubevt());
  }

  m_subevtlist.push_back(&m_realevt);
  m_realevt.p_mom  = &p_int->Momenta().front();
  m_realevt.m_me   = m_realevt.m_result =
    m_realevt.m_last[0] = m_realevt.m_last[1] = 0.0;

  bool trg(false);
  trg=p_tree_process->Selector()->JetTrigger(_mom,&m_subevtlist);
  trg|=!p_tree_process->Selector()->On();

  for (NLO_subevtlist::const_iterator sit(m_subevtlist.begin());
       sit!=--m_subevtlist.end();++sit) (*sit)->p_real=&m_realevt;

  if (trg) {
    m_realevt.m_muf2=p_tree_process->ScaleSetter()->CalculateScale(_mom,mode);
    m_realevt.m_mur2=p_tree_process->ScaleSetter()->Scale(stp::ren);
    m_realevt.m_me = m_realevt.m_mewgt
      = p_tree_process->operator()(&mom.front());
    if (IsBad(m_realevt.m_me)) res=false;
  }

  if (res) m_subevtlist.Mult(p_tree_process->Norm());
  else {
    for (size_t i(0);i<m_subevtlist.size();++i)
      m_subevtlist[i]->m_result=m_subevtlist[i]->m_last[0]=
        m_subevtlist[i]->m_last[1]=m_subevtlist[i]->m_me=0.0;
  }
  m_lastdxs = m_realevt.m_me;
  return m_lastdxs;
}

bool Single_Real_Correction::Trigger(const ATOOLS::Vec4D_Vector &p)
{
//   if (p_tree_process->IsMapped() && p_tree_process->LookUp())
//     return p_tree_process->Selector()->Result();
  return p_tree_process->Selector()->NoJetTrigger(p);
}

void Single_Real_Correction::SetScale(const Scale_Setter_Arguments &args)
{
  p_tree_process->SetScale(args);
  for (size_t i(0);i<m_subtermlist.size();++i) {
    m_subtermlist[i]->SetScale(args);
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

void Single_Real_Correction::FillAmplitudes(METOOLS::Amplitude_Tensor* atensor,double sfactor)
{
  if (p_partner==this) p_tree_process->FillAmplitudes(atensor,sfactor);
  else p_partner->FillAmplitudes(atensor,sfactor*sqrt(m_sfactor));
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
  p_selector=p_tree_process->Selector();
}

void Single_Real_Correction::SetGenerator(ME_Generator_Base *const gen) 
{ 
  p_tree_process->SetGenerator(gen);
  for (size_t i=0;i<m_subtermlist.size();++i) {
    m_subtermlist[i]->GetLOProcess()->SetGenerator(gen);
  }
}

void Single_Real_Correction::SetShower(PDF::Shower_Base *const ps)
{
  p_tree_process->SetShower(ps);
  for (size_t i=0;i<m_subtermlist.size();++i) {
    m_subtermlist[i]->GetLOProcess()->SetShower(ps);
  }
}

void Single_Real_Correction::SetFixedScale(const std::vector<double> &s)
{
  p_tree_process->SetFixedScale(s);
  for (size_t i=0;i<m_subtermlist.size();++i)
    m_subtermlist[i]->GetLOProcess()->SetFixedScale(s);
}

void Single_Real_Correction::SetSelectorOn(const bool on)
{
  p_tree_process->SetSelectorOn(on);
  for (size_t i=0;i<m_subtermlist.size();++i)
    m_subtermlist[i]->GetLOProcess()->SetSelectorOn(on);
}

