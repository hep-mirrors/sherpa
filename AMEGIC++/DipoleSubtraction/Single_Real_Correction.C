#include "AMEGIC++/DipoleSubtraction/Single_Real_Correction.H"
#include "AMEGIC++/Main/Single_Process.H"
#include "AMEGIC++/Main/Single_Process_MHV.H"
#include "AMEGIC++/DipoleSubtraction/Single_DipoleTerm.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "PDF/Main/ISR_Handler.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"

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
  p_partner(this), p_tree_process(NULL)
{
  m_Norm = 1.;  
  m_dalphamax = 1.;
  double helpd;
  Data_Reader reader(" ",";","!","=");
  if (reader.ReadFromFile(helpd,"DIPOLE_ALPHAMAX_CUT")) {
    m_dalphamax = helpd;
    msg_Tracking()<<"Set dipole cut alphamax="<<m_dalphamax<<"."<<std::endl;
  }
  rpa.gen.AddCitation(1,"The automated generation of Catani-Seymour Dipole\
 Terms is published under \\cite{Gleisberg:2007md}.");
}


Single_Real_Correction::~Single_Real_Correction()
{
  if (p_tree_process) delete p_tree_process;
  for (size_t i=0;i<m_subtermlist.size();i++) delete m_subtermlist[i];
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
    vector<ATOOLS::Flavour> flin;
    vector<ATOOLS::Flavour> flout;
    m_pinfo.m_ii.GetExternal(flin);
    m_pinfo.m_fi.GetExternal(flout);
    if (CF.MHVCalculable(flin,flout)) p_tree_process = new Single_Process_MHV();
    if (m_pinfo.m_amegicmhv==2) return 0;
  }
  if (!p_tree_process) p_tree_process = new AMEGIC::Single_Process();

  int status;

  p_tree_process->Get<PHASIC::Process_Base>()->Init(m_pinfo,p_int->Beam(),p_int->ISR());
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
    string partnerID=p_tree_process->Partner()->Name();//+"_REAL";
    for (size_t j=0;j<links.size();j++) if (Type()==links[j]->Type()) {
      if (partnerID==links[j]->Name()) {
	msg_Tracking()<<"Can map full real process: "<<Name()<<" -> "<<partnerID<<" Factor: "<<p_tree_process->GetSFactor()<<endl;
	p_mapproc = p_partner = (Single_Real_Correction*)links[j];
	m_sfactor = p_tree_process->GetSFactor();
	return 1;
      }
    }
  }

  m_realevt.n      = m_nin+m_nout;
  m_realevt.p_fl   = &(p_tree_process->Flavours().front());
  m_realevt.p_mom  = NULL;
  m_realevt.m_ID   = string("Real");

  vector<int> partlist;
  for (size_t i=0;i<m_nin+m_nout;i++) {
    if (m_flavs[i].Strong()) partlist.push_back(i);
  }
  for (size_t i=0;i<partlist.size();i++) {
    for (size_t j=0;j<partlist.size();j++) {
      for (size_t k=0;k<partlist.size();k++) if (k!=i&&k!=j&&i!=j) {
	Single_DipoleTerm *pdummy = new Single_DipoleTerm(m_pinfo,partlist[i],partlist[j],partlist[k],p_int);
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
//   SetScale(m_pinfo.m_scale,m_pinfo.m_mur2tag,m_pinfo.m_muf2tag);
//   SetKFactor(m_pinfo.m_kfactor);

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
  return p_tree_process->FillIntegrator(psh);
}


bool Single_Real_Correction::SetUpIntegrator() 
{  
  if (m_nin==2) {
    if ( (m_flavs[0].Mass() != p_int->ISR()->Flav(0).Mass()) ||
	 (m_flavs[1].Mass() != p_int->ISR()->Flav(1).Mass()) ) p_int->ISR()->SetPartonMasses(&m_flavs.front());
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


double Single_Real_Correction::Differential(const ATOOLS::Vec4D_Vector &moms) { return DSigma(moms); }

double Single_Real_Correction::Differential2() { 
  if (p_int->ISR()->On()==0) return 0.;
  return DSigma2(); 
}


double Single_Real_Correction::DSigma(const ATOOLS::Vec4D_Vector &moms)
{
  m_last = m_lastxs = 0.;
    // So far only massless partons!!!!

  if (p_partner == this) operator()(moms);
  else {
    if (m_lookup) m_lastxs = p_partner->LastXS()*m_sfactor;
    else m_lastxs = p_partner->operator()(moms)*m_sfactor;
    for (size_t i=0;i<m_subevtlist.size();++i) delete m_subevtlist[i];
    m_subevtlist.clear();
    std::vector<NLO_subevt*>* partnerlist=p_partner->GetSubevtList();
    if (partnerlist->size()==0) return 0.;
    for (size_t i=0;i<partnerlist->size();++i) {
      NLO_subevt* cpevt=new NLO_subevt((*partnerlist)[i]);
      m_subevtlist.push_back(cpevt);
    }
    m_subevtlist.Mult(m_sfactor);
  }

  if (m_nin==2) {
    double scale = p_partner->ScaleSetter()->Scale(stp::fac);
    size_t multiscale=0;
    for (size_t i=0;i<m_subevtlist.size();++i) 
      if (m_subevtlist[i]->m_scale!=scale) {
	multiscale++;
	p_int->ISR()->CalculateWeight(m_subevtlist[i]->m_scale);
	(*m_subevtlist[i])*=p_int->ISR()->Weight(&m_flavs.front());
      }
    if (multiscale<m_subevtlist.size()) {
      p_int->ISR()->CalculateWeight(scale);
      m_lastlumi = p_int->ISR()->Weight(&m_flavs.front());
      for (size_t i=0;i<m_subevtlist.size();++i) 
	if (m_subevtlist[i]->m_scale==scale) (*m_subevtlist[i])*=m_lastlumi;
    }
    int    pols[2] = {p_pl[0].type[0],p_pl[1].type[0]};
    double dofs[2] = {p_pl[0].factor[0],p_pl[1].factor[0]};
    if (p_pl[0].num>1) pols[0] = 99;
    if (p_pl[1].num>1) pols[1] = 99;
    m_subevtlist*= m_lastlumi = p_int->Beam()->Weight(pols,dofs);
  }
  for (size_t i=0;i<m_subevtlist.size();++i) m_last+=m_subevtlist[i]->m_result;
  return m_last;
}

double Single_Real_Correction::DSigma2() { 
  if ((m_flavs[0]==m_flavs[1]) || (p_int->ISR()->On()==0) ) return 0.;
  if (m_subevtlist.size()==0) return 0.;
  double scale = p_partner->ScaleSetter()->Scale(stp::fac);
  size_t multiscale=0;
  for (size_t i=0;i<m_subevtlist.size();++i) {
    if (m_subevtlist[i]->m_scale!=scale) {
      multiscale++;
      p_int->ISR()->CalculateWeight2(m_subevtlist[i]->m_scale);
      m_subevtlist[i]->m_result+=m_subevtlist[i]->m_me
	*p_int->ISR()->Weight2(&m_flavs.front());
    }
  }
  double tmp;
  if (multiscale<m_subevtlist.size()) {
    p_int->ISR()->CalculateWeight2(scale);
    tmp = p_int->ISR()->Weight2(&m_flavs.front());
    for (size_t i=0;i<m_subevtlist.size();++i) {
      if (m_subevtlist[i]->m_scale==scale) {
	m_subevtlist[i]->m_result+=m_subevtlist[i]->m_me*tmp;
      }
    }    
  }
  tmp=m_last;
  m_last=0.;
  for (size_t i=0;i<m_subevtlist.size();++i) m_last+=m_subevtlist[i]->m_result;

  return m_last-tmp;
}


double Single_Real_Correction::operator()(const ATOOLS::Vec4D_Vector &mom)
{
  m_subevtlist.clear();
  p_scale->CalculateScale(mom);
  double M2=0.;
  
  if (m_pinfo.m_fi.m_nloqcdtype&nlo_type::rsub) {
  double MS=0.;
  for (size_t i=0;i<m_subtermlist.size();i++) if (m_subtermlist[i]->IsValid()){ 
    MS+=M2=m_subtermlist[i]->operator()(&mom.front())*p_tree_process->Norm();
    if (m_pinfo.m_nlomode==0 && M2==0.0) {
      m_subevtlist.clear();
      return 0.0;
    }
    if (M2>0.||M2<0.) {
      m_subevtlist.push_back(m_subtermlist[i]->GetSubevt());
    }
  }
  if (!(MS>0.)&&!(MS<0.)&&!(MS==0.)) {
    m_subevtlist.clear();
    return 0.;
  }
  }

  if (m_pinfo.m_nlomode==1) {
  bool trg= JetTrigger(Integrator()->PSHandler()->LabPoint(),m_flavs,m_nout);
  if (trg) {
    M2 = p_tree_process->operator()(&mom.front())*KFactor();
    if (M2>0.) {
      m_realevt.p_mom  = &(Integrator()->PSHandler()->LabPoint().front());
      m_realevt.m_me   = m_realevt.m_result = M2;
      m_realevt.m_scale = ScaleSetter()->Scale(stp::fac);
      m_subevtlist.push_back(&m_realevt);
    }
  }
  }
  
  m_subevtlist.Mult(p_tree_process->Norm());
  m_lastxs = M2*p_tree_process->Norm();
  if (!(M2>0.)&&!(M2<0.)) return 0.;
  
  return m_lastxs;
}

void Single_Real_Correction::SetScale(const std::string &scale)
{
  Process_Base::SetScale(scale);
  for (size_t i(0);i<m_subtermlist.size();++i) {
    m_subtermlist[i]->SetScale(scale);
  }
}
 
void Single_Real_Correction::SetKFactor(const std::string &kfactor,
					const size_t &oqcdlo,const size_t &oewlo)
{
  Process_Base::SetKFactor(kfactor,oqcdlo,oewlo);
  for (size_t i(0);i<m_subtermlist.size();++i) {
    m_subtermlist[i]->SetKFactor(kfactor,oqcdlo,oewlo);
  }
}


int Single_Real_Correction::NumberOfDiagrams() { 
  return m_subtermlist.size()+1;
}

Point * Single_Real_Correction::Diagram(int i) { 
  if (p_partner==this) return p_tree_process->Diagram(i); 
  return p_partner->Diagram(i);
} 

void Single_Real_Correction::FillAmplitudes(HELICITIES::Amplitude_Tensor* atensor,double sfactor)
{
  if (p_partner==this) p_tree_process->FillAmplitudes(atensor,sfactor);
  else p_partner->FillAmplitudes(atensor,sfactor*sqrt(m_sfactor));
}

void Single_Real_Correction::AddChannels(std::list<std::string>* list) 
{ 
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
    m_subevtlist[i]->Print();
    for (size_t j=0;j<m_subevtlist[i]->n;++j) cout<<"Mom "<<j<<": "<<m_subevtlist[i]->p_mom[j]<<" ("<<m_subevtlist[i]->p_fl[j]<<")"<<endl; 
  }
}

void Single_Real_Correction::FillAlphaHistogram(ATOOLS::Histogram* histo,double weight)
{
  if (m_last==0.) return;
  for (size_t i=0;i<m_subevtlist.size();++i) 
    histo->InsertMCBIM(m_subevtlist[i]->m_alpha,m_subevtlist[i]->m_result*weight);
}

void Single_Real_Correction::SetSelector(const Selector_Key &key)
{
  Process_Base::SetSelector(key);
  if (m_pinfo.m_nlomode==0) {
    for (size_t i=0;i<GetSubTermNumber();++i)
      GetSubTerm(i)->GetLOProcess()->SetSelector(key);
  }
}
