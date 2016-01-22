#include "AMEGIC++/Main/Single_Process_External.H"

#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/ISR_Handler.H"

#include "ATOOLS/Org/Run_Parameter.H"
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

AMEGIC::Single_Process_External::Single_Process_External():
  p_me2(NULL), p_partner(this)
{
  Data_Reader reader(" ",";","!","=");
  reader.AddComment("#");
  reader.SetInputPath(rpa->GetPath());
  reader.SetInputFile(rpa->gen.Variable("ME_DATA_FILE"));
  m_keep_zero_procs=reader.GetValue<int>("AMEGIC_KEEP_ZERO_PROCS",0);
  m_lastk=1.0;
}

AMEGIC::Single_Process_External::~Single_Process_External()
{
  if (p_me2) delete p_me2;
}

/*------------------------------------------------------------------------------
  
  Generic stuff for initialization of Single_Processes
      
  ------------------------------------------------------------------------------*/


int AMEGIC::Single_Process_External::InitAmplitude(Amegic_Model * model,Topology* top,
					 vector<Process_Base *> & links,
					 vector<Process_Base *> & errs)
{
  Init();
  model->p_model->GetCouplings(m_cpls);
  if (!model->p_model->CheckFlavours(m_nin,m_nout,&m_flavs.front())) return 0;
  m_newlib   = false;
  string ptypename;
  ptypename = "P"+ToString(m_nin)+"_"+ToString(m_nout);
  ATOOLS::MakeDir(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+ptypename);

  m_Norm = SymmetryFactors() * m_pol.Spin_Average(m_nin,&m_flavs.front());
  m_pn=m_flavs.size();
  p_me2 = Tree_ME2_Base::GetME2(m_pinfo);
  if (!p_me2) return 0;
  m_maxcpl[0]=m_mincpl[0]=p_me2->OrderQCD();
  m_maxcpl[1]=m_mincpl[1]=p_me2->OrderEW();
  p_me2->SetCouplings(m_cpls);
  p_me2->FillCombinations(m_ccombs,m_cflavs);
  p_me2->SetNorm(m_Norm);
  
  std::vector<Vec4D> tmoms(p_testmoms,&p_testmoms[m_nin+m_nout]);
  m_iresult=p_me2->Calc(tmoms);
  if (m_iresult==0. && !m_keep_zero_procs) return 0;
  for (size_t j=0;j<links.size();j++) if (Type()==links[j]->Type()) {
    if (m_allowmap && FlavCompare(links[j]) && ATOOLS::IsEqual(links[j]->Result(),Result())) {
      msg_Tracking()<<"AMEGIC::Single_Process_External::InitAmplitude : "<<std::endl
		    <<"   Found a partner for process "<<m_name<<" : "<<links[j]->Name()<<std::endl;
      p_mapproc = p_partner   = (Single_Process_External*)links[j];
      for (size_t i(0);i<m_nin+m_nout;++i)
	AddtoFlavmap(ToString(1<<i),p_partner->Flavours()[i]);
      break;
    } 
  }
  if (p_partner==this) links.push_back(this);
  msg_Info()<<".";
  
  if (p_partner==this && Result()>0.) SetUpIntegrator();
  return 1;
}



int AMEGIC::Single_Process_External::PerformTests()
{
  return 1;
}

/*------------------------------------------------------------------------------
  
  Phase space initialization
  
  ------------------------------------------------------------------------------*/

bool AMEGIC::Single_Process_External::FillIntegrator(PHASIC::Phase_Space_Handler *const psh)
{
  THROW(fatal_error,"No integrator");
  return Process_Base::FillIntegrator(psh);
}

bool AMEGIC::Single_Process_External::SetUpIntegrator() 
{  
  if (m_nin==2) {
    if ( (m_flavs[0].Mass() != p_int->ISR()->Flav(0).Mass()) ||
	 (m_flavs[1].Mass() != p_int->ISR()->Flav(1).Mass()) ) p_int->ISR()->SetPartonMasses(m_flavs);
  }
  return 1;
}

void AMEGIC::Single_Process_External::AddChannels(std::list<std::string>* tlist) 
{
}

/*------------------------------------------------------------------------------
  
  Process management
  
  ------------------------------------------------------------------------------*/
void AMEGIC::Single_Process_External::Minimize()
{
  if (p_partner==this) return;

  if (p_me2) {
    delete p_me2;
    p_me2=NULL;
  }

  m_maxcpl = p_partner->MaxOrders();
  m_mincpl = p_partner->MinOrders();
}

double AMEGIC::Single_Process_External::Partonic(const Vec4D_Vector &moms,const int mode) 
{ 
  if (mode==1) return m_mewgtinfo.m_B=m_lastxs;
  if (!Selector()->Result()) return m_mewgtinfo.m_B=m_lastxs = 0.0;
  if (!(IsMapped() && LookUp())) {
    p_partner->ScaleSetter()->CalculateScale(moms);
  }
  m_mewgtinfo.m_K=p_partner->LastK();
  return m_mewgtinfo.m_B=DSigma(moms,m_lookup);
}

double AMEGIC::Single_Process_External::DSigma(const ATOOLS::Vec4D_Vector &_moms,bool lookup)
{
  m_lastxs = 0.;
  if (p_partner == this) {
    m_lastxs = m_Norm * operator()((ATOOLS::Vec4D*)&_moms.front());
  }
  else {
    if (lookup && p_partner->m_lookup)
      m_lastxs = p_partner->LastXS()*m_sfactor;
    else m_lastxs = m_Norm * p_partner->operator()((ATOOLS::Vec4D*)&_moms.front())*m_sfactor;
  }
  return m_lastxs;
}

double AMEGIC::Single_Process_External::operator()(const ATOOLS::Vec4D* mom)
{
  Vec4D_Vector moms(mom,&mom[m_nin+m_nout]);
  return p_me2->Calc(moms)*(m_lastk=KFactor());
}

bool AMEGIC::Single_Process_External::Combinable
(const size_t &idi,const size_t &idj)
{
  Combination_Set::const_iterator 
    cit(m_ccombs.find(std::pair<size_t,size_t>(idi,idj)));
  return cit!=m_ccombs.end();
}

const Flavour_Vector &AMEGIC::Single_Process_External::
CombinedFlavour(const size_t &idij)
{
  CFlavVector_Map::const_iterator fit(m_cflavs.find(idij));
  if (fit==m_cflavs.end()) THROW(fatal_error,"Invalid request");
  return fit->second;
}
