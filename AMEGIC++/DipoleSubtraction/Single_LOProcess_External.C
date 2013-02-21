#include "AMEGIC++/DipoleSubtraction/Single_LOProcess_External.H"

#include "MODEL/Main/Running_AlphaS.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "AMEGIC++/Phasespace/Phase_Space_Generator.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/ISR_Handler.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"

#include <unistd.h>

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


Single_LOProcess_External::Single_LOProcess_External(const Process_Info &pi,
                                   BEAM::Beam_Spectra_Handler *const beam,
                                   PDF::ISR_Handler *const isr) :   
  Single_LOProcess(pi, beam, isr)
{
}

Single_LOProcess_External::~Single_LOProcess_External()
{
  if (p_me2) delete p_me2;
}

int AMEGIC::Single_LOProcess_External::InitAmplitude(Model_Base * model,Topology* top,
					    vector<Process_Base *> & links,
					    vector<Process_Base *> & errs)
{
  m_type=30;
  Init();
  model->GetCouplings(m_cpls);
  if (!model->CheckFlavours(m_nin,m_nout,&m_flavs.front())) return 0;
  m_newlib   = false;

  m_Norm = SymmetryFactors() * m_pol.Spin_Average(m_nin,&m_flavs.front());
  int oew(m_oew), oqcd(m_oqcd);
  m_pn=m_flavs.size();
  if (oqcd==99) oqcd=m_pn-m_oew-2;
  m_pinfo.m_fi.m_nloqcdtype=nlo_type::lo;  
  p_me2 = Tree_ME2_Base::GetME2(m_pinfo);
  if (!p_me2) return 0;
  p_me2->SetCouplings(m_cpls);
  
  m_oew=oew;
  m_oqcd=oqcd;
  std::vector<Vec4D> tmoms(p_testmoms,&p_testmoms[m_nin+m_nout]);
  m_sfactor=1.0;
  m_iresult=p_me2->Calc(tmoms);
  if (m_iresult==0.) return 0;
  for (size_t j=0;j<links.size();j++) if (Type()==links[j]->Type()) {
    if (ATOOLS::IsEqual(links[j]->Result(),Result())) {
      msg_Tracking()<<"AMEGIC::Single_Process_External::InitAmplitude : "<<std::endl
		    <<"   Found a partner for process "<<m_name<<" : "<<links[j]->Name()<<std::endl;
      p_mapproc = p_partner   = (Single_LOProcess_External*)links[j];
      for (size_t i(0);i<m_nin+m_nout;++i)
	AddtoFlavmap(ToString(1<<i),p_partner->Flavours()[i]);
      break;
    } 
  }
  if (p_partner==this) links.push_back(this);
  msg_Info()<<".";
  
  m_partonlist.clear();
  for (size_t i=0;i<m_nin;i++) if (m_flavs[i].Strong()) m_partonlist.push_back(i);
  for (size_t i=m_nin;i<m_nin+m_nout;i++) if (m_flavs[i].Strong()) m_partonlist.push_back(i);
  return 1;
}

bool AMEGIC::Single_LOProcess_External::PerformTests()
{
  return 1;
}

bool Single_LOProcess_External::SetUpIntegrator() 
{  
  return 0;
}

void Single_LOProcess_External::Minimize()
{
  if (p_partner==this) return;
  m_oqcd      = p_partner->OrderQCD();
  m_oew       = p_partner->OrderEW();
  m_ntchanmin = p_partner->NTchanMin();
}

double Single_LOProcess_External::Partonic(const ATOOLS::Vec4D_Vector& _moms,const int mode)
{
  return 0.;
}

double Single_LOProcess_External::operator()(const ATOOLS::Vec4D_Vector &labmom,const ATOOLS::Vec4D *mom,
				    std::vector<double> * pfactors,std::vector<ATOOLS::Vec4D>* epol,const int mode)
{
  if (p_partner!=this) {
    if (m_lookup) {
      m_lastxs = p_partner->LastXS()*m_sfactor;
      if (m_lastxs!=0.) return m_lastxs;
    }
    return m_lastxs = p_partner->operator()(labmom,mom,pfactors,epol,mode)*m_sfactor;
  }

  p_int->SetMomenta(labmom);
  p_scale->CalculateScale(labmom);
 
  Vec4D_Vector moms(mom,&mom[m_nin+m_nout]);
  m_lastxs = p_me2->Calc(moms);
  return m_lastxs;
}

void Single_LOProcess_External::Calc_AllXS(const ATOOLS::Vec4D_Vector &labmom,
				  const ATOOLS::Vec4D *mom,std::vector<std::vector<double> > &dsij,const int mode) 
{
  if (p_partner!=this) {
    p_partner->Calc_AllXS(labmom,mom,dsij,mode);
    dsij[0][0]*=m_sfactor;
    return;
  }
  p_int->SetMomenta(labmom);
  p_scale->CalculateScale(labmom);

  Vec4D_Vector moms(mom,&mom[m_nin+m_nout]);
  dsij[0][0]=p_me2->Calc(moms);
  for (size_t i=0;i<m_partonlist.size();i++)
    for (size_t k=i+1;k<m_partonlist.size();k++)
      dsij[i][k]=dsij[k][i]=0.0;
}

void Single_LOProcess_External::AddChannels(std::list<std::string>* tlist) 
{
}

int Single_LOProcess_External::NumberOfDiagrams()
{
  return 0;
}
