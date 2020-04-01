#include "PHASIC++/Channels/Single_Channel.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/CXXFLAGS.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

Single_Channel::Single_Channel() :
  m_nin(0),m_nout(0),p_ms(NULL),
  m_rannum(0),p_rans(NULL),
  m_res1(0.),m_res2(0.),m_mres1(0.),m_mres2(0.),
  m_alpha(0.),m_alpha_save(0.),m_weight(1.),
  m_name("no_name")
{ }

Single_Channel::Single_Channel(size_t _nin,size_t _nout,const Flavour * _fl) :
  m_nin(_nin),m_nout(_nout),p_ms(new double[m_nin+m_nout+1]),
  m_rannum(0),p_rans(NULL),
  m_res1(0.),m_res2(0.),m_mres1(0.),m_mres2(0.),
  m_alpha(0.),m_alpha_save(0.),m_weight(1.),
  m_name("no_name")
{ 
  for (int i(0);i<m_nin+m_nout;i++) p_ms[i] = ATOOLS::sqr(_fl[i].Mass());
}

Single_Channel::Single_Channel(Single_Channel * old) :
  m_nin(old->m_nin),m_nout(old->m_nout),p_ms(new double[m_nin+m_nout]),
  m_rannum(old->m_rannum),p_rans(new double[m_rannum]),
  m_res1(0.),m_res2(0.),m_mres1(0.),m_mres2(0.),
  m_alpha(0.),m_alpha_save(0.),m_weight(1.),
  m_name(old->m_name)
{
  for (int i(0);i<m_nin+m_nout;i++) p_ms[i] = old->p_ms[i];
}

Single_Channel::~Single_Channel()
{
  if (p_ms)   delete[] p_ms; 
  if (p_rans) delete[] p_rans; 
}

void Single_Channel::Reset(double value) {
  m_alpha  = m_alpha_save = value;
  m_weight = 0.;
  m_res1   = m_res2 = m_mres1 = m_mres2 =0.;
}

void Single_Channel::ResetOpt() {
  m_res1 = m_res2 = m_mres1 = m_mres2 =0.;
}

void Single_Channel::AddPoint(double Value) {
}


void Single_Channel::GeneratePoint(Vec4D* p,Cut_Data * cuts)
{
  for (int i=0;i<m_rannum;i++) p_rans[i] = ran->Get();
  GeneratePoint(p,cuts,p_rans);
}

void Single_Channel::GeneratePoint(ATOOLS::Vec4D *p,Cut_Data *cuts,double *rans) 
{
  msg_Error()<<"Single_Channel::GeneratePoint(Vec4D *p,Cut_Data *cuts,double *rans): "
	     <<"Virtual Method called !"<<std::endl;
}

void Single_Channel::GenerateWeight(ATOOLS::Vec4D *p,Cut_Data *cuts) 
{
  msg_Error()<<"Single_Channel::GenerateWeight(Vec4D *p,Cut_Data *cuts): "
	     <<"Virtual Method called !"<<std::endl; 
}

void Single_Channel::GeneratePoint(const double * rns)
{
  msg_Error()<<"Single_Channel::GeneratePoint(): "
	     <<"Virtual Method called !"<<std::endl; 
}

void Single_Channel::GenerateWeight() 
{
  msg_Error()<<"Single_Channel::GenerateWeight(): "
	     <<"Virtual Method called !"<<std::endl; 
}

void Single_Channel::CalculateLimits(Info_Key &spkey,Info_Key &ykey) 
{
  msg_Error()<<"Single_Channel::CalculateLimits(..): "
 		     <<"Virtual method called!"<<std::endl;
}

void Single_Channel::CalculateLimits() 
{
  msg_Error()<<"Single_Channel::CalculateLimits(): "
 		     <<"Virtual method called!"<<std::endl;
}

int Single_Channel::ChNumber() 
{
  msg_Error()<<"Method : Single_Channel::ChNumber()"<<std::endl;
  return 0;
}

void Single_Channel::SetChNumber(int) 
{
  msg_Error()<<"Method : Single_Channel::SetChNumber()"<<std::endl;
}

const std::string Single_Channel::ChID() const
{ 
  msg_Error()<<"Virtual Method : Single_Channel::ChID()"<<std::endl;
  return std::string(""); 
}


void Single_Channel::MPISync()
{
#ifdef USING__MPI
  THROW(not_implemented,"Channel not MPI ready");
#endif
}
