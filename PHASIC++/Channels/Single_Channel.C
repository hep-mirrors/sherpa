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
  m_name("no_name"),m_incisr(0),m_status(1)
{ }

Single_Channel::Single_Channel(size_t _nin,size_t _nout,const Flavour * _fl) :
  m_nin(_nin),m_nout(_nout),p_ms(new double[m_nin+m_nout+1]),
  m_rannum(0),p_rans(NULL),
  m_res1(0.),m_res2(0.),m_mres1(0.),m_mres2(0.),
  m_alpha(0.),m_alpha_save(0.),m_weight(1.),
  m_name("no_name"),m_incisr(0),m_status(1)
{ 
  for (int i(0);i<m_nin+m_nout;i++) p_ms[i] = ATOOLS::sqr(_fl[i].Mass());
}

Single_Channel::Single_Channel(Single_Channel * old) :
  m_nin(old->m_nin),m_nout(old->m_nout),p_ms(new double[m_nin+m_nout]),
  m_rannum(old->m_rannum),p_rans(new double[m_rannum]),
  m_res1(0.),m_res2(0.),m_mres1(0.),m_mres2(0.),
  m_alpha(0.),m_alpha_save(0.),m_weight(1.),
  m_name(old->m_name),m_incisr(0),m_status(1)
{
  for (int i(0);i<m_nin+m_nout;i++) p_ms[i] = old->p_ms[i];
  m_incisr=old->IncludesISR();
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
  for (size_t i=0;i<m_rannum;i++) p_rans[i] = ran->Get();
  GeneratePoint(p,cuts,p_rans);
}

void Single_Channel::GeneratePoint(ATOOLS::Vec4D *p,Cut_Data *cuts,double *rans) 
{
  msg_Error()<<"Single_Channel::GeneratePoint(Vec4D *p,Cut_Data *cuts,double *rans): "
	     <<"Virtual Method called !"<<std::endl;
}

void Single_Channel::GenerateWeight(ATOOLS::Vec4D *p,Cut_Data *cuts,bool recompute)
{
  msg_Error()<<"Single_Channel::GenerateWeight(Vec4D *p,Cut_Data *cuts): "
	     <<"Virtual Method called !"<<std::endl; 
}

void Single_Channel::GeneratePoint(const double * rns)
{
  msg_Error()<<"Single_Channel::GeneratePoint(): "
	     <<"Virtual Method called !"<<std::endl; 
}

void Single_Channel::GenerateWeight(const int & mode) 
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

void Single_Channel::MPISync()
{
#ifdef USING__MPI
  THROW(not_implemented,"Channel not MPI ready");
#endif
}

void Single_Channel::Optimize()              {}
void Single_Channel::EndOptimize()           {}
void Single_Channel::WriteOut(std::string)   {}
void Single_Channel::ReadIn(std::string)     {}
bool Single_Channel::OptimizationFinished()  { return false; }

void Single_Channel::ISRInfo(int & type,double & mass,double & width) {
  type = 0; mass = width = 0.0;
}
void Single_Channel::ISRInfo(std::vector<int> &ts,
	     std::vector<double> &ms,std::vector<double> &ws) const {}

void Single_Channel::MPISync();

void Single_Channel::CopyMPIValues() {
  m_res1 += m_mres1;
  m_res2 += m_mres2;
  m_mres1 = m_mres2 = 0.0;
  m_status = 1;
}

size_t Single_Channel::NChannels() { return 1; }
int Single_Channel::OType() { return 0; }
