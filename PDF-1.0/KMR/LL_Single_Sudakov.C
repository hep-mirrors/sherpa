#include "LL_Single_Sudakov.H"

#include "Run_Parameter.H"
#include "Message.H"
#include "Integration_Info.H"
#include <sys/stat.h>

using namespace PDF;

LL_Single_Sudakov::
LL_Single_Sudakov(SHERPA::NLL_Branching_Probability_Base *bp,int mode):
  m_qmax(0.),
  m_qmin(0.),
  m_calcmode((SHERPA::Sudakov::code)mode),
  m_cutmode((SHERPA::Sudakov::code)mode),
  p_bp(bp),
  m_gauss(bp) {}

LL_Single_Sudakov::~LL_Single_Sudakov()
{
  delete p_bp;
}

void LL_Single_Sudakov::AssignKey(ATOOLS::Integration_Info *const info)
{
  m_lastsud.Assign(std::string("Sudakov_")+p_bp->Name(),3,0,info);
}

bool LL_Single_Sudakov::Initialize(double qmin,double qmax) 
{
  m_qmin=qmin;
  m_qmax=qmax;
  m_calcmode=(SHERPA::Sudakov::code)(m_calcmode&896);
  m_cutmode=(SHERPA::Sudakov::code)(m_cutmode&7);
  if (m_qmin==m_qmax) {
    m_qmin=sqrt(ATOOLS::rpa.gen.Ycut())*ATOOLS::rpa.gen.Ecms();
    m_qmax=2.*ATOOLS::rpa.gen.Ecms();
  }
  return true;
}

double LL_Single_Sudakov::Log(double Q, double q) 
{
  if (q==m_lastsud[0] && Q==m_lastsud[1] && 
      m_lastsud.Status()==ATOOLS::si::diced) return m_lastsud[2];
  p_bp->SetQmax(m_lastsud[1]=Q);
  p_bp->SetQmin(m_lastsud[0]=q);
  m_lastsud.SetStatus(ATOOLS::si::diced);
  m_lastsud[2]=m_gauss.Integrate(q,Q,1.e-5,1);
  return m_lastsud[2];
}

double LL_Single_Sudakov::operator()(double Q, double q) 
{
  return exp(-Log(Q,q));
}

double LL_Single_Sudakov::operator()(double Q) 
{
  ATOOLS::msg.Error()<<"LL_Single_Sudakov::operator("<<Q<<"): "
		     <<"Did not set m_qmin yet !"<<std::endl;
  return 1.;
}

