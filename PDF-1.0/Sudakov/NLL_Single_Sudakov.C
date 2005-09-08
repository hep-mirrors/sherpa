#include "NLL_Single_Sudakov.H"

#include "Run_Parameter.H"
#include "Message.H"
#include "Shell_Tools.H"

using namespace SHERPA;
using namespace ATOOLS;

NLL_Single_Sudakov::
NLL_Single_Sudakov(NLL_Branching_Probability_Base * bp,int mode):
  m_calcmode((Sudakov::code)mode), 
  m_cutmode((Sudakov::code)mode), 
  p_bp(bp),
  m_gauss(bp)
{}

NLL_Single_Sudakov::~NLL_Single_Sudakov()
{
  if (p_bp) delete p_bp;
}

bool NLL_Single_Sudakov::Initialize(double _m_qmin,double _m_qmax) 
{
  m_calcmode = (Sudakov::code)(m_calcmode&896);
  m_cutmode  = (Sudakov::code)(m_cutmode&7);
  if (m_calcmode==Sudakov::table) {
    msg.Error()<<"NLL_Single_Sudakov::Initialize("
	       <<m_calcmode<<","<<m_cutmode<<"): "<<std::endl
	       <<"Table mode not implemented. Abort."<<std::endl;
    abort();
  }
  return true;
}

double NLL_Single_Sudakov::Log(double Q, double q) 
{
  if (m_calcmode == Sudakov::analytic) {
    double sum(p_bp->IntGamma(q,Q));
    if (sum!=-1.) return sum;
  }
  p_bp->SetQmax(Q);
  p_bp->SetQmin(q);
  return m_gauss.Integrate(q,Q,1.e-5,1);
}

double NLL_Single_Sudakov::operator()(double Q, double q) 
{
  if (Q<q) return 1.0;
  double sum(Log(Q,q));
  return exp(-sum);
}

double NLL_Single_Sudakov::operator()(double Q) 
{
  return operator()(Q,m_qmin);
}

double NLL_Single_Sudakov::IntGamma(double q,double Q) 
{ 
  return p_bp->IntGamma(q,Q); 
}
