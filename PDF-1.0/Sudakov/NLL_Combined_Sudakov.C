#include "NLL_Combined_Sudakov.H"
#include "Run_Parameter.H"

using namespace SHERPA;
using namespace ATOOLS;

NLL_Combined_Sudakov::NLL_Combined_Sudakov(int mode) 
{
  m_calcmode = (Sudakov::code)(mode&896);
  m_cutmode  = (Sudakov::code)(mode&7);

  m_qlimit = rpa.gen.Ecms() * sqrt(1./3.); 
}

void NLL_Combined_Sudakov::Add(NLL_Single_Sudakov* sud) 
{
  m_suds.push_back(sud);
}

double NLL_Combined_Sudakov::Log(double Q, double q) 
{
  double sum=0.;
  for (size_t i=0;i<m_suds.size();++i) {
    sum+=m_suds[i]->Log(Q,q);
  }
  return sum;
}


double NLL_Combined_Sudakov::operator()(double Q, double q) 
{
  if ((Q>m_qlimit)&&(m_cutmode&Sudakov::cutatkinlim)) Q=m_qlimit;
  if (Q<q) return 1.;

  double sum = Log(Q,q);
  if ((sum<0.)&&(m_cutmode&Sudakov::cutatone)) return 1.;
  return exp(-sum);
  /*
  if (m_calcmode == Sudakov::analytic) {
    double sum=0.;
    int  hit=0;
    for (int i=0;i<m_suds.size();++i) {
      double psum=m_suds[i]->IntGamma(q,Q);
      sum+=psum;
      if (psum==-1.) hit=1;
    }
    if (hit==0) {
      if ((sum<0.)&&(m_cutmode&Sudakov::cutatone)) return 1.;
      return exp(-sum);
    }
  }
  double sum=0.;
  for (int i=0;i<m_suds.size();++i) {
    m_suds[i]->SetQmax(Q);
    m_suds[i]->SetQmin(q);
    sum+=m_suds[i]->GetIntegrator().Integrate(q,Q,1.e-5,1);
  }
  if ((sum<0.)&&(m_cutmode&Sudakov::cutatone)) return 1.;
  return exp(-sum);
  */
}

double NLL_Combined_Sudakov::operator()(double Q) {
  return operator()(Q,m_qmin);
}

NLL_Combined_Sudakov::~NLL_Combined_Sudakov()
{
  for (size_t i=0;i<m_suds.size();++i)
    delete m_suds[i];
  m_suds.clear();
}
