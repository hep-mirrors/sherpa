#include "NLL_Single_Sudakov.H"
#include "Run_Parameter.H"

using namespace SHERPA;
using namespace ATOOLS;

NLL_Single_Sudakov::NLL_Single_Sudakov(NLL_Branching_Probability_Base * bp,int mode):
  p_bp(bp),m_gauss(bp) 
{
  m_calcmode = (Sudakov::code)(mode&896);
  m_cutmode  = (Sudakov::code)(mode&7);

  m_qmin=sqrt(rpa.gen.Ycut())*rpa.gen.Ecms();
  m_qmax=2.*rpa.gen.Ecms();

  m_qlimit = rpa.gen.Ecms() * sqrt(1./3.); 
  if (m_calcmode==Sudakov::table) {
    // initialize table
    std::cout<<" Initializing table "<<std::endl;
    m_calcmode=Sudakov::create_table;
    log_delta.Init(*this,m_qmin,m_qmax,600);
    m_calcmode=Sudakov::table;
    std::cout<<" done "<<std::endl;

    //    log_delta.WriteOut("Sudakov-test.dat");
  }
}

double NLL_Single_Sudakov::Log(double Q, double q) 
{
  // using m_qmin

  if (m_calcmode == Sudakov::analytic) {
    double sum=p_bp->IntGamma(q,Q);
    if (sum!=-1.) {
      return sum;
    }
    else {
      std::cout<<" WARNING integrated branching prob not analytically calculated yet! "<<std::endl;
      std::cout<<"      try numeric integration "<<std::endl;
    }
  }
  else if (m_calcmode == Sudakov::table) {
    if (q!=m_qmin) {
      std::cout<<"ERROR: table calculated with qmin="<<m_qmin<<" but called with "<<q<<std::endl;
      m_qmin=q;
      // reinit
      std::cout<<" Initializing table "<<std::endl;
      m_calcmode=Sudakov::create_table;
      log_delta.Init(*this,m_qmin,m_qmax,600);
      m_calcmode=Sudakov::table;
      std::cout<<" done "<<std::endl;
    }
    return log_delta(Q);
  }

  // m_calcmode == Sudakov::numeric || m_calcmode == Sudakov::create_table
  p_bp->SetQmax(Q);
  p_bp->SetQmin(q);
  return m_gauss.Integrate(q,Q,1.e-5,1);
}

double NLL_Single_Sudakov::operator()(double Q, double q) 
{
  if ((Q>m_qlimit)&&(m_cutmode&Sudakov::cutatkinlim)) Q=m_qlimit;
  if (Q<q) return 1.;

  double sum=Log(Q,q);
  if (m_calcmode==Sudakov::create_table) return sum;
  if ((sum<0.)&&(m_cutmode&Sudakov::cutatone)) return 1.;
  return exp(-sum);
}

double NLL_Single_Sudakov::operator()(double Q) {
  return operator()(Q,m_qmin);
}

NLL_Single_Sudakov::~NLL_Single_Sudakov()
{
  if (p_bp) delete p_bp;
}
