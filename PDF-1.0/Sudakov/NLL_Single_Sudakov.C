#include "NLL_Single_Sudakov.H"

#include "Run_Parameter.H"
#include "Message.H"
#include "Shell_Tools.H"

using namespace SHERPA;
using namespace ATOOLS;

NLL_Single_Sudakov::NLL_Single_Sudakov(NLL_Branching_Probability_Base * bp,int mode):
  m_qmax(0.),
  m_qmin(0.),
  m_calcmode((Sudakov::code)mode),
  m_cutmode((Sudakov::code)mode),
  p_bp(bp),
  m_gauss(bp) {}

bool NLL_Single_Sudakov::Initialize(double _m_qmin,double _m_qmax) 
{
  m_qmin=_m_qmin;
  m_qmax=_m_qmax;

  m_calcmode = (Sudakov::code)(m_calcmode&896);
  m_cutmode  = (Sudakov::code)(m_cutmode&7);

  if (m_qmin==m_qmax) {
    m_qmin=sqrt(rpa.gen.Ycut())*rpa.gen.Ecms();
    m_qmax=2.*rpa.gen.Ecms();
  }

  m_qlimit = rpa.gen.Ecms() * sqrt(1./3.); 
  if (m_calcmode==Sudakov::table) {
    bool create=true;
    if (p_bp->Name().length()>0) {
      if (!m_log_delta.ReadIn((m_outpath+std::string("/")+p_bp->Name()).c_str())) {
	ATOOLS::MakeDir(m_outpath.c_str(),0755);
      }
      else {
	create=false;
      }
    }
    // initialize table
    if (create) {
      m_calcmode=Sudakov::create_table;
      m_log_delta.Init(*this,m_qmin,m_qmax,600);
      m_calcmode=Sudakov::table;
      if (p_bp->Name().length()>0) 
	m_log_delta.WriteOut((m_outpath+std::string("/")+p_bp->Name()).c_str());
    }
  }
  return true;
}

double NLL_Single_Sudakov::Log(double Q, double q) 
{
  if (m_calcmode == Sudakov::analytic) {
    double sum=p_bp->IntGamma(q,Q);
    if (sum!=-1.) {
      return sum;
    }
    else {
      ATOOLS::msg.Error()<<"NLL_Single_Sudakov::Log(..): "
			 <<"Integrated branching probability not analytically calculated yet! "<<std::endl
			 <<"      try numeric integration "<<std::endl;
    }
  }
  else if (m_calcmode == Sudakov::table) {
    if (!IsEqual(q,m_qmin)) {
      // -------------- use numeric routine -------------
      m_calcmode=Sudakov::numeric;
      double logdelta=Log(Q,q);
      m_calcmode=Sudakov::table;
      return logdelta;

      // ---------------- default: recreate table: -------------
      msg.Error()<<"ERROR in NLL_Single_Sudakov : "<<std::endl
		 <<"   Table calculated with qmin="<<m_qmin<<" but called with "<<q<<std::endl;
      m_qmin=q;
      m_calcmode=Sudakov::create_table;
      m_log_delta.Init(*this,m_qmin,m_qmax,600);
      m_calcmode=Sudakov::table;
    }
    return m_log_delta(Q);
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
