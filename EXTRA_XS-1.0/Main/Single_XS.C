#include "Single_XS.H"
#include "Run_Parameter.H"
#include "ISR_Base.H"
#include "Running_AlphaS.H"

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>


using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

int fak(int N)
{
  if (N == 0) return 1;
  if (N < 0) return 0;  
  int res =1;
  for (int i=1;i<=N;i++) res *= i;
  return res;
}

Single_XS::Single_XS(int _nin,int _nout,Flavour * _fl,
		     PDF::ISR_Handler * _isr,BEAM::Beam_Spectra_Handler * _beam,
		     ATOOLS::Selector_Data * _seldata,
		     int _scalescheme,int _kfactorscheme,double _scalefactor) :
  XS_Base(_nin,_nout,_fl,_isr,_beam,_seldata,_scalescheme,_kfactorscheme,_scalefactor)
{
  double mass = 0.;
  for (int i=0;i<m_nin;i++)   mass += _fl[i].PSMass();
  m_thres     = ATOOLS::sqr(mass);
  mass        = 0.;
  for (int i=2;i<2+m_nout;i++) mass += _fl[i].PSMass();
  m_thres     = ATOOLS::Max(m_thres,sqr(mass));

  p_selected = this;
}

Single_XS::Single_XS(int _nin,int _nout,Flavour * _fl) :
  XS_Base(_nin,_nout,_fl)
{
  p_selected = this;
}


bool Single_XS::CalculateTotalXSec() 
{ 
  m_n=0;
  m_last=m_lastlumi=m_lastdxs=0.0;
  m_totalxs=m_totalsum=m_totalsumsqr=m_totalerr=0.0;
  if (p_ps) {
    m_totalxs = p_ps->Integrate()/ATOOLS::rpa.Picobarn();
    if (!(ATOOLS::IsZero((m_n*m_totalxs-m_totalsum)/(m_n*m_totalxs+m_totalsum)))) {
      msg.Error()<<"Result of PS-Integrator and internal summation to not coincide!"<<endl
		 <<"  "<<m_name<<" : "<<m_totalxs<<" vs. "<<m_totalsum/m_n<<endl;
    }
    if (m_totalxs>0.) return 1;
    return 0;
  }
  msg.Error()<<"Error in Single_Process::CalculateTotalXSec()."<<endl
	     <<"   No pointer to Phase_Space_Handler, implies abuse of this Single_XS."<<endl;
  return 0;
}

bool Single_XS::OneEvent()        { return (p_ps->OneEvent()); }

void Single_XS::SetTotalXS() {
  m_totalxs  = m_totalsum/m_n; 
  m_totalerr = sqrt( (m_n*m_totalsumsqr-ATOOLS::sqr(m_totalsum))/(m_n-1))/m_n;
  ATOOLS::msg.Events()<<"      xs for "<<m_name<<" : "
			 <<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb"
			 <<" +/- "<<m_totalerr/m_totalxs*100.<<"%,"<<endl
			 <<"       max : "<<m_max<<endl;
}

void Single_XS::AddPoint(const double value) {
  m_n++;
  m_totalsum             += value;
  m_totalsumsqr          += value*value;
  if (value>m_max) m_max  = value;
}

double Single_XS::Differential(double s,double t,double u)
{
  m_lastdxs = operator()(s,t,u);
  if (m_lastdxs <= 0.) return m_lastdxs = m_last = 0.;
  if ((p_isr) && m_nin==2) m_lastlumi = p_isr->Weight(p_fl);
                     else  m_lastlumi = 1.;

  return m_last = m_lastdxs * m_lastlumi;
}



double Single_XS::Differential2() {
  if ((p_isr) && m_nin==2) {
    if ((p_fl[0]==p_fl[1]) || (p_isr->On()==0) ) return 0.;
    double tmp = m_lastdxs * p_isr->Weight2(p_fl); 
    m_last    += tmp;
    return tmp;
  }
  return 0;
}


