#include "Single_XS.H"
#include "ISR_Handler.H"
#include "Phase_Space_Handler.H"
#include "Run_Parameter.H"
#include "Message.H"

using namespace EXTRAXS;

Single_XS::Single_XS(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours,
		     const int scalescheme,const int kfactorscheme,const double scalefactor,
		     BEAM::Beam_Spectra_Handler *const beamhandler,
		     PDF::ISR_Handler *const isrhandler,
		     ATOOLS::Selector_Data *const selectordata):
  XS_Base(nin,nout,flavours,scalescheme,kfactorscheme,scalefactor,
	  beamhandler,isrhandler,selectordata)
{
  p_selected=this;
}

Single_XS::Single_XS(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours):
  XS_Base(nin,nout,flavours)
{
  p_selected=this;
}

void Single_XS::WriteOutXSecs(std::ofstream &outfile)
{
  outfile.precision(12);
  outfile<<m_name<<"  "<<m_totalxs<<"  "<<m_max<<"  "<<m_totalerr<<" "
	 <<m_totalsum<<" "<<m_totalsumsqr<<" "<<m_n<<" "
	 <<m_ssum<<" "<<m_ssumsqr<<" "<<m_ssigma2<<" "<<m_sn<<std::endl; 
}

bool Single_XS::CalculateTotalXSec(const std::string &resultpath) 
{ 
  m_n=0;
  m_last=m_lastlumi=m_lastdxs=0.0;
  m_totalxs=m_totalsum=m_totalsumsqr=m_totalerr=0.0;
  if (p_pshandler) {
    m_totalxs=p_pshandler->Integrate()/ATOOLS::rpa.Picobarn();
    if (!(ATOOLS::IsZero((m_n*m_totalxs-m_totalsum)/(m_n*m_totalxs+m_totalsum)))) {
      ATOOLS::msg.Error()<<"Single_XS::CalculateTotalXSec(..): ("<<this<<")"<<std::endl
			 <<ATOOLS::om::red<<"   Integrator result and internal summation "
			 <<"do not coincide!"<<ATOOLS::om::reset<<std::endl
			 <<"  "<<m_name<<" : "<<m_totalxs<<" vs. "<<m_totalsum/m_n<<std::endl;
    }
    if (m_totalxs>0.) return true;
    return false;
  }
  ATOOLS::msg.Error()<<"Single_XS::CalculateTotalXSec(): ("<<this<<")"<<ATOOLS::om::red
		     <<"No pointer to Phase_Space_Handler ! Abort."<<ATOOLS::om::reset<<std::endl;
  return false;
}

bool Single_XS::OneEvent() 
{ 
  return p_activepshandler->OneEvent(); 
}

ATOOLS::Blob_Data_Base *Single_XS::WeightedEvent(const int mode) 
{ 
  return p_activepshandler->WeightedEvent(mode); 
}

void Single_XS::SetTotal() 
{
  m_totalxs=TotalResult();//m_totalsum/m_n; 
  m_totalerr=TotalVar();//sqrt((m_n*m_totalsumsqr-ATOOLS::sqr(m_totalsum))/(m_n-1))/m_n;
  msg_Info()<<"      xs for "<<ATOOLS::om::bold<<m_name<<" : "
		    <<ATOOLS::om::blue<<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb"
		    <<ATOOLS::om::reset<<" +/- ( "<<ATOOLS::om::red<<m_totalerr*ATOOLS::rpa.Picobarn()<<" pb = "
		    <<m_totalerr/m_totalxs*100.<<" %"<<ATOOLS::om::reset<<" )"<<std::endl
		    <<"       max : "<<m_max<<std::endl;
}

void Single_XS::AddPoint(const double value) 
{
  m_n++;
  m_sn++;
  m_ssum    += value;
  m_ssumsqr += value*value;
  if (value>m_max)  m_max  = value;
  if (value>m_smax) m_smax = value;
}

void Single_XS::OptimizeResult()
{
  double ssigma2 = (m_ssumsqr/m_sn - ATOOLS::sqr(m_ssum/m_sn))/(m_sn-1);
  m_ssigma2  += 1./ssigma2; 
  m_totalsum += m_ssum/ssigma2/m_sn;
  m_totalsumsqr+= m_ssumsqr/ssigma2/m_sn;
  m_ssum     = 0.;
  m_ssumsqr  = 0.;
  m_sn       = 0;
}

void Single_XS::ResetMax(int flag)
{
  if (flag==0) {
    if (m_vsmax.size()>1) {
      m_vsmax.erase(m_vsmax.begin());
      m_vsn.erase(m_vsn.begin());
    }
    m_vsmax.back() = ATOOLS::Max(m_smax,m_vsmax.back());
    m_vsn.back()   = m_n;
  }
  else {
    if (flag==2 && m_vsmax.size()==4) {
      m_vsmax.erase(m_vsmax.begin());
      m_vsn.erase(m_vsn.begin());
    }
    m_vsmax.push_back(m_smax);
    m_vsn.push_back(m_n);
    if (flag==2) m_smax = 0.;
  }
  m_max  = 0.;
//   cout<<"Single_Process::ResetMax "<<flag<<": "<<endl;
//   for (size_t i=0;i<m_vsmax.size();i++) cout<<m_vsmax[i]<<" "<<m_vsn[i]<<endl;
  for (size_t i=0;i<m_vsmax.size();i++) m_max=ATOOLS::Max(m_max,m_vsmax[i]);
}

double Single_XS::Differential(const double s,const double t,const double u)
{
  m_lastdxs=(*p_regulator)(operator()(s,t,u));
  if (m_lastdxs<=0.) return m_lastdxs=m_last=0.;
  if (p_isrhandler && m_nin==2) { 
    if (p_isrhandler->On()) m_lastlumi=p_isrhandler->Weight(p_flavours); 
    else m_lastlumi=1.;
  }
  else m_lastlumi=1.;
  return m_last=m_lastdxs*m_lastlumi;
}

double Single_XS::Differential2() 
{
  if (p_isrhandler && m_nin==2) {
    if (p_flavours[0]==p_flavours[1] || p_isrhandler->On()==0) return 0.;
    double tmp=m_lastdxs*p_isrhandler->Weight2(p_flavours); 
    m_last+=tmp;
    return tmp;
  }
  return 0;
}

double Single_XS::operator()(const double s,const double t,const double u) 
{
  ATOOLS::msg.Error()<<"Single_XS::operator()("<<s<<","<<t<<","<<u<<"): "
		     <<"Virtual method called!"<<std::endl;
  return 0.;
}

void Single_XS::SetISR(PDF::ISR_Handler *const isrhandler) 
{ 
  p_isrhandler=isrhandler; 
}

XS_Base *const Single_XS::operator[](const size_t i) const 
{
  const XS_Base *value=this;
  return (XS_Base*)value;
}

size_t Single_XS::Size() const
{ 
  return 1; 
}

