#include "Single_XS.H"

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
  outfile<<m_name<<"  "<<m_totalxs<<"  "<<m_max<<"  "<<m_totalerr<<" "
	 <<m_totalsum<<" "<<m_totalsumsqr<<" "<<m_n<<std::endl;
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

ATOOLS::Blob_Data_Base *Single_XS::WeightedEvent() 
{ 
  return p_activepshandler->WeightedEvent(); 
}

void Single_XS::SetTotal() 
{
  m_totalxs=m_totalsum/m_n; 
  m_totalerr=sqrt((m_n*m_totalsumsqr-ATOOLS::sqr(m_totalsum))/(m_n-1))/m_n;
  msg_Info()<<"      xs for "<<ATOOLS::om::bold<<m_name<<" : "
		    <<ATOOLS::om::blue<<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb"
		    <<ATOOLS::om::reset<<" +/- ( "<<ATOOLS::om::red<<m_totalerr<<" pb = "
		    <<m_totalerr/m_totalxs*100.<<" %"<<ATOOLS::om::reset<<" )"<<std::endl
		    <<"       max : "<<m_max<<std::endl;
}

void Single_XS::AddPoint(const double value) 
{
  m_n++;
  m_totalsum+=value;
  m_totalsumsqr+=value*value;
  if (value>m_max) m_max=value;
}

double Single_XS::Differential(const double s,const double t,const double u)
{
  m_lastdxs=operator()(s,t,u);
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

size_t Single_XS::Size() 
{ 
  return 1; 
}

void Single_XS::SetISR(PDF::ISR_Handler *const isrhandler) 
{ 
  p_isrhandler=isrhandler; 
}

