#include "PDF/Main/PDF_Base.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PDF::PDF_Base
#define PARAMETER_TYPE PDF::PDF_Arguments
#include "ATOOLS/Org/Getter_Function.C"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"

using namespace PDF;
using namespace ATOOLS;

PDF_Base::PDF_Base():
  m_type("none"), m_member(0), m_exponent(1.),
  m_rescale(1.) {

  Data_Reader dr(" ",";","!","=");
  m_lhef_number = dr.GetValue<int>("LHEF_PDF_NUMBER",-1);
}

PDF_Base::~PDF_Base()
{
}

double PDF_Base::AlphaSPDF(const double &q2)
{
  return -1.0;
}

PDF_Base *PDF_Base::GetBasicPDF() 
{
  return this;
}

void PDF_Base::SetBounds()
{
  m_rq2min=m_q2min;
  m_rq2max=m_q2max;
}

void PDF_Base::SetPDFMember()
{
}

void PDF_Base::Calculate(double x,double Q2) 
{
  if(Q2<m_q2min) {
    static double lasterr(-1.0);
    if (Q2!=lasterr)
    msg_Error()<<METHOD<<"(): Q-range violation Q = "<<sqrt(Q2)
	       <<" < "<<sqrt(m_q2min)<<". Set Q -> "
	       <<sqrt(m_q2min)<<"."<<std::endl;
    lasterr=Q2;
    Q2=1.000001*m_q2min;
  }
  if(Q2>m_q2max) {
    static double lasterr(-1.0);
    if (Q2!=lasterr)
    msg_Error()<<METHOD<<"(): Q-range violation Q = "<<sqrt(Q2)
	       <<" > "<<sqrt(m_q2max)<<". Set Q -> "
	       <<sqrt(m_q2max)<<"."<<std::endl;
    lasterr=Q2;
    Q2=0.999999*m_q2max;
  }
  if(x<m_xmin*m_rescale) {
    static double lasterr(-1.0);
    if (x!=lasterr)
    msg_Error()<<METHOD<<"(): x = "<<x<<" ("<<m_rescale
	       <<") < "<<m_xmin<<". Set x -> "
	       <<m_xmin<<"."<<std::endl;
    lasterr=x;
    x=1.000001*m_xmin*m_rescale;
  }
  if(x>m_xmax*m_rescale) {
    static double lasterr(-1.0);
    if (x!=lasterr)
    msg_Error()<<METHOD<<"(): x = "<<x<<" ("<<m_rescale
	       <<") > "<<m_xmax<<". Set x -> "
	       <<m_xmax<<"."<<std::endl;
    lasterr=x;
    x=0.999999*m_xmax*m_rescale;
  }
  return CalculateSpec(x,Q2);
}

void PDF_Base::SingleExtract(const ATOOLS::Flavour flavour,const double x) 
{
  m_rescale-=x;
  m_extracted.push_back(std::pair<ATOOLS::Flavour,double>(flavour,x));
  for (std::vector<PDF_Base*>::iterator cit(m_copies.begin());
       cit!=m_copies.end();++cit) {
    (*cit)->SingleExtract(flavour,x);
  }
}

void PDF_Base::SingleReset()
{ 
  m_rescale=1.;
  m_extracted.clear();
  for (std::vector<PDF_Base*>::iterator cit(m_copies.begin());
       cit!=m_copies.end();++cit) {
    (*cit)->SingleReset();
  }
}

void PDF_Base::ShowSyntax(const size_t i)
{
  if (!msg_LevelIsInfo() || i==0) return;
  msg_Out()<<METHOD<<"(): {\n\n"
	   <<"   // available PDF sets (specified by PDF_SET=<value>)\n\n";
  PDF_Getter_Function::PrintGetterInfo(msg->Out(),25);
  msg_Out()<<"\n}"<<std::endl;
}
