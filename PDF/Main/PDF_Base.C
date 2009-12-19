#include "PDF/Main/PDF_Base.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PDF::PDF_Base
#define PARAMETER_TYPE PDF::PDF_Arguments
#include "ATOOLS/Org/Getter_Function.C"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Info_Key.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace PDF;
using namespace ATOOLS;

PDF_Base::PDF_Base():
  m_type("none"), m_orderas(-1), m_member(0), m_exponent(1.),
  m_rescale(1.), m_asmz(0.0), m_fac_scale_factor(1.) {

  if (rpa.gen.Variable("FACTORIZATION_SCALE_FACTOR")!="")
    m_fac_scale_factor = ToType<double>(rpa.gen.Variable("FACTORIZATION_SCALE_FACTOR"));
  if (m_fac_scale_factor!=1.0) 
    msg_Debugging()<<METHOD<<"(): Setting scale factor "<<m_fac_scale_factor<<"\n";
}

PDF_Base::~PDF_Base()
{
}

bool PDF_Base::Collinear(const double kp2) const
{
  return true;
}

PDF_Base *PDF_Base::GetBasicPDF() 
{
  return this;
}

double PDF_Base::Cut(const std::string &type)
{
  return ATOOLS::UNDEFINED_LOWER;
}

void PDF_Base::SetBounds()
{
  m_rq2min=m_q2min/m_fac_scale_factor;
  m_rq2max=m_q2max/m_fac_scale_factor;
}

void PDF_Base::SetPDFMember()
{
}

void PDF_Base::Calculate(double x,double Q2) 
{
  Q2*=m_fac_scale_factor;
  if(Q2<m_q2min) {
    msg_Error()<<METHOD<<"(): Q-range violation Q = "<<sqrt(Q2)
	       <<" < "<<sqrt(m_q2min)<<". Set Q -> "
	       <<sqrt(m_q2min)<<"."<<std::endl;
    Q2=1.000001*m_q2min;
  }
  if(Q2>m_q2max) {
    msg_Error()<<METHOD<<"(): Q-range violation Q = "<<sqrt(Q2)
	       <<" > "<<sqrt(m_q2max)<<". Set Q -> "
	       <<sqrt(m_q2max)<<"."<<std::endl;
    Q2=0.999999*m_q2max;
  }
  if(x<m_xmin*m_rescale) {
    msg_Error()<<METHOD<<"(): x = "<<x<<" ("<<m_rescale
	       <<") < "<<m_xmin<<". Set x -> "
	       <<m_xmin<<"."<<std::endl;
    x=1.000001*m_xmin*m_rescale;
  }
  if(x>m_xmax*m_rescale) {
    msg_Error()<<METHOD<<"(): x = "<<x<<" ("<<m_rescale
	       <<") > "<<m_xmax<<". Set x -> "
	       <<m_xmax<<"."<<std::endl;
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

void PDF_Base::AssignKeys(ATOOLS::Integration_Info *const info)
{
}

void PDF_Base::ShowSyntax(const size_t i)
{
  if (!msg_LevelIsInfo() || i==0) return;
  msg_Out()<<METHOD<<"(): {\n\n"
	   <<"   // available PDF sets (specified by PDF_SET=<value>)\n\n";
  PDF_Getter_Function::PrintGetterInfo(msg->Out(),25);
  msg_Out()<<"\n}"<<std::endl;
}

