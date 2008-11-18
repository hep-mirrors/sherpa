#include "PDF_Base.H"

#include "Message.H"
#include "Info_Key.H"
#include "Run_Parameter.H"
#include "Exception.H"
#include "MyStrStream.H"

using namespace PDF;
using namespace ATOOLS;



PDF_Base::Box PDF_Base::s_box=PDF_Base::Box();

PDF_Base::Box::~Box() {
  if (exh->LastSignal()!=0 || exh->LastException()!=NULL) return;
  for(unsigned i=0; i<v_pdfp.size(); ++i) if(v_pdfp[i]) delete v_pdfp[i];
}

std::size_t PDF_Base::Box::TrueEntryNumber() const {
  std::size_t n=0;
  for(unsigned i=0; i<v_pdfp.size(); ++i) if(v_pdfp[i]) ++n;
  return n;
}



PDF_Base::PDF_Base()
  : m_type("none"), m_exponent(1.), m_rescale(1.), m_fac_scale_factor(1.) {

  s_box.v_pdfp.push_back(this);
  if (rpa.gen.Variable("FACTORIZATION_SCALE_FACTOR")!="")
    m_fac_scale_factor = ToType<double>(rpa.gen.Variable("FACTORIZATION_SCALE_FACTOR"));
  msg_Tracking()<<s_box.v_pdfp.size()<<"|"<<s_box.TrueEntryNumber()
                <<"    PDF_Base CONSTRUCT "<<m_copies.size()<<" "<<this
                <<std::endl;
  if (m_fac_scale_factor!=1.0) 
    msg_Debugging()<<METHOD<<"(): Setting scale factor "<<m_fac_scale_factor<<"\n";
}

PDF_Base::~PDF_Base() {
  for(std::size_t i=0; i<s_box.v_pdfp.size(); ++i)
    if(this==s_box.v_pdfp[i]) s_box.v_pdfp[i]=NULL;
  msg_Tracking()<<s_box.v_pdfp.size()<<"|"<<s_box.TrueEntryNumber()
		<<"    PDF_Base DESTRUCT "<<m_copies.size()<<" "<<this
		<<std::endl;
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
