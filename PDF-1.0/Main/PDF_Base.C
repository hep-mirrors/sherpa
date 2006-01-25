#include "PDF_Base.H"

#include "Message.H"
#include "Info_Key.H"
#include "Run_Parameter.H"

using namespace PDF;

PDF_Base::PDF_Base():
  m_type("none"),
  m_exponent(1.), m_rescale(1.),
  m_ren_scale_factor(1.)
{
  m_ren_scale_factor = ATOOLS::rpa.gen.RenormalizationScaleFactor();
  if(m_ren_scale_factor<=0.0) {
    m_ren_scale_factor=1.0;
    ATOOLS::msg.Error()<<"{ "<<__PRETTY_FUNCTION__<<": Confused about the "
		       <<"!Renormalization! scale factor setting: "
		       <<"Set it to 1.0 and go on.}\n";
  }
}

PDF_Base::~PDF_Base() {
  for(size_t i=0; i<m_copies.size(); ++i) {
    if(m_copies[i]) delete m_copies[i];
  }
  //std::cout<<"          PDF_BASE DESTRUCT "<<m_copies.size()<<std::endl;///////
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

void PDF_Base::Calculate(double x,double Q2)
{
  Calculate(x,0.,0.,Q2*m_ren_scale_factor);
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
