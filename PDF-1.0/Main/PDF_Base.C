#include "PDF_Base.H"

#include "Info_Key.H"

using namespace PDF;

PDF_Base::PDF_Base():
  m_type("none"),
  m_exponent(1.),
  m_rescale(1.) {}

PDF_Base::~PDF_Base() {}

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

void PDF_Base::Extract(const ATOOLS::Flavour flavour,const double x) 
{ 
  m_rescale-=x;
  m_extracted.push_back(std::pair<ATOOLS::Flavour,double>(flavour,x));
  for (std::vector<PDF_Base*>::iterator cit=m_copies.begin();
       cit!=m_copies.end();++cit) {
    (*cit)->Extract(flavour,x);
  }
}

void PDF_Base::Reset()
{ 
  m_rescale=1.;
  m_extracted.clear();
  for (std::vector<PDF_Base*>::iterator cit=m_copies.begin();
       cit!=m_copies.end();++cit) {
    (*cit)->Reset();
  }
}
