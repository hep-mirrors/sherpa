#include "PDF_Base.H"

#include "Info_Key.H"
#include "Run_Parameter.H"

using namespace PDF;

PDF_Base::PDF_Base():
  m_type("none"),
  m_exponent(1.),
  m_rescale(1.),
  m_renormalization_scale_factor(1.)
{
  m_renormalization_scale_factor = ATOOLS::rpa.gen.RenormalizationScaleFactor();
  //  std::cout<<" mu_R-fac in PDF:"<<m_renormalization_scale_factor;
}

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

void PDF_Base::Calculate(double x,double Q2)
{
  //  std::cout<<"PDF_Base::Calculate("<<Q2<<") Called"<<std::endl;
  Calculate(x,0.,0.,Q2*m_renormalization_scale_factor);
}

void PDF_Base::Extract(const ATOOLS::Flavour flavour,const double x) 
{ 
  PDF_Base *pdf=GetBasicPDF();
  pdf->m_rescale-=x;
  pdf->m_extracted.push_back(std::pair<ATOOLS::Flavour,double>(flavour,x));
  for (std::vector<PDF_Base*>::iterator cit=pdf->m_copies.begin();
       cit!=pdf->m_copies.end();++cit) {
    (*cit)->Extract(flavour,x);
  }
}

void PDF_Base::Reset()
{ 
  PDF_Base *pdf=GetBasicPDF();
  pdf->m_rescale=1.;
  pdf->m_extracted.clear();
  for (std::vector<PDF_Base*>::iterator cit=pdf->m_copies.begin();
       cit!=pdf->m_copies.end();++cit) {
    (*cit)->Reset();
  }
}
