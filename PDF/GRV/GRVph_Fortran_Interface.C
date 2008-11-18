#include "GRVph_Fortran_Interface.H"
#include "Run_Parameter.H"
#include "Model_Base.H"
#include "Message.H"
#include <unistd.h> 

using namespace PDF;
using namespace ATOOLS;


extern "C" {
  void  grvglo_(float &,float &,float&,float&,float&,float&,float&,float&);
}

GRVph_Fortran_Interface::GRVph_Fortran_Interface(const ATOOLS::Flavour _bunch) 
{
  m_xmin=1.e-5;
  m_xmax=1.;
  m_q2min=.25;
  m_q2max=1.e6;

  m_bunch = _bunch;
  m_d = m_u = m_s = m_c = m_b = m_g = 0.;
  
  for (int i=1;i<6;i++) {
    m_partons.push_back(Flavour((kf_code)(i)));
    m_partons.push_back(Flavour((kf_code)(i)).Bar());
  }
  m_partons.push_back(Flavour(kf_gluon));
  m_partons.push_back(Flavour(kf_jet));
  m_partons.push_back(Flavour(kf_quark));
  m_partons.push_back(Flavour(kf_quark).Bar());                               
}

PDF_Base * GRVph_Fortran_Interface::GetCopy()
{
  return new GRVph_Fortran_Interface(m_bunch);
}

void GRVph_Fortran_Interface::Output() {}

void GRVph_Fortran_Interface::Calculate(double _x,double _Q2) 
{
  float x = _x, Q2 = _Q2*m_fac_scale_factor;
  
  grvglo_(x,Q2,m_u,m_d,m_s,m_c,m_b,m_g);
}

double GRVph_Fortran_Interface::GetXPDF(const ATOOLS::Flavour infl) 
{
  double value = 0.;

  if (infl == Flavour(kf_gluon)) value = m_g;
  if (infl.Kfcode()==1)           value = m_d;
  if (infl.Kfcode()==2)           value = m_u;
  if (infl.Kfcode()==3)           value = m_s;
  if (infl.Kfcode()==4)           value = m_c;
  if (infl.Kfcode()==5)           value = m_b;
  
  value  *= MODEL::s_model->ScalarFunction(std::string("alpha_QED"),sqr(rpa.gen.Ecms()));
  
  return value;
}

void GRVph_Fortran_Interface::AssignKeys(ATOOLS::Integration_Info *const info)
{
}
