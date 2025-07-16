#include "AMEGIC++/Amplitude/Zfunctions/Basic_Func.H"
#include "AMEGIC++/Amplitude/Zfunctions/Basic_Sfuncs.H"
#include "AMEGIC++/String/String_Generator.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Unitarityfunc::Unitarityfunc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS)
  : Basic_Func(_sgen,_BS)  
{
  m_n = m_lambda2 = 0.;
  if (!(MODEL::s_model->Name()=="SM+AGC")) return;
  m_n         = MODEL::s_model->ScalarConstant(std::string("UNITARIZATION_N"));
  m_m         = MODEL::s_model->ScalarConstant(std::string("UNITARIZATION_M"));
  m_lambda2   = sqr(MODEL::s_model->ScalarConstant(std::string("UNITARIZATION_SCALE")));
  m_n3        = MODEL::s_model->ScalarConstant(std::string("UNITARIZATION_N3"));
  m_m3        = MODEL::s_model->ScalarConstant(std::string("UNITARIZATION_M3"));
  m_lambda2_3 = sqr(MODEL::s_model->ScalarConstant(std::string("UNITARIZATION_SCALE3")));
  m_n4        = MODEL::s_model->ScalarConstant(std::string("UNITARIZATION_N4"));
  m_m4        = MODEL::s_model->ScalarConstant(std::string("UNITARIZATION_M4"));
  m_lambda2_4 = sqr(MODEL::s_model->ScalarConstant(std::string("UNITARIZATION_SCALE4")));
}

Kabbala Unitarityfunc::U(const int & n)
{ 
  Complex uf=Ucalc(n);
  return (n==3?sgen->GetSFnumber(uf,11):sgen->GetSFnumber(uf,12));
}

Complex Unitarityfunc::Ucalc(const int & n)
{
  DEBUG_FUNC((n==3?"aTGC":(n==4?"aQGC":"undefined")));
  double nn(m_n), mm(m_m), lambda2(m_lambda2);
  if (n==3) { nn = m_n3; mm = m_m3; lambda2 = m_lambda2_3; }
  if (n==4) { nn = m_n4; mm = m_m4; lambda2 = m_lambda2_4; }
  if (!(nn>0.)||!(lambda2>0.)) return Complex(1.,0.);
  Vec4D h = BS->Momentum(0);
  if (BS->Sign(1)==BS->Sign(0)) h+= BS->Momentum(1);
  double ff(pow(1.+pow(h.Abs2()/lambda2,mm),-nn));
  msg_Debugging()<<"n = "<<nn<<" ,  m = "<<mm
                 <<" ,  \\Lambda^2 = "<<lambda2
                 <<" ,  shat = "<<h.Abs2()
                 <<" => f(shat) = "<<ff<<std::endl;
  return Complex(ff,0.);
}




