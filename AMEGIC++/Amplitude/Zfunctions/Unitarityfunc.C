#include "AMEGIC++/Amplitude/Zfunctions/Basic_Func.H"
#include "AMEGIC++/Amplitude/Zfunctions/Basic_Sfuncs.H"
#include "AMEGIC++/String/String_Generator.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/MathTools.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Unitarityfunc::Unitarityfunc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS)
  : Basic_Func(_sgen,_BS)  
{
  m_n = m_lambda2 = 0.;
  if (!(MODEL::s_model->Name()=="SM+AGC")) return;
  m_n        = MODEL::s_model->ScalarConstant(std::string("UNITARIZATION_N"));
  m_lambda2  = sqr(MODEL::s_model->ScalarConstant(std::string("UNITARIZATION_SCALE")));
} 

Kabbala Unitarityfunc::U()
{ 
  Complex uf=Ucalc();
  return sgen->GetSFnumber(uf,11);
}

Complex Unitarityfunc::Ucalc()
{ 
  if (!m_n>0.||!m_lambda2>0.) return Complex(1.,0.);
  Vec4D h = BS->Momentum(0);
  if (BS->Sign(1)==BS->Sign(0)) h+= BS->Momentum(1);
  return Complex(pow(1.+h.Abs2()/m_lambda2,-m_n),0.);
}




