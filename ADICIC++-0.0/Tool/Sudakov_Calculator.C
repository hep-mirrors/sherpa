//bof
//Version: 1 ADICIC++-0.0/2004/06/08

//Implementation of Sudakov_Calculator.H.



#include "Sudakov_Calculator.H"
#include "Sudakov_Calculator.dat.cc"

#include "Running_AlphaS.H"////////////////////////////////////////////////////





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





#include "Sudakov_Calculator.tpt.cc"





//=============================================================================



//So far there is no static Sudakov_Calculator.
int Sudakov_Calculator::s_count=0;
const int& Sudakov_Calculator::InStore=Sudakov_Calculator::s_count;

double Sudakov_Calculator::s_approx=Sudakov_Calculator::s_alphasfix;

Function_Base* Sudakov_Calculator::s_pas=NULL;

Sudakov_Calculator::AlphaSCorr_Func
Sudakov_Calculator::GetAlphaSCorr=&Sudakov_Calculator::FixAlphaSCorr;



//=============================================================================



const bool Sudakov_Calculator::Init(MODEL::Model_Base* pmod) {
  //Static method.
  if(s_isalphasrun==false) return false;
  if(pmod) s_pas=pmod->GetScalarFunction("alpha_S");
  else {
    //static MODEL::Running_AlphaS as(0.1188,8315.25,1);
    static MODEL::Running_AlphaS as(0.118,8315.0,1);
    s_pas=&as;
  }
  assert(s_pas);
  s_approx=(*s_pas)(s_k2tmin);
  GetAlphaSCorr=&RunAlphaSCorr;
  return true;
}



//=============================================================================



template class Sudakov<Dipole::qqbar>;
template class Sudakov<Dipole::qg>;
template class Sudakov<Dipole::gqbar>;
template class Sudakov<Dipole::gg>;



//=============================================================================





//eof
