//bof
//Version: 2 ADICIC++-0.0/2004/08/06

//Implementation of Sudakov_Calculator.H.



#include "Running_AlphaS.H"////////////////////////////////////////////////////
#include "Sudakov_Calculator.H"





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





#include "Sudakov_Calculator.tpt.cc"





//=============================================================================



template struct Sudakov_Info<Dipole::qqbar>;
template struct Sudakov_Info<Dipole::qg>;
template struct Sudakov_Info<Dipole::gqbar>;
template struct Sudakov_Info<Dipole::gg>;



template<Dipole::Type D>
const Dipole::Type Sudakov_Info<D>::Dipoletype=Dipole::qqbar;
template<Dipole::Type D> const short Sudakov_Info<D>::X1power=2;
template<Dipole::Type D> const short Sudakov_Info<D>::X3power=2;
template<Dipole::Type D>
const double Sudakov_Info<D>::Colourfactor=1.5/*0.75*/*M_PI;

template<>
const Dipole::Type Sudakov_Info<Dipole::qg>::Dipoletype=Dipole::qg;
template<> const short Sudakov_Info<Dipole::qg>::X1power=2;
template<> const short Sudakov_Info<Dipole::qg>::X3power=3;
template<>
const double Sudakov_Info<Dipole::qg>::Colourfactor=4.0/*2.0*/*M_PI/3.0;

template<>
const Dipole::Type Sudakov_Info<Dipole::gqbar>::Dipoletype=Dipole::gqbar;
template<> const short Sudakov_Info<Dipole::gqbar>::X1power=3;
template<> const short Sudakov_Info<Dipole::gqbar>::X3power=2;
template<>
const double Sudakov_Info<Dipole::gqbar>::Colourfactor=4.0/*2.0*/*M_PI/3.0;

template<>
const Dipole::Type Sudakov_Info<Dipole::gg>::Dipoletype=Dipole::gg;
template<> const short Sudakov_Info<Dipole::gg>::X1power=3;
template<> const short Sudakov_Info<Dipole::gg>::X3power=3;
template<>
const double Sudakov_Info<Dipole::gg>::Colourfactor=4.0/*2.0*/*M_PI/3.0;



//=============================================================================



//So far there is no static Sudakov_Calculator.
int Sudakov_Calculator::s_count=0;
const int& Sudakov_Calculator::InStore=Sudakov_Calculator::s_count;

bool   Sudakov_Calculator::s_isalphasrun=true;    //AlphaS treatment flag
double Sudakov_Calculator::s_k2tmin=1.0;          //GeV^2
double Sudakov_Calculator::s_k2tmax=8100.0;       //GeV^2
double Sudakov_Calculator::s_alphasfix=0.12;      //coupling

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
