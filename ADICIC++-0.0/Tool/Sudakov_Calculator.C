//bof
//Version: 2 ADICIC++-0.0/2004/08/10

//Implementation of Sudakov_Calculator.H.



#include "Running_AlphaS.H"////////////////////////////////////////////////////
#include "Dipole_Parameter.H"
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



const bool Sudakov_Calculator::sf_start=Dipole_Parameter::ForceFirstInit();

//So far there is no static Sudakov_Calculator.
int Sudakov_Calculator::s_count=0;
const int& Sudakov_Calculator::InStore=Sudakov_Calculator::s_count;

//AlphaS treatment flag and fixed coupling
bool   Sudakov_Calculator::s_isalphasrun=Dipole_Parameter::IsAlphaSRunning();
double Sudakov_Calculator::s_alphasfix=Dipole_Parameter::AlphaSFix();
double Sudakov_Calculator::s_k2tmin=Dipole_Parameter::MinOfK2t();    //GeV^2
double Sudakov_Calculator::s_k2tmax=Dipole_Parameter::MaxOfK2t();    //GeV^2

double Sudakov_Calculator::s_approx=Sudakov_Calculator::s_alphasfix;

Function_Base* Sudakov_Calculator::s_pas=NULL;

Sudakov_Calculator::AlphaSCorr_Func
Sudakov_Calculator::GetAlphaSCorr=&Sudakov_Calculator::FixAlphaSCorr;



//=============================================================================



void Sudakov_Calculator::ShowParameters() {    //Static.
  cout<<endl;
  cout<<"==========================================="<<endl;
  cout<<"    Valid Sudakov_Calculator parameters"<<endl;
  cout<<"-------------------------------------------"<<endl;
  cout<<"Running AlphaS is allowed for....."<<s_isalphasrun<<".\n";
  cout<<"Running AlphaS is initialized....."<<(s_isalphasrun && s_pas)<<".\n";
  cout<<"Fixed AlphaS is set to............"<<s_alphasfix<<"."<<endl;
  cout<<"AlphaS approximation is set to...."<<s_approx<<"."<<endl;
  cout<<"Dipole shower cut-off scale is set to...."<<s_k2tmin<<" GeV^2."<<endl;
  cout<<"Dipole shower maximum scale is set to...."<<s_k2tmax<<" GeV^2."<<endl;
  cout<<"==========================================="<<endl;
}





const bool Sudakov_Calculator::AdjustParameters() {    //Static.
  s_alphasfix=Dipole_Parameter::AlphaSFix();
  s_approx=s_alphasfix;
  s_k2tmin=Dipole_Parameter::MinOfK2t();
  s_k2tmax=Dipole_Parameter::MaxOfK2t();
  if(s_isalphasrun != Dipole_Parameter::IsAlphaSRunning()) {
    if(s_isalphasrun) {
      s_pas=NULL;
      GetAlphaSCorr=&FixAlphaSCorr;
    }
    s_isalphasrun=Dipole_Parameter::IsAlphaSRunning();
  }
  return true;
}





const bool Sudakov_Calculator::Init(MODEL::Model_Base* pmod) {    //Static.
  if(s_isalphasrun==false) return false;
  if(pmod) {
    //The Running_AlphaS object physically resides in the initialized Model.
    //So the following is simply an assignment (no new-operator is needed).
    s_pas=pmod->GetScalarFunction("alpha_S");
  }
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
