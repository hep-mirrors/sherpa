//bof
//Version: 2 ADICIC++-0.0/2004/10/28

//Implementation of Sudakov_Calculator.H.



#include "Running_AlphaS.H"////////////////////////////////////////////////////
#include "Dipole_Parameter.H"
#include "Sudakov_Calculator.H"





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





#include "Sudakov_Calculator.tpt.cc"





//=============================================================================



template struct Sudakov_Info<Dipole::qqbar,Radiation::gluon>;
template struct Sudakov_Info<Dipole::qg,Radiation::gluon>;
template struct Sudakov_Info<Dipole::gqbar,Radiation::gluon>;
template struct Sudakov_Info<Dipole::gg,Radiation::gluon>;

template struct Sudakov_Info<Dipole::qg,Radiation::quark>;
template struct Sudakov_Info<Dipole::gqbar,Radiation::quark>;
template struct Sudakov_Info<Dipole::gg,Radiation::quark>;




template<Dipole::Type D> const Radiation::Group
Sudakov_Info<D,Radiation::gluon>::Radiationgroup=Radiation::gluon;
template<Dipole::Type D> const Radiation::Group
Sudakov_Info<D,Radiation::quark>::Radiationgroup=Radiation::quark;




template<> const Dipole::Type
Sudakov_Info<Dipole::qqbar,Radiation::gluon>::Dipoletype=Dipole::qqbar;
template<> const short
Sudakov_Info<Dipole::qqbar,Radiation::gluon>::X1power=2;
template<> const short
Sudakov_Info<Dipole::qqbar,Radiation::gluon>::X3power=2;
template<> const double
Sudakov_Info<Dipole::qqbar,Radiation::gluon>::Colourfactor=1.5*M_PI;


template<> const Dipole::Type
Sudakov_Info<Dipole::qg,Radiation::gluon>::Dipoletype=Dipole::qg;
template<> const short
Sudakov_Info<Dipole::qg,Radiation::gluon>::X1power=2;
template<> const short
Sudakov_Info<Dipole::qg,Radiation::gluon>::X3power=3;
template<> const double
Sudakov_Info<Dipole::qg,Radiation::gluon>::Colourfactor=4.0*M_PI/3.0;


template<> const Dipole::Type
Sudakov_Info<Dipole::gqbar,Radiation::gluon>::Dipoletype=Dipole::gqbar;
template<> const short
Sudakov_Info<Dipole::gqbar,Radiation::gluon>::X1power=3;
template<> const short
Sudakov_Info<Dipole::gqbar,Radiation::gluon>::X3power=2;
template<> const double
Sudakov_Info<Dipole::gqbar,Radiation::gluon>::Colourfactor=4.0*M_PI/3.0;


template<> const Dipole::Type
Sudakov_Info<Dipole::gg,Radiation::gluon>::Dipoletype=Dipole::gg;
template<> const short
Sudakov_Info<Dipole::gg,Radiation::gluon>::X1power=3;
template<> const short
Sudakov_Info<Dipole::gg,Radiation::gluon>::X3power=3;
template<> const double
Sudakov_Info<Dipole::gg,Radiation::gluon>::Colourfactor=4.0*M_PI/3.0;




template<> const Dipole::Type
Sudakov_Info<Dipole::qg,Radiation::quark>::Dipoletype=Dipole::qg;
template<> const double
Sudakov_Info<Dipole::qg,Radiation::quark>::Colourfactor=8.0*M_PI;


template<> const Dipole::Type
Sudakov_Info<Dipole::gqbar,Radiation::quark>::Dipoletype=Dipole::gqbar;
template<> const double
Sudakov_Info<Dipole::gqbar,Radiation::quark>::Colourfactor=8.0*M_PI;


template<> const Dipole::Type
Sudakov_Info<Dipole::gg,Radiation::quark>::Dipoletype=Dipole::gg;
template<> const double
Sudakov_Info<Dipole::gg,Radiation::quark>::Colourfactor=4.0*M_PI;



//=============================================================================



const bool Sudakov_Calculator::sf_start=Dipole_Parameter::ForceFirstInit();

//Mimic Ariadne.
const bool Sudakov_Calculator::sf_ariadne=false;//true;//false;
const bool& Sudakov_Calculator::Ariadne=Sudakov_Calculator::sf_ariadne;

//So far there is no static Sudakov_Calculator.
int Sudakov_Calculator::s_count=0;
const int& Sudakov_Calculator::InStore=Sudakov_Calculator::s_count;

//AlphaS treatment flag and fixed coupling.
bool   Sudakov_Calculator::s_isalphasrun=Dipole_Parameter::IsAlphaSRunning();
double Sudakov_Calculator::s_alphasfix=Dipole_Parameter::AlphaSFix();
double Sudakov_Calculator::s_k2tmin=Dipole_Parameter::MinOfK2t();    //GeV^2
double Sudakov_Calculator::s_k2tmax=Dipole_Parameter::MaxOfK2t();    //GeV^2

double Sudakov_Calculator::s_approx=Sudakov_Calculator::s_alphasfix;

//
const int Sudakov_Calculator::s_nffix=5;

//
Function_Base* Sudakov_Calculator::s_pas=NULL;

Sudakov_Calculator::Double_Double_Func
Sudakov_Calculator::GetAlphaSCorr=&Sudakov_Calculator::FixAlphaSCorr;
Sudakov_Calculator::Int_Double_Func
Sudakov_Calculator::GetNf=&Sudakov_Calculator::FixNf;



//=============================================================================



int Sudakov_Base::s_count=0;
const int& Sudakov_Base::InStore=Sudakov_Base::s_count;



//=============================================================================



Sudakov_Calculator::~Sudakov_Calculator() {    //Virtual.
  --s_count;
#ifdef TEMP_OUTPUT
  cout<<"~Sudakov_Calculator"<<endl;///////////////////////////////////////////
#endif
}





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
      GetNf=&FixNf;
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
  GetNf=&RunNf;
  return true;
}





void Sudakov_Calculator::Which() const {    //Virtual.
  cout<<"Incomplete Sudakov_Group object!"<<endl;
}



//=============================================================================



Sudakov_Base::~Sudakov_Base() {
  --s_count;
#ifdef TEMP_OUTPUT
  std::cout<<"~Sudakov_Base"<<std::endl;//////////////////////////
#endif
}





void Sudakov_Base::Which() const {
  cout<<"Incomplete Sudakov object!"<<endl;
}



//=============================================================================



template class Sudakov_Group<Dipole::qqbar>;
template class Sudakov_Group<Dipole::qg>;
template class Sudakov_Group<Dipole::gqbar>;
template class Sudakov_Group<Dipole::gg>;



//=============================================================================



template class Sudakov<Dipole::qqbar,Radiation::gluon>;
template class Sudakov<Dipole::qg,Radiation::gluon>;
template class Sudakov<Dipole::gqbar,Radiation::gluon>;
template class Sudakov<Dipole::gg,Radiation::gluon>;

template class Sudakov<Dipole::qg,Radiation::quark>;
template class Sudakov<Dipole::gqbar,Radiation::quark>;
template class Sudakov<Dipole::gg,Radiation::quark>;



//=============================================================================





//eof
