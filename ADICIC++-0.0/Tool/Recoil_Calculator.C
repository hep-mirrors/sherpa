//bof
//Version: 2 ADICIC++-0.0/2004/08/06

//Implementation of Recoil_Calculator.H.



#include "Recoil_Calculator.H"





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





#include "Recoil_Calculator.tpt.cc"





//=============================================================================



//So far there is no static Recoil_Calculator.
int Recoil_Calculator::s_count=0;
const Vec4D Recoil_Calculator::s_zaxis=Vec4D(1.0,0.0,0.0,1.0);

const int&   Recoil_Calculator::InStore=Recoil_Calculator::s_count;
const Vec4D& Recoil_Calculator::ZAxis=Recoil_Calculator::s_zaxis;


//=============================================================================



//If a compiler crash is preferred then comment it out.
template class Recoil<Recoil_Strategy::Unknown>;

template class Recoil<Recoil_Strategy::Kleiss>;
template class Recoil<Recoil_Strategy::FixDir1>;
template class Recoil<Recoil_Strategy::FixDir3>;
template class Recoil<Recoil_Strategy::MinimizePt>;
template class Recoil<Recoil_Strategy::Lonnblad>;
template class Recoil<Recoil_Strategy::OldAdicic>;
template class Recoil<Recoil_Strategy::Test>;



//=============================================================================





//eof
