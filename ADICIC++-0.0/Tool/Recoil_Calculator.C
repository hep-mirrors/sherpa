//bof
//Version: 1 ADICIC++-0.0/2004/05/11

//Implementation of Recoil_Calculator.H.



#include "Recoil_Calculator.H"
//#include "Recoil_Calculator.dat.cc"





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



template class Recoil<Kleiss_Strategy>;
template class Recoil<FixDir1_Strategy>;
template class Recoil<FixDir3_Strategy>;
template class Recoil<MinimizePt_Strategy>;
template class Recoil<Lonnblad_Strategy>;
template class Recoil<OldAdicic_Strategy>;
template class Recoil<Test_Strategy>;



//=============================================================================



bool TEMP::CPTEST=false;



//=============================================================================





//eof
