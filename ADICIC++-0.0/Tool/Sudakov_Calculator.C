//bof
//Version: 1 ADICIC++-0.0/2004/04/30

//Implementation of Sudakov_Calculator.H.



#include "Sudakov_Calculator.H"
#include "Sudakov_Calculator.dat.cc"





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





#include "Sudakov_Calculator.tpt.cc"





//=============================================================================



//So far there is no static Sudakov_Calculator.
int Sudakov_Calculator::s_count=0;
const int& Sudakov_Calculator::InStore=Sudakov_Calculator::s_count;



//=============================================================================



template class Sudakov<Dipole::qqbar,Alpha_S_Fix>;
template class Sudakov<Dipole::qg,Alpha_S_Fix>;
template class Sudakov<Dipole::gqbar,Alpha_S_Fix>;
template class Sudakov<Dipole::gg,Alpha_S_Fix>;



//=============================================================================





//eof
