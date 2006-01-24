//bof
//Version: 3 ADICIC++-0.0/2005/07/25

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



Recoil_Calculator::~Recoil_Calculator() {    //Virtual.
  --s_count;
#ifdef RECOIL_CALCULATOR_OUTPUT
  cout<<"  ~Recoil_Calculator\n";
#endif
}



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

template class Recoil<Recoil_Strategy::Ktii>;



//=============================================================================



void ADICIC::MakeRecos(const std::vector<bool>& gate,
		       std::vector<Recoil_Calculator*>& reco) {
  //Has to have the same order as Recoil_Strategy::List.
  assert(gate.size()==reco.size() &&
	 gate.size()>=Recoil_Strategy::NumberOfTypes-1);
  if(gate[0]) assert(reco[0]=new Recoil<Recoil_Strategy::Unknown>);
  if(gate[1]) assert(reco[1]=new Recoil<Recoil_Strategy::Kleiss>);
  if(gate[2]) assert(reco[2]=new Recoil<Recoil_Strategy::FixDir1>);
  if(gate[3]) assert(reco[3]=new Recoil<Recoil_Strategy::FixDir3>);
  if(gate[4]) assert(reco[4]=new Recoil<Recoil_Strategy::MinimizePt>);
  if(gate[5]) assert(reco[5]=new Recoil<Recoil_Strategy::Lonnblad>);
  if(gate[6]) assert(reco[6]=new Recoil<Recoil_Strategy::OldAdicic>);
  if(gate[7]) assert(reco[7]=new Recoil<Recoil_Strategy::Test>);
  if(gate[8]) assert(reco[8]=new Recoil<Recoil_Strategy::Ktii>);
  if(gate[9]) assert(reco[9]=new Recoil<Recoil_Strategy::stop>);
}



//=============================================================================





//eof
