//bof
//Version: 3 ADICIC++-0.0/2005/09/21

//Implementation of Recoil_Strategy.hpp.



#include "Recoil_Strategy.hpp"





//using namespace





ADICIC::Recoil_Result::Recoil_Result()
  : Poc(both), Vec() {}


void ADICIC::Recoil_Result::Print() const {
  std::cout<<"Partner(top..1,both..0,bot..-1)="<<Poc<<"\n";
  if(Vec.empty()) return;
  std::cout<<"Four vectors: p"<<rr::code(0)+1<<"="<<Vec[0]<<"\n";
  size_t stop1=rr::stop;
  size_t stop2=Vec.size();
  if(stop2<stop1) stop1=stop2;
  for(size_t i=1; i<stop1; ++i)
    std::cout<<"              p"<<rr::code(i)+1<<Vec[i]<<"\n";
  for(size_t i=stop1; i<stop2; ++i)
    std::cout<<"               "<<rr::code(i)+1<<Vec[i]<<"\n";
}





ADICIC::Recoil_Tool::Recoil_Tool(rdt::code c)
  : m_mode(c), p_fly(NULL), p_flyprime(NULL) {}


ADICIC::Recoil_Tool::~Recoil_Tool() {
  if(p_fly) { delete p_fly; p_fly=NULL;}
  if(p_flyprime) { delete p_flyprime; p_flyprime=NULL;}
}





bool ADICIC::TEMP::CPTEST=false;





//eof
