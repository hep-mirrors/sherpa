//bof
//Version: 4 ADICIC++-0.0/2006/05/26

//Implementation of Recoil_Utils.hpp.



#include "Recoil_Utils.hpp"





using namespace ADICIC;





Multipoincare::Multipoincare(rdt::code c)
  : m_mode(c), v_trafs() {}


Multipoincare::~Multipoincare() {
  Destruct();
}


const bool Multipoincare::Apply(ATOOLS::Vec4D& pvec) const {
  if(v_trafs.empty()) return true;
  ATOOLS::Vec4D temp(pvec);
  for(size_t i=v_trafs.size()-1; true; --i) {
    ATOOLS::Poincare& traf=*(v_trafs[i].second);
    switch(v_trafs[i].first) {
    case  0:
      if(traf.CheckBoost()) { traf.Boost(temp); break;}
      else return false;
    case  1:
      if(traf.CheckBoost()) { traf.BoostBack(temp); break;}
      else return false;
    case  2:
      if(traf.CheckRotation()) { traf.Rotate(temp); break;}
      else return false;
    case  3:
      if(traf.CheckRotation()) { traf.RotateBack(temp); break;}
      else return false;
    default: return false;
    }
    if(i==0) break;
  }
  pvec=temp;
  return true;
}


void Multipoincare::Clear() {
  //Maintain the mode.
  Destruct();
}


void Multipoincare::Clear(rdt::code c) {
  m_mode=c;
  Destruct();
}



//=============================================================================



Recoil_Result::Recoil_Result(rdt::code c)
  : Poc(both), Vec(), Mup(c) {}


void Recoil_Result::Print() const {
  std::cout<<"Partner(top..1,both..0,bot..-1)="<<Poc<<"\n";
  if(Vec.empty()) std::cout<<"No momenta.\n";
  else {
    std::cout<<"Four vectors: p"<<rr::code(0)+1<<"="
	     <<Vec[0]<<"  "<<Vec[0].Abs2()<<"\n";
    size_t stop1=rr::stop;
    size_t stop2=Vec.size();
    if(stop2<stop1) stop1=stop2;
    for(size_t i=1; i<stop1; ++i)
      std::cout<<"              p"<<rr::code(i)+1<<"="
	       <<Vec[i]<<"  "<<Vec[i].Abs2()<<"\n";
    for(size_t i=stop1; i<stop2; ++i)
      std::cout<<"               "<<rr::code(i)+1<<"="
	       <<Vec[i]<<"  "<<Vec[i].Abs2()<<"\n";
  }
  if(Mup.GetItsVec().empty()) std::cout<<"No Poincares kept."<<std::endl;
  else std::cout<<"Number of Poincares kept: "<<Mup.GetItsVec().size()
		<<"   using mode: "<<Mup.Mode()<<std::endl;
}





bool TEMP::CPTEST=false;





//eof
