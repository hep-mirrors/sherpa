//bof
//Version: 4 ADICIC++-0.0/2006/06/02

//Implementation of Sudakov_Utils.hpp.



#include <iomanip>
#include "Sudakov_Utils.hpp"
#include "Dipole_Parameter.H"





//using namespace





ADICIC::sr::stringmap::stringmap() : m_map() {
  m_map[xpini]  = "x+ini=";
  m_map[xmini]  = "x-ini=";
  m_map[fasc]   = "inimuf=";
  m_map[mdip]   = "mdip=";
  m_map[shatmax]= "-/shatmax=";
  m_map[expydip]= "e^ydip=";
  m_map[shat]   = "-/shat=";
  m_map[kt]     = "kt=";
  m_map[mt]     = "mt=";
  m_map[expy]   = "e^y=";
  m_map[xpfin]  = "x+fin=";
  m_map[xmfin]  = "x-fin=";
  m_map[stop]   = "stop !";
}

const std::string& ADICIC::sr::stringmap::operator[](code c) {
  return m_map[c];
}





ADICIC::Sudakov_Flavour::Sudakov_Flavour()
  : Glu(NULL), Qua(NULL), Aqu(NULL) {}


void ADICIC::Sudakov_Flavour::Print() const {
  if(Glu) std::cout<<(*Glu)()<<","; else std::cout<<"0,";
  if(Qua) std::cout<<(*Qua)()<<","; else std::cout<<"0,";
  if(Aqu) std::cout<<(*Aqu)(); else std::cout<<"0";
}





ADICIC::Sudakov_Result::Sudakov_Result()
  : Rad(Radiation::incorrect), Sfc(),
    P2t(0.0), Y(0.0), X1(1.0), X3(1.0), Isr(), Dir(true) {}


void ADICIC::Sudakov_Result::Print() const {
  std::cout<<setiosflags(std::ios::left)
	   <<"Rad="<<std::setw(5)<<Rad
	   <<"    P2t="<<std::setw(12)<<P2t
	   <<"    Y="<<std::setw(12)<<Y
	   <<"    X1="<<std::setw(8)<<X1
	   <<"    X3="<<std::setw(8)<<X3
	   <<"    Sfc="; Sfc.Print();
  std::cout<<"\n"<<resetiosflags(std::ios::left);
  if(Isr.empty()) return;
  std::cout<<setiosflags(std::ios::left)<<"       ";
  size_t stop1=sr::stop;
  size_t stop2=Isr.size();
  if(stop2<stop1) stop1=stop2;
  for(size_t i=0; i<stop1; ++i) {
    if(i==sr::shatmax || i==sr::shat)
      std::cout<<sr::name[sr::code(i)]<<std::setw(18)<<sqrt(Isr[i]);
    else
      std::cout<<sr::name[sr::code(i)]<<std::setw(18)<<Isr[i];
  }
  for(size_t i=stop1; i<stop2; ++i)
    std::cout<<std::setw(18)<<Isr[i];
  std::cout<<"       Dir="<<std::setw(1)<<Dir
	   <<"\n"<<resetiosflags(std::ios::left);
}





//eof
