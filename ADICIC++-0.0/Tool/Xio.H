//bof
//Version: ADICIC++-0.0/2004/01/29

//A collection of little helper classes providing additional in-/out-put tools.



#ifndef _Xio_H_
#define _Xio_H _Xio_H_


#include <iostream>
#include <string>





//namespace ... {



//=========================================================
//  The class sformat - a standard formatted-output class
//=========================================================


class sformat {
  //Formatted-output class converting double to string.
private:
  unsigned short wth;    //width
  bool adj;    //adjustfield(false..left, true..right)
  unsigned short dgt;    //number of digits and ...
  bool cmm;    //ignoring the comma (false), otherwise (true)
public:
  sformat(unsigned short w=8, unsigned short d=6, bool a=false, bool c=false);
  ~sformat() {}

  bool& adjust;

  inline void status() const;
  unsigned short& width(unsigned short w=0);
  unsigned short& digits(short d=-1);
  bool& comma(char c=-1);
  std::string operator()(const char c) const;
  std::string operator()(const double& D) const;
};



inline void sformat::status() const {
  std::cout<<"sformat element with:  width="<<wth<<"  digits="<<dgt
	   <<"  adjust="<<adj<<"  comma="<<cmm<<std::endl;
}



//The standard formatted-output elements sformat:
//===============================================

extern sformat sform;
extern sformat sform65;
extern sformat sform44;
extern sformat sform54;
extern sformat sform55;
extern sformat sform74;
extern sformat sform66;
extern sformat sform75;










//===========================================================================
//  An istream manipulator without parameter (a function for class istream)
//===========================================================================


std::istream& enter(std::istream&);





#endif    //eo _Xio_H_



//eof





