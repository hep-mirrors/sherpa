//bof
//Version: 1 ADICIC++-0.0/2005/01/31

//Implementation of Xio.H.



#include "Xio.H"


#ifdef __GNUC__
#if __GNUC__ > 2
#include <ios>
#include <iomanip>
#include <sstream>
using std::ios;
using std::istream;
using std::string;
//using std::ostrstream;
//using std::ostringstream;
using std::stringstream;
using std::cin;
using std::cout;
using std::endl;
using std::ends;
#else
#define GCC_295 GCC_295
#include <iomanip.h>
#include <strstream.h>
#endif
#endif





//Implementation of formatted-output class sformat:
//=================================================

sformat::sformat(unsigned short w, unsigned short d, bool a, bool c)
  : wth(w), adj(a), dgt(d), cmm(c), adjust(adj) {

  //digits is the crucial option
  if(wth>80) wth=80; if(dgt>77) dgt=77;
  if(!dgt) cmm=true;  
  if(cmm) { if(wth<dgt+3) wth=dgt+3;} else { if(wth<dgt+2) wth=dgt+2;}
}





unsigned short& sformat::width(unsigned short w) {
  if(!w) return wth;
  if(w>80) wth=80; else wth=w;
  if(cmm) { if(wth<dgt+3) wth=dgt+3; return wth;} 
  else { if(wth<dgt+2) wth=dgt+2; return wth;}
}





unsigned short& sformat::digits(short d) {
  if(d<0) return dgt;
  if(d>77) dgt=77; else dgt=d;
  if(!dgt) cmm=true;
  if(cmm) { if(wth<dgt+3) wth=dgt+3;} 
  else { if(wth<dgt+2) wth=dgt+2;}
  return dgt;
}





bool& sformat::comma(char c) {
  if(c<0) return cmm;
  if(!c && dgt) cmm=false; else cmm=true;
  if(cmm) { if(wth<dgt+3) wth=dgt+3;}
  return cmm;
}





string sformat::operator()(const char c) const {
  string S; S.reserve(80);

  if(cmm) {
    S.append(dgt,c);
    if(adj) S.insert(0, string(wth-dgt, '°'));
    else S.append(wth-dgt, '°');
  }
  else S=string(wth,c);
  
  return S;
}





string sformat::operator()(const double& D) const {

#ifdef GCC_295
  const short size=80;
  char outbuf[size];
  ostrstream ostr(outbuf,size);
#else
  stringstream ostr;
#endif

  if(cmm) ostr.setf(ios::fixed|ios::showpoint,ios::floatfield);
  ostr<<std::setprecision(dgt)<<D<<ends;

#ifdef GCC_295
  string S(outbuf);
  int white=wth-S.length();
#else
  string S; ostr>>S;
  int white=wth-S.length()+1;
#endif

  if(white>0) {
    if(adj) S.insert(0, string(white,'°'));
    else S.append(white, '°');
  }

  return S;

}





//Global definitions assigned to formatted-output class sformat:
//==============================================================

//the standard sformat element
sformat sform(10,6,false,false);

//the long standard sformat element
sformat sform65(11,6,false,false);

//the short version sformat element
sformat sform44(8,4,false,false);

//the middle version sformat element
sformat sform54(9,5,false,false);

//the long middle version sformat element
sformat sform55(10,5,false,false);

//the long version sformat element
sformat sform74(11,7,false,false);

//the long long standard sformat element
sformat sform66(12,6,false,false);

//the long long version sformat element
sformat sform75(12,7,false,false);










//An istream manipulator without parameter:
//=========================================

istream& enter(istream& is) {
  char lin[2];
  do { 
    cout<<" \e[1mpress enter\e[0m ";
    is.getline(lin,1,'\n');
  }  
  while(*lin);
  return is;
}





//eof

