#include "Zfunc.H"
#include "Message.H"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace std;

void Zfunc::Print() 
{
  msg.Out()<<"Z(["<<type<<"],";
  msg.Out()<<"[";
  for (int i=0;i<narg-1;i++) msg.Out()<<arg[i]<<";";
  
  if (narg>0) msg.Out()<<arg[narg-1];
  msg.Out()<<"][";
  msg.Out().precision(2);
  for (int i=0;i<ncoupl-1;i++) {
    if ( !AMATOOLS::IsZero(real(coupl[i])) &&
	    AMATOOLS::IsZero(imag(coupl[i])) )
      msg.Out()<<real(coupl[i])<<";";
    if (  AMATOOLS::IsZero(real(coupl[i])) &&
	  !AMATOOLS::IsZero(imag(coupl[i])) )
      msg.Out()<<imag(coupl[i])<<" I;";
    if ( !AMATOOLS::IsZero(real(coupl[i])) &&
	 !AMATOOLS::IsZero(imag(coupl[i])) )
      msg.Out()<<real(coupl[i])<<"+"<<imag(coupl[i])<<" I;";
    if (  AMATOOLS::IsZero(real(coupl[i])) &&
	  AMATOOLS::IsZero(imag(coupl[i])) )
      msg.Out()<<"0;";
  }
  if ( !AMATOOLS::IsZero(real(coupl[ncoupl-1])) &&
       AMATOOLS::IsZero(imag(coupl[ncoupl-1])) )
      msg.Out()<<real(coupl[ncoupl-1])<<"])";
  if (  AMATOOLS::IsZero(real(coupl[ncoupl-1])) &&
	!AMATOOLS::IsZero(imag(coupl[ncoupl-1])) )
    msg.Out()<<imag(coupl[ncoupl-1])<<" I])";
  if ( !AMATOOLS::IsZero(real(coupl[ncoupl-1])) &&
       !AMATOOLS::IsZero(imag(coupl[ncoupl-1])) )
    msg.Out()<<real(coupl[ncoupl-1])<<"+"<<imag(coupl[ncoupl-1])<<" I])";
  if (  AMATOOLS::IsZero(real(coupl[ncoupl-1])) &&
	AMATOOLS::IsZero(imag(coupl[ncoupl-1])) )
	msg.Out()<<"0])";
  msg.Out()<<endl;
  msg.Out().precision(6);
}

void Zfunc_Group::Print() 
{
  msg.Out()<<"SZ(["<<type<<"],";
  msg.Out()<<"[";
  for (int i=0;i<narg-1;i++) msg.Out()<<arg[i]<<";";
  
  if (narg>0) msg.Out()<<arg[narg-1];
  msg.Out()<<"])";
  msg.Out()<<endl;
  
  for (int i=0;i<zlist.size();i++) {
    if (zsign[i]==-1) {cout<<"   - "<<zlist[i]->psnew[0].numb<<" * ";zlist[i]->Print();}
                 else {cout<<"   + "<<zlist[i]->psnew[0].numb<<" * ";zlist[i]->Print();}
  }
}
