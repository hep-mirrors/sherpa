#include "AMEGIC++/Phasespace/Channel_Generator_Base.H"
#include "AMEGIC++/Main/Topology.H"
#include "AMEGIC++/Main/Point.H"
#ifdef __GNUC__
#if __GNUC__ < 3
#include <stdio.h>
#endif
#endif

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Channel_Generator_Base::Channel_Generator_Base(int _nin,int _nout,
					       Point * _plist) 
  : nin(_nin), nout(_nout), m_valid(1)
{
  Topology top;
  plist  = new Point[2*(nout+1)];
  int ll = 0;
  top.Copy(_plist,plist,ll);
}

Channel_Generator_Base::~Channel_Generator_Base() { delete[] plist; }


string Channel_Generator_Base::GetMassIndex(string &str)
{
  char c = str[0];
  c<58 ? c-=48 : c-=55;
  char hc[4];
  sprintf(hc,"%i",c);
  return string(hc);
}

string Channel_Generator_Base::GetMassIndex(char &str)
{
  char c = str;
  c<58 ? c-=48 : c-=55;
  char hc[4];
  sprintf(hc,"%i",c);
  return string(hc);
}

