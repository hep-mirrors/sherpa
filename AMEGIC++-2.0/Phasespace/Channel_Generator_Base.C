#include "Channel_Generator_Base.H"
#include "Topology.H"
#include "Point.H"

using namespace AMEGIC;
using namespace ATOOLS;

Channel_Generator_Base::Channel_Generator_Base(int _nin,int _nout,Flavour * _fl,
                                     Point * _plist) 
  : nin(_nin), nout(_nout)
{
  Topology top;
  plist  = new Point[2*(nout+1)];
  int ll = 0;
  top.Copy(_plist,plist,ll);
}

Channel_Generator_Base::~Channel_Generator_Base() { delete[] plist; }
