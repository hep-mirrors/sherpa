#include "XS_Base.H"
#include "Run_Parameter.H"
#include "Message.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace AORGTOOLS;
using namespace BEAM;
using namespace ISR;
using namespace std;

XS_Base::~XS_Base() {
  if (fl)      { delete [] fl; fl = 0; }
  if (colours) { 
    for (int i=0;i<nin+nout;i++) delete colours[i];
    delete [] colours; colours = 0;
  }
};

void XS_Base::Init(int _nin,int _nout,Flavour * _fl) {
  msg.Debugging()<<"In XS_Base::Init("<<_nin<<","<<_nout<<")"<<endl;
  nin     = _nin; nout = _nout;
  fl      = new Flavour[nin+nout];
  if (_fl) {
    for (short int i=0;i<nin+nout;i++) {
      fl[i]  = _fl[i];
      msg.Debugging()<<fl[i]<<" ";
    }
  }
  msg.Debugging()<<endl;

  colours = new int*[nin+nout];
  for (int i=0;i<nin+nout;i++) colours[i] = new int[2];
}

double XS_Base::operator()(double,double,double) {
  AORGTOOLS::msg.Error()<<"Virtual Method : XS_Base::operator()."<<std::endl; 
  return 0.;
}

bool XS_Base::SetColours(double,double,double) {
  AORGTOOLS::msg.Error()<<"Virtual Method : XS_Base::SetColours()."<<std::endl; 
  return 0.; 
}

bool XS_Base::SetColours(AMATOOLS::Vec4D * p) {
  s = (p[0]+p[1]).Abs2();
  t = (p[0]-p[2]).Abs2();
  u = (p[0]-p[3]).Abs2();
  return SetColours(s,t,u);
}

XS_Base * XS_Base::Selected()  { 
  msg.Debugging()<<"Error : Virtual method XS_Base::Selected()"<<endl;
}



