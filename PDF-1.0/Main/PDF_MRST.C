#include "PDF_MRST.H"
#include "Message.H"

using namespace PDF;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
c_mrst * PDF_MRST::proton = 0;


PDF_MRST::PDF_MRST(int isanti) {
  if (isanti) {
    msg.Tracking()<<"initialising Anti-Proton"<<std::endl;
    anti=-1;
  }
  else {
    msg.Tracking()<<"initialising Proton"<<std::endl;
    anti=+1;
  }
  if (proton) {
    msg.Tracking()<<" using already initialised C++ PDF Library MRST99... "<<std::endl;
  } 
  else {  
    msg.Tracking()<<" initialising C++ PDF Library MRST99... "<<std::endl;
    proton = new c_mrst;
    
    msg.Tracking()<<" done "<<std::endl;
  }
  for (int i=1;i<6;i++) {
    partons.push_back(Flavour(kf::code(i)));
    partons.push_back(Flavour(kf::code(i)).bar());
  }
  partons.push_back(Flavour(kf::gluon));
  partons.push_back(Flavour(kf::jet));
  partons.push_back(Flavour(kf::quark));
  partons.push_back(Flavour(kf::quark).bar());
};
