#include "PDF_Handler.H"
#include "PDF_Electron.H"
#include "PDF_MRST.H"
#include "Message.H"

using namespace PDF;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace std;

PDF_Base * PDF_Handler::GetPDFLib(const Flavour fl) {
  int kfc = int(fl);
  msg.Tracking()<<"PDF_Handler::GetPDFLib called with "<<fl<<" : "<<kfc<<endl;
  switch (kfc) {
  case (kf::e) : 
    return new PDF_Electron(0);
  case (-kf::e) :  
    return new PDF_Electron(1);
  case (kf::code(2212)) :        // P- = anti-proton
    return new PDF_MRST(0);
  case -(kf::code(2212)) :       // P+ = proton
    return new PDF_MRST(1);
  default :                    // proton
    msg.Tracking()<<" Error: no PDF Library for "<<fl<<" found, assuming proton!"<<endl;
    return new PDF_MRST();
  }
}
