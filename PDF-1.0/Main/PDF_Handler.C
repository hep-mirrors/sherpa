#include "PDF_Handler.H"
#include "PDF_Electron.H"
#include "PDF_MRST99.H"
#include "LHAPDF_Fortran_Interface.H"
#include "Message.H"
#include <stdio.h>

using namespace PDF;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace std;

PDF_Base * PDF_Handler::GetPDFLib(Data_Read * dataread,Flavour & bunch_particle,
				  const int num) {
  char help[20];
  sprintf(help,"%i",num+1);
  std::string number         = string(help); 
  int         flav           = dataread->GetValue<int>("BUNCH_"+number);  
  bunch_particle             = Flavour(kf::code(abs(flav)));
  if (flav<0) bunch_particle = bunch_particle.Bar();

  if (dataread->GetValue<Switch::code>("ISR_"+number)) {    
    if (bunch_particle.IsLepton()) {
      if (bunch_particle.IntCharge()!=0) {
	return new PDF_Electron(bunch_particle,
				dataread->GetValue<int>("ISR_E_ORDER",1),
				dataread->GetValue<int>("ISR_E_SCHEME",2));
      }
      msg.Error()<<"Error in PDF_Handler::GetPDFLib :"<<endl
		 <<"   Tried to initialize a structure function for an uncharged particle."<<endl
		 <<"   Will abort the program."<<endl;
      abort();
    }
    if ((bunch_particle==Flavour(kf::p_plus) || (bunch_particle==Flavour(kf::p_plus).Bar()))) {
      std::string set       = dataread->GetValue<string>("PDF_SET",std::string("MRST99"));
      std::string grid_path = dataread->GetValue<string>("PDF_GRID_PATH",std::string("MRST99Grid"));
      int         version   = dataread->GetValue<int>("PDF_SET_VERSION",1);
      if (set==std::string("MRST99")) {
	msg.Tracking()<<"Initialize MRST99 : "<<version<<" from "<<grid_path<<endl;
	return new PDF_MRST99(bunch_particle,version,grid_path);
      }
      else if (set==std::string("Alekhin_100") ||
	       set==std::string("Alekhin_1000") ||
	       set==std::string("Botje_100") ||
	       set==std::string("Botje_1000") ||
	       set==std::string("Fermi_2002_100") ||
	       set==std::string("Fermi_2002__100") ||
	       set==std::string("MRST2001") ||
	       set==std::string("MRST98") ||
	       set==std::string("cteq6")  ||
	       set==std::string("cteq6l")  ||
	       set==std::string("cteq6ll")) {
	msg.Tracking()<<"Initialize "<<set<<" : "<<version<<" from "<<grid_path<<endl;
	return new LHAPDF_Fortran_Interface(bunch_particle,set,version,grid_path,m_initlhapdf);
      }
      msg.Error()<<"Error in PDF_Handler::GetPDFLib :"<<endl
		 <<"   Combination of set/member/path for proton not properly specified :"<<endl
		 <<"   ("<<set<<"/"<<version<<"/"<<grid_path<<")"<<endl
		 <<"   Will abort the program."<<endl;
      abort();
    }
    msg.Error()<<"Error in PDF_Handler::GetPDFLib :"<<endl
	       <<"   So far no PDF for the bunch_particle :"<<bunch_particle<<endl
	       <<"   Will abort the program."<<endl;
    abort();
  }

  return NULL;
}
