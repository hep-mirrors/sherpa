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
  std::string number     = string(help); 
  bunch_particle = dataread->GetValue<Flavour>("BUNCH_"+number);

  if (dataread->GetValue<Switch::code>("ISR_"+number)) {    
    if (bunch_particle.IsLepton()) {
      if (bunch_particle.IntCharge()!=0) {
	return new PDF_Electron(bunch_particle.IsAnti(),
				dataread->GetValue<int>("ISR_E_ORDER"),
				dataread->GetValue<int>("ISR_E_SCHEME"));
      }
      msg.Error()<<"Error in PDF_Handler::GetPDFLib :"<<endl
		 <<"   Tried to initialize a structure function for an uncharged particle."<<endl
		 <<"   Will abort the program."<<endl;
      abort();
    }
    if ((bunch_particle==Flavour(kf::p_plus) || (bunch_particle==Flavour(kf::p_plus).Bar()))) {
      std::string set = dataread->GetValue<string>("PDF_SET");
      if (set==std::string("MRST99")) {
	return new PDF_MRST99(bunch_particle,
			      dataread->GetValue<int>("PDF_SET_VERSION"),
			      dataread->GetValue<string>("PDF_GRID_PATH"));
      }
      else if (set==std::string("Alekhin_100") ||
	       set==std::string("Alekhin_1000") ||
	       set==std::string("Botje_100") ||
	       set==std::string("Botje_1000") ||
	       set==std::string("Fermi_2002_100") ||
	       set==std::string("Fermi_2002__100") ||
	       set==std::string("MRST2001") ||
	       set==std::string("MRST98") ||
	       set==std::string("cteq6") ) {
	return new LHAPDF_Fortran_Interface(bunch_particle,set,
					    dataread->GetValue<int>("PDF_SET_VERSION"),
					    dataread->GetValue<string>("PDF_GRID_PATH"));
      }
      msg.Error()<<"Error in PDF_Handler::GetPDFLib :"<<endl
		 <<"   Combination of set/member/path for proton not properly specified :"<<endl
		 <<"   ("<<set<<"/"<<dataread->GetValue<int>("PDF_SET_VERSION")
		 <<"/"<<dataread->GetValue<string>("PDF_GRID_PATH")<<")"<<endl
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
