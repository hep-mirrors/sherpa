#include "PDF_Handler.H"
#include "PDF_Electron.H"
#include "PDF_MRST99.H"
#include "GRVph_Fortran_Interface.H"
#include "LHAPDF_Fortran_Interface.H"
#include "CTEQ6_Fortran_Interface.H"
#include "Exception.H"
#include "Running_AlphaS.H"
#include "Doubly_Unintegrated_PDF.H"
#include <stdio.h>

using namespace PDF;
using namespace ATOOLS;
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
      throw(ATOOLS::Exception(ATOOLS::ex::fatal_error,
			      "Tried to initialize a structure function for an uncharged particle.",
			      "PDF_Handler","GetPDFLib"));
    }
    if (bunch_particle.IsPhoton()) {
      msg.Out()<<"PDF_Handler::GetPDFLib : Try to initialize photon PDF."<<endl;
      return new GRVph_Fortran_Interface(bunch_particle);
    }
    if ((bunch_particle==Flavour(kf::p_plus) || (bunch_particle==Flavour(kf::p_plus).Bar()))) {
      PDF_Base *pdfbase=NULL;
      Switch::code kmr      = dataread->GetValue<Switch::code>("KMR_DUPDF");
      std::string set       = dataread->GetValue<string>("PDF_SET",std::string("MRST99"));
      std::string grid_path = dataread->GetValue<string>("PDF_GRID_PATH",std::string("MRST99Grid"));
      int         version   = dataread->GetValue<int>("PDF_SET_VERSION",1);
      if (set==std::string("MRST99")) {
	msg.Tracking()<<"Initialize MRST99 : "<<version<<" from "<<grid_path<<endl;
	pdfbase = new PDF_MRST99(bunch_particle,version,grid_path);
      }
      else if ((set==std::string("cteq6m") ||
	  set==std::string("cteq6d") ||
	  set==std::string("cteq6l") ||
	  set==std::string("cteq6l1")) && grid_path==std::string("CTEQ6Grid") ) {
	  
	  msg.Tracking()<<"Initialize CTEQ6 : "<<version<<" from "<<grid_path<<endl;
	  msg.Tracking()<<"Initialize CTEQ6_Fortran_Interface : "<<set<<"/"<<version<<" from "<<grid_path<<endl;
	  pdfbase = new CTEQ6_Fortran_Interface(bunch_particle,set,version,grid_path);
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
	pdfbase = new LHAPDF_Fortran_Interface(bunch_particle,set,version,grid_path,m_initlhapdf);
      }
      if (pdfbase!=NULL) {
	double mu0=dataread->GetValue("KMR_KPERP_CUT",(double)1.0);
	int kpscheme=dataread->GetValue("KMR_KPERP_SCHEME",(int)0);
 	if (kmr==Switch::On) return new Doubly_Unintegrated_PDF(pdfbase,MODEL::as,
								mu0*mu0,(kps::type)kpscheme);
	return pdfbase;
      }
      throw(ATOOLS::Exception(ATOOLS::ex::fatal_error,
			      "Combination of set/member/path for proton not properly specified.",
			      "PDF_Handler","GetPDFLib"));
    }
    throw(ATOOLS::Exception(ATOOLS::ex::fatal_error,"So far no PDF for the bunch_particle.",
			    "PDF_Handler","GetPDFLib"));
  }

  return NULL;
}
