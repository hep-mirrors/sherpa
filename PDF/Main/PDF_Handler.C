#include "PDF_Handler.H"
#include "PDF_Electron.H"
#include "GRVph_Fortran_Interface.H"

#include "CXXFLAGS_PACKAGES.H"
#ifdef USING__LHAPDF
#include "LHAPDF_Fortran_Interface.H"
#else
#include "CTEQ6_Fortran_Interface.H"
#include "PDF_MRST99.H"
#include "PDF_MRST01LO.H"
#endif

#include "Exception.H"
#include "Running_AlphaS.H"
#include "Continued_PDF.H"
#include "Run_Parameter.H"
#include "Data_Reader.H"
#include <stdio.h>

using namespace PDF;
using namespace ATOOLS;
using namespace std;

PDF_Base * PDF_Handler::GetPDFLib(Data_Reader * dataread,Flavour & bunch_particle,
				  const int num) {
  char help[20];
  sprintf(help,"%i",num+1);
  std::string number         = string(help); 
  int defaultflav(0);
  if (num==0) {
    defaultflav=rpa.gen.Beam1().IsAnti()? 
      -rpa.gen.Beam1().Kfcode() : rpa.gen.Beam1().Kfcode();
  }
  else if (num==1) {
    defaultflav=rpa.gen.Beam2().IsAnti()? 
      -rpa.gen.Beam2().Kfcode() : rpa.gen.Beam2().Kfcode();
  }
  int flav = dataread->GetValue<int>("BUNCH_"+number,defaultflav);
  bunch_particle             = Flavour((kf_code)abs(flav));
  if (flav<0) bunch_particle = bunch_particle.Bar();

  string defaultisr("Off");
  if (bunch_particle.Kfcode()==kf_p_plus) defaultisr="On";
  if (dataread->GetValue<std::string>("ISR_"+number,defaultisr)=="On") {    
    if (bunch_particle.IsLepton()) {
      if (bunch_particle.IntCharge()!=0) {
	return new PDF_Electron(bunch_particle,
				dataread->GetValue<int>("ISR_E_ORDER",1),
				dataread->GetValue<int>("ISR_E_SCHEME",2));
      }
      THROW(fatal_error,"Tried to initialize a structure function for an uncharged particle.");
    }
    if (bunch_particle.IsPhoton()) {
      msg_Out()<<"PDF_Handler::GetPDFLib : Try to initialize photon PDF."<<endl;
      return new GRVph_Fortran_Interface(bunch_particle);
    }
    if ((bunch_particle==Flavour(kf_p_plus) || (bunch_particle==Flavour(kf_p_plus).Bar()))) {
      PDF_Base *pdfbase=NULL;
      std::string cont     = dataread->GetValue<std::string>("CONTINUE_PDF","Off");
#ifdef USING__LHAPDF
      std::string set       = dataread->GetValue<string>("PDF_SET",std::string("cteq6l.LHpdf"));
#else
      std::string set       = dataread->GetValue<string>("PDF_SET",std::string("cteq6l"));
#endif
      std::string grid_path = dataread->GetValue<string>("PDF_GRID_PATH",std::string("CTEQ6Grid"));
      int         version   = dataread->GetValue<int>("PDF_SET_VERSION",1);
      grid_path=ATOOLS::rpa.gen.Variable("SHERPA_SHARE_PATH")+std::string("/")+grid_path;
      if (set==std::string("MRST99")) {
#ifdef USING__LHAPDF
	msg_Error()<<"ERROR : Cannot initialize MRST interface when LHAPDF "
                   <<"is enabled ! "<<std::endl;
#else
	msg_Tracking()<<"Initialize MRST99 : "<<version<<" from "<<grid_path<<endl;
	pdfbase = new PDF_MRST99(bunch_particle,version,grid_path);
#endif
      }
      else if (set==std::string("MRST01LO")) {
#ifdef USING__LHAPDF
	msg_Error()<<"ERROR : Cannot initialize MRST interface when LHAPDF "
                   <<"is enabled ! "<<std::endl;
#else
	msg_Tracking()<<"Initialize MRST01LO from "<<grid_path<<endl;
	pdfbase = new PDF_MRST01LO(bunch_particle,grid_path);
#endif
      }
      else if ((set==std::string("cteq6m") ||
		set==std::string("cteq6d") ||
		set==std::string("cteq6l") ||
		set==std::string("cteq6l1")) && grid_path.find("CTEQ6Grid") ) {
	  
#ifdef USING__LHAPDF
	msg_Error()<<"ERROR : Cannot initialize CTEQ6 interface when LHAPDF is enabled ! "<<std::endl;
#else	
	msg_Tracking()<<"Initialize CTEQ6 : "<<version<<" from "<<grid_path<<endl;
	msg_Tracking()<<"Initialize CTEQ6_Fortran_Interface : "<<set<<"/"<<version<<" from "<<grid_path<<endl;
	pdfbase = new CTEQ6_Fortran_Interface(bunch_particle,set,version,grid_path);
#endif
      }
      else if (set.find(std::string("LHpdf")) || 
	        set.find(std::string("LHgrid"))) {
#ifdef USING__LHAPDF
	msg_Tracking()<<"Initialize LHAPDF "<<set<<" : "<<version<<endl;
	pdfbase = new LHAPDF_Fortran_Interface(bunch_particle,set,version,m_initlhapdf);
#else 
	msg_Error()<<"ERROR : USING__LHAPDF is not enabled ! "<<std::endl;
	pdfbase = NULL;
#endif
      }
      if (pdfbase!=NULL) {
 	if (cont=="On") {
	  int pcscheme=dataread->GetValue("CONTINUATION",(int)1);
	  return new Continued_PDF(pdfbase,(pcs::type)pcscheme);
	}
	return pdfbase;
      }
      THROW(fatal_error,"Combination of set/member/path for proton not properly specified.");
    }
    THROW(fatal_error,"So far no PDF for the bunch_particle.");
  }

  return NULL;
}
