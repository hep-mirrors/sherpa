#include "PDF_Handler.H"
#include "PDF_Electron.H"
#include "PDF_MRST99.H"
#include "PDF_MRST01LO.H"
#include "GRVph_Fortran_Interface.H"

#ifdef USING__LHAPDF
#include "LHAPDF_Fortran_Interface.H"
#else
#include "CTEQ6_Fortran_Interface.H"
#endif

#include "Exception.H"
#include "Running_AlphaS.H"
#include "Doubly_Unintegrated_PDF.H"
#include "Continued_PDF.H"
#include "Run_Parameter.H"
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
      THROW(fatal_error,"Tried to initialize a structure function for an uncharged particle.");
    }
    if (bunch_particle.IsPhoton()) {
      msg.Out()<<"PDF_Handler::GetPDFLib : Try to initialize photon PDF."<<endl;
      return new GRVph_Fortran_Interface(bunch_particle);
    }
    if ((bunch_particle==Flavour(kf::p_plus) || (bunch_particle==Flavour(kf::p_plus).Bar()))) {
      PDF_Base *pdfbase=NULL;
      Switch::code kmr      = dataread->GetValue<Switch::code>("KMR_DUPDF");
      Switch::code cont     = dataread->GetValue<Switch::code>("CONTINUE_PDF");
      std::string set       = dataread->GetValue<string>("PDF_SET",std::string("MRST99"));
      std::string grid_path = dataread->GetValue<string>("PDF_GRID_PATH",std::string("MRST99Grid"));
      int         version   = dataread->GetValue<int>("PDF_SET_VERSION",1);
      grid_path=ATOOLS::rpa.gen.Variable("SHERPA_PDF_PATH")+std::string("/")+grid_path;
      if (set==std::string("MRST99")) {
	msg_Tracking()<<"Initialize MRST99 : "<<version<<" from "<<grid_path<<endl;
	pdfbase = new PDF_MRST99(bunch_particle,version,grid_path);
      }
      else if (set==std::string("MRST01LO")) {
	msg_Tracking()<<"Initialize MRST01LO from "<<grid_path<<endl;
	pdfbase = new PDF_MRST01LO(bunch_particle,grid_path);
      }
      else if ((set==std::string("cteq6m") ||
		set==std::string("cteq6d") ||
		set==std::string("cteq6l") ||
		set==std::string("cteq6l1")) && grid_path.find("CTEQ6Grid") ) {
	  
#ifdef USING__LHAPDF
	msg.Error()<<"ERROR : Cannot initialize CTEQ6 interface when LHAPDF is enabled ! "<<std::endl;
#else	
	msg_Tracking()<<"Initialize CTEQ6 : "<<version<<" from "<<grid_path<<endl;
	msg_Tracking()<<"Initialize CTEQ6_Fortran_Interface : "<<set<<"/"<<version<<" from "<<grid_path<<endl;
	pdfbase = new CTEQ6_Fortran_Interface(bunch_particle,set,version,grid_path);
#endif
      }
      else if (set==std::string("Alekhin_100.LHpdf") ||
	       set==std::string("Alekhin_1000.LHpdf") ||
	       set==std::string("a02_lo_v.LHgrid") ||
	       set==std::string("a02_nlo_v.LHgrid") ||
	       set==std::string("a02_nnlo_v.LHgrid") ||
	       set==std::string("Botje_100.LHpdf") ||
	       set==std::string("Botje_1000.LHpdf") ||
	       set==std::string("Fermi2002_100.LHpdf") ||
	       set==std::string("Fermi2002_1000.LHpdf") ||
	       set==std::string("MRST2001.LHpdf") ||
	       set==std::string("MRST98.LHpdf") ||
	       set==std::string("MRST2001E.LHgrid") ||
	       set==std::string("MRST2001E.LHpdf") ||
	       set==std::string("MRST2001lo.LHgrid") ||
	       set==std::string("MRST2001nlo.LHgrid") ||
	       set==std::string("MRST2001nlo.LHpdf") ||
	       set==std::string("MRST2001nnlo.LHgrid") ||
	       set==std::string("MRST2002nlo.LHgrid") ||
	       set==std::string("MRST2002nlo.LHpdf") ||
	       set==std::string("MRST2002nnlo.LHgrid") ||
	       set==std::string("MRST2003cnlo.LHgrid") ||
	       set==std::string("MRST2003cnlo.LHpdf") ||
	       set==std::string("MRST2003cnnlo.LHgrid") ||
	       set==std::string("cteq6.LHpdf")  ||
	       set==std::string("cteq6l.LHpdf")  ||
	       set==std::string("cteq6ll.LHpdf") ||
	        set==std::string("cteq61.LHgrid") ||
	       set==std::string("cteq61.LHpdf") ||
	       set==std::string("cteq6mE.LHgrid") ||
	       set==std::string("cteq6m.LHpdf") ||
	       set==std::string("cteq5l.LHgrid") ||
	       set==std::string("cteq5d.LHgrid") ||
	       set==std::string("cteq5m1.LHgrid") ||
	       set==std::string("cteq5m.LHgrid") ||
	       set==std::string("cteq4l.LHgrid") ||
	       set==std::string("cteq4d.LHgrid") ||
	       set==std::string("cteq4l.LHgrid") ||
	       set==std::string("cteq4m.LHgrid") ||
	       set==std::string("GRV98lo.LHgrid") ||
	       set==std::string("GRV98nlo.LHgrid") ||
	       set==std::string("H12000disE.LHgrid") ||
	       set==std::string("H12000dis.LHgrid") ||
	       set==std::string("H12000lo2E.LHgrid") ||
	       set==std::string("H12000lo2.LHgrid") ||
	       set==std::string("H12000loE.LHgrid") ||
	       set==std::string("H12000lo.LHgrid") ||
	       set==std::string("H12000msE.LHgrid") ||
	       set==std::string("H12000ms.LHgrid") ||
	       set==std::string("ZEUS2002_FF.LHpdf") ||
	       set==std::string("ZEUS2002_TR.LHpdf") ||
	       set==std::string("ZEUS2002_ZM.LHpdf")) {
#ifdef USING__LHAPDF
	msg_Tracking()<<"Initialize LHAPDF "<<set<<" : "<<version<<" from "<<grid_path<<endl;
	pdfbase = new LHAPDF_Fortran_Interface(bunch_particle,set,version,grid_path,m_initlhapdf);
#else 
	msg.Error()<<"ERROR : USING__LHAPDF is not enabled ! "<<std::endl;
	pdfbase = NULL;
#endif
      }
      if (pdfbase!=NULL) {
	if (kmr==Switch::On) {
	  double mu0=dataread->GetValue("KMR_KPERP_CUT",(double)1.0);
	  Doubly_Unintegrated_PDF *dupdf =
	    new Doubly_Unintegrated_PDF(pdfbase,MODEL::as,mu0*mu0);
	  int kpscheme=dataread->GetValue("KMR_KPERP_SCHEME",(int)0);
	  dupdf->SetKPerpScheme((kps::type)kpscheme);
	  double exponent=dataread->GetValue("KMR_FIXED_EXPONENT",(double)0.0);
	  dupdf->SetFixedKtExponent(exponent);
	  int mode=dataread->GetValue("KMR_MODE",(int)0);
	  dupdf->SetMode(mode);
	  dupdf->Initialize();
	  return dupdf;
	}
 	if (cont==Switch::On) {
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
