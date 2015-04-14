#include "ATOOLS/Org/Exception.H"
#include "MODEL/UFO/UFO_Model.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"

namespace UFO{

  UFO_Model::UFO_Model(std::string path, std::string file, bool elementary) : Model_Base(path, file, elementary) 
  {
    p_numbers          = new MODEL::ScalarNumbersMap();
    p_constants        = new MODEL::ScalarConstantsMap();
    p_complexconstants = new MODEL::ComplexConstantsMap();
    p_functions        = new MODEL::ScalarFunctionsMap();

    ATOOLS::Data_Reader* run_read = new ATOOLS::Data_Reader(" ",";","#","=");
    run_read->SetInputPath(path);
    run_read->SetInputFile(file);
    p_dataread = new UFO::UFO_Param_Reader(run_read->GetValue<std::string>("UFO_PARAM_CARD",""));
    delete run_read;
    ATOOLS::rpa->gen.AddCitation(1,"Sherpa's BSM features are published under \\cite{Hoche:2014kca}.");
    ATOOLS::rpa->gen.AddCitation(1,"The UFO model format is published under \\cite{Degrande:2011ua}.");
  }

  UFO_Model::~UFO_Model(){
    delete p_dataread;
  }

  bool UFO_Model::ModelInit(const PDF::ISR_Handler_Map& isr)
  { 
    std::string widthscheme = MODEL::Model_Base::p_dataread->GetValue<std::string>("WIDTH_SCHEME","Fixed");
    p_numbers->insert(make_pair(std::string("WidthScheme"), widthscheme=="CMS"));

    SetAlphaQCD(isr);
    // set default value to UFO input such that
    // we recover standard cross sections for fixed QCD coupling
    double alphaSU = p_dataread->GetEntry<double>("SMINPUTS",3);
    MODEL::as->SetDefault(alphaSU);

    double alphaU = 1./p_dataread->GetEntry<double>("SMINPUTS",1);
    double alpha  = 1./MODEL::Model_Base::p_dataread->GetValue<double>("1/ALPHAQED(0)",137.03599976);
    MODEL::aqed = new MODEL::Running_AlphaQED(alpha);
    MODEL::aqed->SetDefault(alphaU);
    p_functions->insert(make_pair(std::string("alpha_QED"),MODEL::aqed));
    p_constants->insert(make_pair(std::string("alpha_QED"),alphaU));

    return true;
  }

  const Complex UFO_Model::complexconjugate(const Complex& arg) { return conj(arg); }
  const Complex UFO_Model::re(const Complex& arg) { return real(arg); }
  const Complex UFO_Model::im(const Complex& arg) { return imag(arg); }
  const Complex UFO_Model::complex(double real, double imag) { return Complex(real, imag); }

}
