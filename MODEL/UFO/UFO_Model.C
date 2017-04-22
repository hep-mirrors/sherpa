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

  // Overwrite masses of SM particles if they are
  // zero in UFO. Necessary for hadronization,
  // running couplings etc. Respect zero UFO masses
  // at ME level by setting 'massive' to zero in SetMassiveFlags.
  void UFO_Model::SetSMMass(const kf_code &kf,const double &m)
  {
    if (ATOOLS::s_kftable.find(kf)==ATOOLS::s_kftable.end())
      THROW(fatal_error,"SM particle not in model");
    if (ATOOLS::s_kftable[kf]->m_mass) return;
    ATOOLS::s_kftable[kf]->m_mass=m;
    ATOOLS::s_kftable[kf]->m_hmass=m;
  }

  void UFO_Model::SetSMMasses(){
    SetSMMass(kf_d,0.01);
    SetSMMass(kf_u,0.005);
    SetSMMass(kf_s,0.2);
    SetSMMass(kf_c,1.42);
    SetSMMass(kf_b,4.8);
    SetSMMass(kf_t,173.21);
    SetSMMass(kf_e,0.000511);
    SetSMMass(kf_mu,.105);
    SetSMMass(kf_tau,1.777);
  }

  // Set the massive flag consistent with UFO input.
  // Needs to be called AFTER ParamInit.
  void UFO_Model::SetMassiveFlags(){
    for (ATOOLS::KF_Table::iterator it=ATOOLS::s_kftable.begin(); it!=ATOOLS::s_kftable.end(); ++it)
      if (it->second->m_mass==0.0)
	it->second->m_massive=0;
  }

  // Set the stable flag consistent with UFO input.
  // Needs to be called AFTER ParamInit.
  void UFO_Model::SetStableFlags(){
    for (ATOOLS::KF_Table::iterator it=ATOOLS::s_kftable.begin(); it!=ATOOLS::s_kftable.end(); ++it)
      if (it->second->m_width==0.)
	it->second->m_stable=1;
  }

  bool UFO_Model::ModelInit(const PDF::ISR_Handler_Map& isr)
  { 
    std::string widthscheme = MODEL::Model_Base::p_dataread->GetValue<std::string>("WIDTH_SCHEME","Fixed");
    p_numbers->insert(make_pair(std::string("WidthScheme"), widthscheme=="CMS"));

    // set default value to UFO input such that
    // we recover standard cross sections for fixed QCD coupling
    SetAlphaQCD(isr,p_dataread->GetEntry<double>("SMINPUTS",3));

    // set default value to UFO input such that
    // we recover standard cross sections for fixed QED coupling
    SetAlphaQED(1./p_dataread->GetEntry<double>("SMINPUTS",1));
    
    return true;
  }

  Complex UFO_Model::complexconjugate(const Complex& arg) { return conj(arg); }
  Complex UFO_Model::re(const Complex& arg) { return real(arg); }
  Complex UFO_Model::im(const Complex& arg) { return imag(arg); }
  Complex UFO_Model::complex(double real, double imag) { return Complex(real, imag); }
  // Need to resolve the complex std::sqrt() /  double std::sqrt() ambiguity
  // to avoid 'nans' when double std::sqrt() is called with negative double arg
  Complex UFO_Model::sqrt(const double& arg) { return std::sqrt(Complex(arg));}
  Complex UFO_Model::sqrt(const Complex& arg) { return std::sqrt(arg);}
  // Initializing doubles with expressions involving the above sqrt
  // then requires explicit conversion
  double  UFO_Model::ToDouble(const Complex& arg){
    if (arg.imag()!=0.0)
      THROW(fatal_error, "Initializing double from complex with nonzero imaginary part");
    return arg.real();
  }

}
