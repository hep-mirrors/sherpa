#include "ATOOLS/Org/Exception.H"
#include "MODEL/UFO/UFO_Model.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Flavour.H"
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
    // this part hopes the UFO model has the SM with the same kfcodes included
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
    bool cms(widthscheme=="CMS");
    p_numbers->insert(make_pair(std::string("WidthScheme"), cms));

    // set default value to UFO input such that
    // we recover standard cross sections for fixed QCD coupling
    // warning is printed to user to check value if consistent with UFO model
    msg_Info()<<METHOD<<"(): Trying to read in \\alpha_s(m_Z) as parameter "
              <<"3 in SMINPUTS block."<<std::endl;
    SetAlphaQCD(isr,p_dataread->GetEntry<double>("SMINPUTS",3,0.118,false));

    // set default value to UFO input such that
    // we recover standard cross sections for fixed QED coupling
    // warning is printed to user to check value if consistent with UFO model
    msg_Info()<<METHOD<<"(): Trying to read in \\alpha_QED as parameter "
              <<"1 in SMINPUTS block."<<std::endl;
    SetAlphaQED(1./p_dataread->GetEntry<double>("SMINPUTS",1,127.952,false));

    // set default value for sin(theta), cos(theta), vev if not available 
    // because the parameter is needed for the parton shower and beyond; 
    // it will be incorrect for most EFT parameter points, but 
    // as UFO does not have a canonical name for it, this is the 
    // best we can do
    // warning is printed to user to check value if consistent with UFO model
    // only fill if the W and Z are defined by the model with the usual kfcodes
    if (ATOOLS::s_kftable.find(kf_Wplus)!=ATOOLS::s_kftable.end() &&
        ATOOLS::s_kftable.find(kf_Z)!=ATOOLS::s_kftable.end()) {
      Complex I(0.,1.);
      double MW(ATOOLS::Flavour(kf_Wplus).Mass()),  MZ(ATOOLS::Flavour(kf_Z).Mass());
      double GW(ATOOLS::Flavour(kf_Wplus).Width()), GZ(ATOOLS::Flavour(kf_Z).Width());
      Complex muW2(MW*(MW-(cms?I*GW:0.))), muZ2(MZ*(MZ-(cms?I*GZ:0.)));
      Complex ccos2thetaW=muW2/muZ2;
      Complex csin2thetaW=1.-ccos2thetaW;
      Complex cvev=2.*sqrt(muW2*csin2thetaW/(4.*M_PI*MODEL::aqed->Default()));
      if (p_complexconstants->find("ccos2_thetaW")==p_complexconstants->end()) {
        msg_Info()<<METHOD<<"(): Trying to read in cos(\\theta_W) as parameter "
                  <<"10 in SMINPUTS block."<<std::endl;
        ccos2thetaW=p_dataread->GetEntry<Complex>("SMINPUTS",10,ccos2thetaW,false);
        p_complexconstants->insert(make_pair(std::string("ccos2_thetaW"),ccos2thetaW));
      }
      if (p_complexconstants->find("csin2_thetaW")==p_complexconstants->end()) {
        msg_Info()<<METHOD<<"(): Trying to read in sin(\\theta_W) as parameter "
                  <<"11 in SMINPUTS block."<<std::endl;
        csin2thetaW=p_dataread->GetEntry<Complex>("SMINPUTS",11,csin2thetaW,false);
        p_complexconstants->insert(make_pair(std::string("csin2_thetaW"),csin2thetaW));
      }
      if (p_complexconstants->find("cvev")==p_complexconstants->end()) {
        msg_Info()<<METHOD<<"(): Trying to read in vev as parameter "
                  <<"12 in SMINPUTS block."<<std::endl;
        cvev=p_dataread->GetEntry<Complex>("SMINPUTS",12,cvev,false);
        p_complexconstants->insert(make_pair(std::string("cvev"), cvev));
      }
    }
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
