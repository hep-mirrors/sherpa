#ifndef AddOns_BlackHat_BlackHat_Interface_H
#define AddOns_BlackHat_BlackHat_Interface_H

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "AddOns/BlackHat/BlackHat_Virtual.H"
#include "AddOns/BlackHat/BlackHat_Tree.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"

namespace BLACKHAT {

  class BlackHat_Interface: public PHASIC::ME_Generator_Base {
  private:

    BH::BH_interface  *p_interface;
    MODEL::Model_Base *p_model;

  public :

    // constructor
    BlackHat_Interface();

    // destructor
    ~BlackHat_Interface();

    // member functions
    bool Initialize(const std::string &path,const std::string &file,
		    MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beam,
		    PDF::ISR_Handler *const isr);
    PHASIC::Process_Base *InitializeProcess(const PHASIC::Process_Info &pi, bool add);
    int  PerformTests();
    bool NewLibraries();

  }; // end of class BlackHat_Interface

} // end of namespace WHITEHAT

#endif

#include "MODEL/Main/Model_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"

using namespace BLACKHAT;
using namespace PHASIC;
using namespace ATOOLS;

BlackHat_Interface::BlackHat_Interface(): 
  ME_Generator_Base("BlackHat"), p_interface(NULL)
{
}

BlackHat_Interface::~BlackHat_Interface() 
{
  if (p_interface) delete p_interface;
}

bool BlackHat_Interface::Initialize
(const std::string &path,const std::string &file,MODEL::Model_Base *const model,
 BEAM::Beam_Spectra_Handler *const beam,PDF::ISR_Handler *const isrhandler)
{
  if (p_interface==NULL) {
    rpa->gen.AddCitation(1,"The BlackHat library is described in \\cite{Berger:2008sj}.");
    msg_Info()<<"Initialising BlackHat interface {"<<std::endl;
    p_model=model;
    BlackHat_Tree::SetModel(p_model);
    BlackHat_Virtual::SetModel(p_model);
    Data_Reader reader(" ",";","!","=");
    p_interface=new BH::BH_interface
      (reader.GetValue<std::string>("BH_SETTINGS_FILE",std::string("")));
    p_interface->set("Z_mass",Flavour(kf_Z).Mass());
    p_interface->set("Z_width",Flavour(kf_Z).Width());
    p_interface->set("W_mass",Flavour(kf_Wplus).Mass());
    p_interface->set("W_width",Flavour(kf_Wplus).Width());
#ifdef INCLUDE_COUPLINGS_IN_VIRTUAL
    p_interface->set("H_mass",Flavour(kf_h0).Mass());
    p_interface->set("H_width",Flavour(kf_h0).Width());
#endif
    double sin_th_2=std::abs(model->ComplexConstant(std::string("csin2_thetaW")));
    p_interface->set("sin_th_2",sin_th_2);
    p_interface->set("alpha_S",model->ScalarConstant("alpha_S"));
    p_interface->set("alpha_QED",model->ScalarConstant("alpha_QED"));
#ifdef INCLUDE_COUPLINGS_IN_VIRTUAL
    p_interface->set("YUK2",1.0/std::abs(sqr(model->ComplexConstant("cvev"))));
#endif
    msg_Info()<<"}"<<std::endl;
    BlackHat_Tree::SetInterface(p_interface);
    BlackHat_Virtual::SetInterface(p_interface);
  }
  return true;
}

Process_Base *BlackHat_Interface::InitializeProcess(const Process_Info &pi, bool add)
{
  return NULL;
}

int BlackHat_Interface::PerformTests()
{
  return 1;
}
  
bool BlackHat_Interface::NewLibraries()
{
  return false;
}

DECLARE_GETTER(BlackHat_Interface,"BlackHat",ME_Generator_Base,ME_Generator_Key);

ME_Generator_Base *ATOOLS::Getter<ME_Generator_Base,ME_Generator_Key,
				  BlackHat_Interface>::
operator()(const ME_Generator_Key &key) const
{
  return new BlackHat_Interface();
}

void ATOOLS::Getter<ME_Generator_Base,ME_Generator_Key,BlackHat_Interface>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Interface to the BlackHat loop ME generator"; 
}
