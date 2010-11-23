#ifndef AddOns_WhiteHat_WhiteHat_Interface_H
#define AddOns_WhiteHat_WhiteHat_Interface_H

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "AddOns/WhiteHat/WhiteHat_Virtual.H"
#include "AddOns/WhiteHat/WhiteHat_Tree.H"

namespace WHITEHAT {

  class WhiteHat_Interface: public PHASIC::ME_Generator_Base {
  private:

    BH::BH_interface  *p_interface;
    MODEL::Model_Base *p_model;

  public :

    // constructor
    WhiteHat_Interface();

    // destructor
    ~WhiteHat_Interface();

    // member functions
    bool Initialize(const std::string &path,const std::string &file,
		    MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beam,
		    PDF::ISR_Handler *const isr);
    PHASIC::Process_Base *InitializeProcess(const PHASIC::Process_Info &pi, bool add);
    bool PerformTests();

    void SetClusterDefinitions(PDF::Cluster_Definitions_Base *const defs);

    ATOOLS::Cluster_Amplitude *ClusterConfiguration
    (PHASIC::Process_Base *const proc,const size_t &mode,const double &kt2);

  }; // end of class WhiteHat_Interface

} // end of namespace WHITEHAT

#endif

#include "MODEL/Main/Model_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"

using namespace WHITEHAT;
using namespace PHASIC;
using namespace ATOOLS;

WhiteHat_Interface::WhiteHat_Interface(): 
  ME_Generator_Base("WhiteHat"), p_interface(NULL)
{
}

WhiteHat_Interface::~WhiteHat_Interface() 
{
  if (p_interface) delete p_interface;
}

bool WhiteHat_Interface::Initialize
(const std::string &path,const std::string &file,MODEL::Model_Base *const model,
 BEAM::Beam_Spectra_Handler *const beam,PDF::ISR_Handler *const isrhandler)
{
  if (p_interface==NULL) {
    msg_Info()<<"Initialising WhiteHat interface {"<<std::endl;
    p_model=model;
    WhiteHat_Tree::SetModel(p_model);
    WhiteHat_Virtual::SetModel(p_model);
    Data_Reader reader(" ",";","!","=");
    p_interface=new BH::BH_interface
      (reader.GetValue<std::string>("BH_SETTINGS_FILE",std::string("")));
    p_interface->set("Z_mass",Flavour(kf_Z).Mass());
    p_interface->set("Z_width",Flavour(kf_Z).Width());
    p_interface->set("W_mass",Flavour(kf_Wplus).Mass());
    p_interface->set("W_width",Flavour(kf_Wplus).Width());
    double sin_th_2=model->ScalarConstant(std::string("sin2_thetaW"));
    p_interface->set("sin_th_2",sin_th_2);
    p_interface->set("sin_2th",sin(2.*asin(sqrt(sin_th_2))));
    p_interface->set("alpha_S",model->ScalarFunction(std::string("alpha_S")));
    p_interface->set("alpha_QED",model->ScalarFunction(std::string("alpha_QED")));
    msg_Info()<<"}"<<std::endl;
    WhiteHat_Tree::SetInterface(p_interface);
    WhiteHat_Virtual::SetInterface(p_interface);
  }
  return true;
}

Process_Base *WhiteHat_Interface::InitializeProcess(const Process_Info &pi, bool add)
{
  return NULL;
}

bool WhiteHat_Interface::PerformTests()
{
  return true;
}
  
void WhiteHat_Interface::SetClusterDefinitions
(PDF::Cluster_Definitions_Base *const defs)
{
}

Cluster_Amplitude *WhiteHat_Interface::ClusterConfiguration
(Process_Base *const proc,const size_t &mode,const double &kt2)
{
  return NULL;
}

namespace PHASIC {

  DECLARE_GETTER(WhiteHat_Interface_Getter,"WhiteHat",ME_Generator_Base,ME_Generator_Key);

  ME_Generator_Base *WhiteHat_Interface_Getter::operator()(const ME_Generator_Key &key) const
  {
    return new WhiteHat_Interface();
  }

  void WhiteHat_Interface_Getter::PrintInfo(std::ostream &str,const size_t width) const
  { 
    str<<"Interface to the WhiteHat loop ME generator"; 
  }

}
