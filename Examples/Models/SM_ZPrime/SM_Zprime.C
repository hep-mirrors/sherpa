#include "MODEL/Main/Model_Base.H"

#include "MODEL/Main/Standard_Model.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"

#include <iostream>
#include <iomanip>

using namespace std;
using namespace ATOOLS;

namespace MODEL {
  
  class SM_Zprime : public Model_Base {
  protected :      
    MODEL::Model_Base * p_sm;

    void ParticleInit();
    void FillSpectrum(const PDF::ISR_Handler_Map& isr);
    bool ModelInit(const PDF::ISR_Handler_Map& isr);

  public :
    SM_Zprime(std::string,std::string,bool);
  };
  
  
  SM_Zprime::SM_Zprime(std::string dir, std::string file, bool elementary) :
      Model_Base(dir,file,elementary)
  {
    p_sm = new Standard_Model(m_dir,m_file,false);

    ParticleInit();
    if (m_elementary) {
      ATOOLS::OutputParticles(msg->Info());
      ATOOLS::OutputContainers(msg->Info());
    }
  }


  void SM_Zprime::ParticleInit() {
    // Add Z' particle
    // kf_code,mass,width,3*charge,icharge,strong,2*spin,majorana,take,stable,massive,idname,tex_name
    kf_code kfZp=32;
    s_kftable[kfZp] = new Particle_Info
      (kfZp, 1000., 10., 0, 0, 0, 2, -1, 1, 0, 1, "Zprime", "Z^{\\prime}");

    // Allow overriding the Zprime data from run card or command line
    ReadParticleData();
  }

  
  void SM_Zprime::FillSpectrum(const PDF::ISR_Handler_Map& isr) {
    // Constants needed for Z'
    p_dataread = new Data_Reader(" ",";","!","=");
    p_dataread->AddWordSeparator("\t");
    p_dataread->SetInputPath(m_dir);
    p_dataread->SetInputFile(m_file);
    p_constants->insert(make_pair(string("Zp_cpl_L"),
                                  p_dataread->GetValue<double>("Zp_cpl_L",1.)));
    p_constants->insert(make_pair(string("Zp_cpl_R"),
                                  p_dataread->GetValue<double>("Zp_cpl_R",1.)));
  }

  bool SM_Zprime::ModelInit(const PDF::ISR_Handler_Map& isr)
  {
    if (m_elementary) {
      msg_Info()<<"Initialize the Standard Model plus Zprime from "
          <<m_dir<<" / "<<m_file<<std::endl;
    }
    m_name      = std::string("SM+Zprime");

    p_sm->ModelInit(isr);
    p_numbers   = p_sm->ExtractScalarNumbers();
    p_constants = p_sm->ExtractScalarConstants();
    p_complexconstants = p_sm->ExtractComplexConstants();
    p_functions = p_sm->ExtractScalarFunctions();
    p_matrices  = p_sm->ExtractComplexMatrices();

    FillSpectrum(isr);

    return true;
  }
}


// Now follows some magic to make the model known to Sherpa and print a summary

using namespace MODEL;

DECLARE_GETTER(SM_Zprime,"SM+Zprime",Model_Base,Model_Arguments);

Model_Base *ATOOLS::Getter
<Model_Base,Model_Arguments,SM_Zprime>::
operator()(const Model_Arguments &args) const
{
  return new SM_Zprime(args.m_path,args.m_file,args.m_elementary);
}

void ATOOLS::Getter<Model_Base,Model_Arguments,SM_Zprime>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The Standard Model + Zprime\n"
     <<std::setw(width+4)<<" "<<"{\n"
     <<std::setw(width+7)<<" "<<"parameter specification [keyword=value]\n"
     <<std::setw(width+7)<<" "<<"- all the SM parameters\n"
     <<std::setw(width+7)<<" "<<"- mass/width of the Zprime, via MASS[32] & WIDTH[32]\n"
     <<std::setw(width+7)<<" "<<"- multiplicative coupling parameter Zp_cpl_L\n"
     <<std::setw(width+7)<<" "<<"- multiplicative coupling parameter Zp_cpl_R\n"
     <<std::setw(width+4)<<" "<<"}\n";
}
