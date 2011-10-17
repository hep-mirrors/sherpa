#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/MyComplex.H"
#include "MODEL/Main/Standard_Model.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "MODEL/Main/Effective_Higgs_Coupling.H"
#include "ATOOLS/Org/Message.H"
#include <iomanip>

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

namespace MODEL {
  class Standard_Model;
  class SM_TauPi : public Model_Base {
  private :
    Standard_Model *p_sm;
    void ParticleInit();
    void FillSpectrum(PDF::ISR_Handler *const isr);
  public :
    SM_TauPi(std::string,std::string,bool);
    ~SM_TauPi();
    bool ModelInit(PDF::ISR_Handler *const isr);
  };
}


DECLARE_GETTER(SM_TauPi_Getter,"SM+TauPi",Model_Base,Model_Arguments);

Model_Base *SM_TauPi_Getter::operator()(const Model_Arguments &args) const
{
  return new SM_TauPi(args.m_path,args.m_file,args.m_elementary);
}

void SM_TauPi_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The Standard Model + tau->pi nu vertex\n"
     <<std::setw(width+4)<<" "<<"{\n"
     <<std::setw(width+7)<<" "<<"parameter specification [keyword=value]\n"
     <<std::setw(width+7)<<" "<<"- all the SM parameters\n"
     <<std::setw(width+7)<<" "<<"- F_PI \n"
     <<std::setw(width+7)<<" "<<"- ...\n"
     <<std::setw(width+4)<<" "<<"}\n";
}


SM_TauPi::SM_TauPi(std::string _dir,std::string _file,bool _elementary) :
  Model_Base(_dir,_file,_elementary)
{
  p_sm = new Standard_Model(m_dir,m_file,false);
  ParticleInit();
  if (m_elementary) {
    ATOOLS::OutputParticles(msg->Info());
    ATOOLS::OutputContainers(msg->Info());
  }
}

bool SM_TauPi::ModelInit(PDF::ISR_Handler *const isr)
{
  if (m_elementary)
    msg_Info()<<"Initialize the Standard Model plus TauPi from "
	      <<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("SM+TauPi");

  p_sm->ModelInit(isr);
  p_numbers          = p_sm->ExtractScalarNumbers();
  p_constants        = p_sm->ExtractScalarConstants();
  p_complexconstants = p_sm->ExtractComplexConstants();
  p_functions        = p_sm->ExtractScalarFunctions();
  p_matrices         = p_sm->ExtractComplexMatrices();

  delete p_sm;

  FillSpectrum(isr);

  return true;
}

SM_TauPi::~SM_TauPi() 
{ }

void SM_TauPi::ParticleInit() {
  s_kftable[211] = new Particle_Info(211, 0.13957, 2.5242e-17,
                                     3, 0, 0, 1, 1, "pi+", "pi+");
  ReadParticleData();
}

void SM_TauPi::FillSpectrum(PDF::ISR_Handler *const isr) {
  p_dataread->RereadInFile();
  p_constants->insert(make_pair(string("F_PI"),
				p_dataread->GetValue<double>("F_PI",0.0924)));
}
