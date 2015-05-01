#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "MODEL/UFO/UFO_Model.H"
#include "MODEL/Main/Single_Vertex.H"
#include "MODEL/UFO/UFO_Color_Functions.H"
#include <iomanip>

namespace MODEL{

  class ${model_name} : public UFO::UFO_Model
  {
  public:
    ${model_name}(std::string path, std::string file, bool elementary) : UFO_Model(path, file, elementary) 
    {
      m_name = std::string("${model_name}");
      ParticleInit();
      ParamInit();
      CustomContainerInit();
      ATOOLS::OutputParticles(ATOOLS::msg->Info());
      ATOOLS::OutputContainers(ATOOLS::msg->Info());
    }
  protected:
    void ParticleInit()
    {
      ATOOLS::s_kftable[kf_none] = new ATOOLS::Particle_Info(kf_none,-1,0,0,0,0,-1,0,1,0,"no_particle","no_particle","no_particle","no_particle",1,1);
      ${particle_init}
      
      double jet_mass_threshold=MODEL::Model_Base::p_dataread->GetValue<double>("JET_MASS_THRESHOLD", 10.0);
      for (ATOOLS::KF_Table::iterator it=ATOOLS::s_kftable.begin(); it!=ATOOLS::s_kftable.end(); ++it) {
        if (it->second->m_mass<jet_mass_threshold) {
          it->second->m_massive=0;
          it->second->m_yuk=0;
        }
        if (it->second->m_width==0.) {
          it->second->m_stable=1;
        }
      }

    }
    void ParamInit()
    {
      DEBUG_FUNC(this);
      msg_Debugging() << std::setprecision(20);
      ${param_init}
      msg_Debugging() << std::setprecision(6);
      ReadParticleData();
      SetHadMass(kf_d,0.01);
      SetHadMass(kf_u,0.005);
      SetHadMass(kf_s,0.2);
      SetHadMass(kf_c,1.42);
      SetHadMass(kf_b,4.8);
      SetHadMass(kf_t,173.21);
      SetHadMass(kf_e,0.000511);
      SetHadMass(kf_mu,.105);
      SetHadMass(kf_tau,1.777);
      AddStandardContainers();
    }
    ${declarations}
    void InitVertices()
    {
      ${calls}
    }
  };
  
}

using namespace MODEL;

DECLARE_GETTER(${model_name},"${model_name}",Model_Base,Model_Arguments);

Model_Base *ATOOLS::Getter<Model_Base,Model_Arguments,${model_name}>::operator()(const Model_Arguments &args) const 
{
  return new ${model_name}(args.m_path, args.m_file, args.m_elementary);
}

void ATOOLS::Getter<Model_Base,Model_Arguments,${model_name}>::PrintInfo(std::ostream &str,const size_t width) const 
{
  str<<"Automatically generated model \"${model_name}\" based on UFO output"<<std::endl;
}
