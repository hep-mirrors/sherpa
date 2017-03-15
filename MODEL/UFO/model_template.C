#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "MODEL/UFO/UFO_Model.H"
#include "MODEL/Main/Single_Vertex.H"
#include "MODEL/UFO/UFO_Color_Functions.H"
#include <iomanip>

using namespace std;

namespace MODEL{

  class ${model_name} : public UFO::UFO_Model
  {
  public:
    ${model_name}(string path, string file, bool elementary) : UFO_Model(path, file, elementary) 
    {
      m_name = string("${model_name}");
      ParticleInit();
      ParamInit();
      SetMassiveFlags();
      SetStableFlags();
      SetSMMasses();
      // Massive and Stable flags
      // are set consistently with
      // UFO above, ReadParticleData
      // allows to overwrite these flags.
      // Need this before AddStandardContainers.
      ReadParticleData();
      AddStandardContainers();
      CustomContainerInit();
      ATOOLS::OutputParticles(ATOOLS::msg->Info());
      ATOOLS::OutputContainers(ATOOLS::msg->Info());
    }
  protected:
    void ParticleInit()
    {
      ATOOLS::s_kftable[kf_none] = new ATOOLS::Particle_Info(kf_none,-1,0,0,0,0,-1,0,1,0,"no_particle","no_particle","no_particle","no_particle",1,1);
      ${particle_init}
    }
    void ParamInit()
    {
      DEBUG_FUNC(this);
      msg_Debugging() << setprecision(20);
      ${param_init}
      msg_Debugging() << setprecision(6);
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

void ATOOLS::Getter<Model_Base,Model_Arguments,${model_name}>::PrintInfo(ostream &str,const size_t width) const 
{
  str<<"Automatically generated model \"${model_name}\" based on UFO output"<<endl;
}
