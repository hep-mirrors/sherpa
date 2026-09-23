#include "MODEL/SM/LowEnergy_Model.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "MODEL/Main/Single_Vertex.H"
#include "METOOLS/HadronCurrents/FormFactors/FormFactor_Base.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(MODEL::LowEnergy_Model,"LowEnergy",MODEL::Model_Base,MODEL::Model_Arguments);

Model_Base *ATOOLS::Getter<MODEL::Model_Base,MODEL::Model_Arguments,MODEL::LowEnergy_Model>::
operator()(const Model_Arguments &args) const
{
  return new LowEnergy_Model();
}

void ATOOLS::Getter<MODEL::Model_Base,MODEL::Model_Arguments,MODEL::LowEnergy_Model>::
PrintInfo(ostream &str,const size_t width) const
{
  str<<"The LowEnergy Model\n";
  str<<setw(width+4)<<" "<<"{\n"
     <<setw(width+7)<<" "<<"# possible parameters in yaml configuration [usage: \"keyword: value\"]\n"
     <<setw(width+7)<<" "<<"- 1/ALPHAQED(0) (alpha QED Thompson limit)\n"
     <<setw(width+7)<<" "<<"- SIN2THETAW (weak mixing angle)\n"
     <<setw(width+4)<<" "<<"}";
}

LowEnergy_Model::LowEnergy_Model() : Model_Base(true)
{
  msg_Out()<<METHOD<<" starts initialising.\n";
  m_name="LowEnergy";
  ParticleInit();
  RegisterDefaults();
}

bool LowEnergy_Model::ModelInit()
{
  Settings& s = Settings::GetMainSettings();
  m_alpha     = 1./s["1/ALPHAQED(0)"].Get<double>();
  m_sinthetaW = sqrt(s["SIN2THETAW"].Get<double>());
  msg_Out()<<METHOD<<": 1/alpha = "<<(1./m_alpha)<<"\n";
  return true;
}

void LowEnergy_Model::ParticleInit() {
  // kf_code, mass, radius, width,3*charge,2*spin,on,stable,idname,texname
  if (s_kftable.find(kf_p_plus)==s_kftable.end())
    s_kftable[kf_p_plus] =
      new Particle_Info(kf_p_plus,0.938272,0.8783,.0,3,1,1,1,"P+","P^{+}");
  if (s_kftable.find(kf_n)==s_kftable.end())
    s_kftable[kf_n] =
      new Particle_Info(kf_n,0.939566,0.8783,7.424e-28,0,1,1,1,"n","n");
}


void LowEnergy_Model::InitVertices() {
  InitQEDVertices();
  InitEWVertices();
}

void LowEnergy_Model::InitQEDVertices() {
  Kabbala g1("g_1",sqrt(4.*M_PI*m_alpha));
  Kabbala cpl=g1*Kabbala("i",Complex(0.,1.));
  Flavour flav;
  for (map<kf_code,Particle_Info*>::iterator kfit=s_kftable.begin();
       kfit!=s_kftable.end();kfit++) {
    Flavour flav = Flavour(kfit->first);
    // only create vertices for hadrons that are switched on.
    if (!flav.IsOn() || !flav.IsHadron()) continue;
    Kabbala Q("Q_{"+flav.TexName()+"}",flav.Charge());
    if (flav.IntSpin()==1 && flav.IsBaryon()) {
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle(flav.Bar());
      m_v.back().AddParticle(flav);
      m_v.back().AddParticle(Flavour(kf_photon));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().order[1]=1;
      m_v.back().cpl.push_back(cpl*Q);
      m_v.back().Lorentz.push_back("FFV");
    }
  }
}


void LowEnergy_Model::InitEWVertices() {}
  
